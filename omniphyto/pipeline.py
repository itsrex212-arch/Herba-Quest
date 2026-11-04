"""CSV-to-ligand-prep pipeline with optional, explicit AutoDock Vina stage."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from .docking import DockingError, VinaConfig, dock_ligand
from .ligand_prep import LigandPreparationError, dependency_status, prepare_smiles


def _safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    return cleaned or "unnamed"


def run_pipeline(
    input_csv: Path,
    output_dir: Path,
    *,
    id_column: str,
    smiles_column: str,
    seed: int,
    overwrite: bool,
    vina_config: VinaConfig | None = None,
) -> tuple[int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "manifest.csv"
    if manifest_path.exists() and not overwrite:
        raise LigandPreparationError(
            f"Refusing to overwrite existing manifest: {manifest_path}. Use --overwrite."
        )
    if vina_config is not None:
        vina_config.validate()
    with input_csv.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise LigandPreparationError("Input CSV has no header row")
        missing = [col for col in (id_column, smiles_column) if col not in reader.fieldnames]
        if missing:
            raise LigandPreparationError(
                "Input CSV is missing configured column(s): " + ", ".join(missing)
                + ". Set --id-column and --smiles-column to match your file."
            )
        rows = list(reader)

    records = []
    succeeded = 0
    for index, row in enumerate(rows, start=1):
        compound_id = (row.get(id_column) or "").strip() or f"row_{index}"
        smiles = row.get(smiles_column) or ""
        stem = _safe_name(compound_id)
        ligand_path = output_dir / f"{stem}.pdbqt"
        record = {
            "compound_id": compound_id,
            "smiles": smiles,
            "ligand_pdbqt": "",
            "ligand_status": "failed",
            "docking_status": "not_requested" if vina_config is None else "not_run",
            "score": "",
            "docking_output": "",
            "docking_log": "",
            "message": "",
        }
        try:
            prepare_smiles(smiles, ligand_path, seed=seed, overwrite=overwrite)
            succeeded += 1
            record["ligand_pdbqt"] = ligand_path.name
            record["ligand_status"] = "prepared"
            record["message"] = "Ligand PDBQT prepared; docking not requested"
            if vina_config is not None:
                docked_path = output_dir / f"{stem}_docked.pdbqt"
                log_path = output_dir / f"{stem}_vina.log"
                record["docking_output"] = docked_path.name
                record["docking_log"] = log_path.name
                try:
                    score, message = dock_ligand(
                        ligand_path, vina_config, docked_path, log_path, overwrite=overwrite
                    )
                    if score is not None:
                        record["score"] = score
                        record["docking_status"] = "completed"
                    else:
                        record["docking_status"] = "failed"
                    record["message"] = message
                except (DockingError, OSError) as exc:
                    record["docking_status"] = "failed"
                    record["message"] = f"Ligand prepared; docking failed: {exc}"
        except (LigandPreparationError, OSError) as exc:
            record["message"] = str(exc)
        records.append(record)

    fields = ["compound_id", "smiles", "ligand_pdbqt", "ligand_status",
              "docking_status", "score", "docking_output", "docking_log", "message"]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    return succeeded, len(rows)


def _triple(values: list[float], label: str) -> tuple[float, float, float]:
    if len(values) != 3:
        raise argparse.ArgumentTypeError(f"{label} requires exactly three values: x y z")
    return (values[0], values[1], values[2])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Prepare SMILES as PDBQT; optional real docking with AutoDock Vina."
    )
    parser.add_argument("input_csv", type=Path, help="CSV with configurable ID and SMILES columns")
    parser.add_argument("--output-dir", type=Path, default=Path("data/output"))
    parser.add_argument("--id-column", default="compound_id", help="provisional default: compound_id")
    parser.add_argument("--smiles-column", default="smiles", help="provisional default: smiles")
    parser.add_argument("--seed", type=int, default=2026, help="deterministic RDKit embedding seed")
    parser.add_argument("--overwrite", action="store_true", help="allow replacing existing outputs")
    parser.add_argument("--receptor", type=Path, help="receptor PDBQT; supplying this enables Vina")
    parser.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="Vina search-box center in angstroms")
    parser.add_argument("--size", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="Vina search-box dimensions in angstroms")
    parser.add_argument("--vina-executable", default="vina", help="Vina executable or full path")
    parser.add_argument("--vina-timeout", type=int, default=600, help="per-ligand timeout in seconds")
    args = parser.parse_args(argv)

    docking_args = (args.receptor, args.center, args.size)
    vina_config = None
    if any(value is not None for value in docking_args):
        if not all(value is not None for value in docking_args):
            parser.error("--receptor, --center, and --size must be supplied together")
        vina_config = VinaConfig(args.receptor, tuple(args.center), tuple(args.size),
                                 args.vina_executable, args.vina_timeout)
        try:
            vina_config.validate()
        except DockingError as exc:
            parser.error(str(exc))

    status = dependency_status()
    missing = [name for name, available in status.items() if not available]
    if missing:
        print("missing dependencies: " + ", ".join(missing), file=sys.stderr)
        print("see README.md for Windows 11 setup", file=sys.stderr)
        return 2
    try:
        ok, total = run_pipeline(args.input_csv, args.output_dir,
                                 id_column=args.id_column,
                                 smiles_column=args.smiles_column,
                                 seed=args.seed,
                                 overwrite=args.overwrite,
                                 vina_config=vina_config)
    except (LigandPreparationError, DockingError, OSError) as exc:
        print(f"pipeline error: {exc}", file=sys.stderr)
        return 2
    print(f"ligand prep: {ok}/{total} succeeded; manifest: {args.output_dir / 'manifest.csv'}")
    if vina_config is None:
        print("docking/scoring: NOT RUN (Vina was not configured; no scores claimed)")
    else:
        print(f"docking/scoring: requested with Vina; per-ligand statuses and logs are in {args.output_dir / 'manifest.csv'}")
    return 0 if ok == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
