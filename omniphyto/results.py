"""Read actual pipeline outputs, verifying prepared artifacts before results display."""
from __future__ import annotations

import csv
from pathlib import Path

RANKING_NAMES = ("ranking.csv", "ranked_results.csv", "rankings.csv")
POSE_SUFFIXES = {".pdbqt", ".pdb"}


def load_validated_preparations(output_dir: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Return verified prepared ligand rows and all manifest rows for errors/status review."""
    output_dir = Path(output_dir).resolve()
    manifest = output_dir / "manifest.csv"
    if not manifest.is_file():
        return [], []
    with manifest.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return [], []
        rows = [{str(key): (value or "") for key, value in row.items() if key is not None} for row in reader]
    valid = []
    for row in rows:
        artifact = (row.get("ligand_pdbqt") or "").strip()
        candidate = (output_dir / artifact).resolve() if artifact else None
        if (row.get("ligand_status", "").strip().lower() == "prepared"
                and candidate is not None and candidate.is_relative_to(output_dir)
                and candidate.is_file() and candidate.stat().st_size > 0):
            valid.append(row)
    return valid, rows


def load_results(output_dir: Path) -> tuple[list[dict[str, str]], str | None]:
    """Return only actual manifest rows; hide unverified scores and no demo fallback."""
    output_dir = Path(output_dir)
    ranking = next((output_dir / name for name in RANKING_NAMES if (output_dir / name).is_file()), None)
    path = ranking or (output_dir / "manifest.csv" if (output_dir / "manifest.csv").is_file() else None)
    if path is None:
        return [], None
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = [{str(key): (value or "") for key, value in row.items() if key is not None} for row in reader] if reader.fieldnames else []
    if path.name == "manifest.csv":
        for row in rows:
            if row.get("docking_status", "").strip().lower() != "completed":
                row["score"] = ""
    return rows, path.name


def available_pose_files(project_root: Path) -> list[Path]:
    """List existing docked output files only when manifest marks them completed, plus data/poses files."""
    root = Path(project_root).resolve()
    output = root / "data" / "output"
    candidates: set[Path] = set()
    manifest = output / "manifest.csv"
    completed_names = set()
    if manifest.is_file():
        with manifest.open("r", newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                if row.get("docking_status", "").strip().lower() == "completed":
                    name = (row.get("docking_output") or "").strip()
                    if name and Path(name).name == name and name.lower().endswith("_docked.pdbqt"):
                        completed_names.add(name)
    for path in output.iterdir() if output.is_dir() else ():
        if path.is_file() and path.name in completed_names:
            candidates.add(path.resolve())
    folder = root / "data" / "poses"
    if folder.is_dir():
        candidates.update(path.resolve() for path in folder.rglob("*") if path.is_file() and path.suffix.lower() in POSE_SUFFIXES and path.resolve().is_relative_to(root))
    return sorted(candidates, key=lambda p: str(p).lower())
