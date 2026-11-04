"""Convert a SMILES string into a 3D PDBQT ligand.

RDKit builds and minimizes a 3D conformer, Meeko writes AutoDock-compatible
PDBQT atom types/charges, and the Open Babel CLI performs a parse/round-trip
sanity check. This module does not dock ligands or calculate binding scores.
"""
from __future__ import annotations

import importlib
import shutil
import subprocess
import tempfile
from pathlib import Path


class LigandPreparationError(RuntimeError):
    """Raised when a molecule cannot be prepared or validated."""


def dependency_status() -> dict[str, bool]:
    """Return availability of required Python packages and Open Babel CLI."""
    status = {}
    for package in ("rdkit", "meeko"):
        try:
            importlib.import_module(package)
            status[package] = True
        except ImportError:
            status[package] = False
    status["obabel"] = shutil.which("obabel") is not None
    return status


def require_dependencies() -> None:
    missing = [name for name, available in dependency_status().items() if not available]
    if missing:
        raise LigandPreparationError(
            "Missing required ligand-prep dependencies: " + ", ".join(missing)
            + ". See README.md for Windows 11 setup instructions."
        )


def _validate_with_openbabel(pdbqt_path: Path) -> None:
    """Check that Open Babel can parse the generated PDBQT and emit PDB."""
    with tempfile.TemporaryDirectory(prefix="omniphyto_obabel_") as temp_dir:
        check_path = Path(temp_dir) / "roundtrip.pdb"
        result = subprocess.run(
            ["obabel", "-ipdbqt", str(pdbqt_path), "-opdb", "-O", str(check_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0 or not check_path.exists() or check_path.stat().st_size == 0:
            detail = (result.stderr or result.stdout or "no diagnostic returned").strip()
            raise LigandPreparationError(f"Open Babel could not validate PDBQT: {detail}")


def prepare_smiles(
    smiles: str,
    output_path: str | Path,
    *,
    seed: int = 2026,
    overwrite: bool = False,
    validate_openbabel: bool = True,
) -> Path:
    """Build one 3D ligand and write PDBQT; returns the output path.

    The output is only committed after RDKit, Meeko, and (by default) Open
    Babel validation succeed. Existing outputs are protected unless overwrite
    is explicitly enabled.
    """
    require_dependencies()
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    smiles = (smiles or "").strip()
    if not smiles:
        raise LigandPreparationError("SMILES is empty")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise LigandPreparationError("RDKit could not parse the supplied SMILES")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    embedded = AllChem.EmbedMolecule(mol, params)
    if embedded != 0:
        raise LigandPreparationError("RDKit could not generate a 3D conformer")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol)
    else:
        raise LigandPreparationError("RDKit has no MMFF or UFF parameters for this molecule")

    try:
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            raise LigandPreparationError("Meeko returned no prepared molecule")
        pdbqt_text, success, error = PDBQTWriterLegacy.write_string(setups[0])
        if not success or not pdbqt_text:
            raise LigandPreparationError(f"Meeko could not write PDBQT: {error}")
    except LigandPreparationError:
        raise
    except Exception as exc:
        raise LigandPreparationError(f"Meeko preparation failed: {exc}") from exc

    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not overwrite:
        raise LigandPreparationError(f"Refusing to overwrite existing file: {target}")
    # Write to a neighboring temporary file so failed validation cannot leave
    # a partial final artifact behind.
    temp_path = target.with_name(target.name + ".tmp.pdbqt")
    try:
        temp_path.write_text(pdbqt_text, encoding="utf-8")
        if validate_openbabel:
            _validate_with_openbabel(temp_path)
        if target.exists() and not overwrite:
            raise LigandPreparationError(f"Refusing to overwrite existing file: {target}")
        temp_path.replace(target)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return target
