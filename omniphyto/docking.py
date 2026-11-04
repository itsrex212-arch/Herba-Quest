"""Safe, opt-in AutoDock Vina execution for prepared ligand PDBQT files."""
from __future__ import annotations

import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class DockingError(RuntimeError):
    """Raised for invalid Vina settings or inputs before a docking run."""


@dataclass(frozen=True)
class VinaConfig:
    receptor: Path
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    executable: str = "vina"
    timeout_seconds: int = 600

    def validate(self) -> str:
        receptor = Path(self.receptor)
        if not receptor.is_file():
            raise DockingError(f"Receptor PDBQT does not exist or is not a file: {receptor}")
        if len(self.center) != 3 or len(self.size) != 3:
            raise DockingError("Vina box center and size each require x, y, and z values")
        if not all(math.isfinite(float(value)) for value in self.center):
            raise DockingError("Vina box center values must be finite numbers")
        if not all(math.isfinite(float(value)) and float(value) > 0 for value in self.size):
            raise DockingError("Vina box size values must be finite numbers greater than zero")
        if self.timeout_seconds <= 0:
            raise DockingError("Vina timeout must be greater than zero")
        resolved = shutil.which(self.executable)
        if resolved is None:
            candidate = Path(self.executable)
            if candidate.is_file():
                resolved = str(candidate.resolve())
        if resolved is None:
            raise DockingError(f"AutoDock Vina executable was not found: {self.executable}")
        return resolved


def _read_score(text: str) -> str | None:
    """Extract the best-mode score from Vina's output remarks, if present."""
    match = re.search(r"^REMARK VINA RESULT:\s*(-?(?:\d+(?:\.\d*)?|\.\d+))", text, re.MULTILINE)
    return match.group(1) if match else None


def dock_ligand(
    ligand: Path,
    config: VinaConfig,
    output_pdbqt: Path,
    log_path: Path,
    *,
    overwrite: bool = False,
) -> tuple[str | None, str]:
    """Run Vina without a shell; return (best score, diagnostic message).

    Scores are read only from Vina's actual output PDBQT REMARK. A successful
    process without that record is not reported as a scored docking result.
    """
    executable = config.validate()
    ligand = Path(ligand)
    if not ligand.is_file():
        raise DockingError(f"Ligand PDBQT does not exist: {ligand}")
    for artifact in (Path(output_pdbqt), Path(log_path)):
        if artifact.exists() and not overwrite:
            raise DockingError(f"Refusing to overwrite existing docking output: {artifact}")
        artifact.parent.mkdir(parents=True, exist_ok=True)
    out_path = Path(output_pdbqt)
    if out_path.exists() and overwrite:
        # Avoid treating a previous run's output/score as this run's result.
        out_path.unlink()

    command = [
        executable,
        "--receptor", str(Path(config.receptor).resolve()),
        "--ligand", str(ligand.resolve()),
        "--center_x", str(float(config.center[0])),
        "--center_y", str(float(config.center[1])),
        "--center_z", str(float(config.center[2])),
        "--size_x", str(float(config.size[0])),
        "--size_y", str(float(config.size[1])),
        "--size_z", str(float(config.size[2])),
        "--out", str(Path(output_pdbqt).resolve()),
        "--log", str(Path(log_path).resolve()),
    ]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False,
                                timeout=config.timeout_seconds, shell=False)
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        Path(log_path).write_text(
            f"command: {command!r}\nstatus: timeout after {config.timeout_seconds} seconds\n"
            f"stdout:\n{stdout}\nstderr:\n{stderr}\n", encoding="utf-8")
        return None, f"Vina timed out after {config.timeout_seconds} seconds"
    except OSError as exc:
        Path(log_path).write_text(f"command: {command!r}\nexecution error: {exc}\n", encoding="utf-8")
        return None, f"Could not execute Vina: {exc}"

    Path(log_path).write_text(
        f"command: {command!r}\nexit_code: {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}\n", encoding="utf-8")
    if result.returncode != 0:
        return None, f"Vina exited with code {result.returncode}; see {Path(log_path).name}"
    out_path = Path(output_pdbqt)
    if not out_path.is_file():
        return None, "Vina returned success but did not create its output PDBQT"
    score = _read_score(out_path.read_text(encoding="utf-8", errors="replace"))
    if score is None:
        return None, "Vina output contained no REMARK VINA RESULT score"
    return score, "Vina docking completed; score read from output PDBQT"
