"""Project settings persistence and validation, with no scientific defaults."""
from __future__ import annotations

import json
import math
from pathlib import Path

FIELDS = ("project_name", "target_name", "receptor_path", "box_center_x", "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z", "protocol")
DEFAULTS = {key: "" for key in FIELDS}


def load_settings(path: Path) -> dict[str, str]:
    if not Path(path).is_file():
        return dict(DEFAULTS)
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return dict(DEFAULTS)
    return {key: str(data.get(key, "")) for key in FIELDS}


def validate_settings(settings: dict[str, str], *, require_receptor: bool = False) -> list[str]:
    errors = []
    for key in ("project_name", "target_name", "protocol"):
        if not str(settings.get(key, "")).strip():
            errors.append(f"{key.replace('_', ' ').capitalize()} is required")
    box_keys = FIELDS[3:9]
    values = [str(settings.get(key, "")).strip() for key in box_keys]
    if any(values) and not all(values):
        errors.append("Box center and size require all six numeric coordinates")
    elif all(values):
        try:
            numbers = [float(value) for value in values]
            if not all(math.isfinite(n) for n in numbers):
                errors.append("Box coordinates must be finite numbers")
            elif any(value <= 0 for value in numbers[3:]):
                errors.append("Box sizes must be greater than zero")
        except ValueError:
            errors.append("Box coordinates must be numeric")
    if require_receptor:
        raw = str(settings.get("receptor_path", "")).strip()
        if not raw:
            errors.append("Receptor PDBQT path is required for docking")
        elif not Path(raw).expanduser().is_file():
            errors.append("Receptor file does not exist")
        if not all(values):
            errors.append("Complete target-specific box center and size are required for docking")
    return errors


def save_settings(path: Path, settings: dict[str, str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    clean = {key: str(settings.get(key, "")).strip() for key in FIELDS}
    Path(path).write_text(json.dumps(clean, indent=2), encoding="utf-8")
