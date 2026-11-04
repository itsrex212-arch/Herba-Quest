"""CSV intake validation for real compound inputs."""
from __future__ import annotations

import csv
import io

MAX_UPLOAD_BYTES = 5 * 1024 * 1024
PREVIEW_ROWS = 20
REQUIRED_COLUMNS = ("compound_id", "smiles")


def validate_compound_csv(content: bytes, *, id_column: str = "compound_id", smiles_column: str = "smiles", max_bytes: int = MAX_UPLOAD_BYTES) -> tuple[list[str], list[dict[str, str]]]:
    """Validate full CSV structure, required configured columns, unique IDs and nonempty SMILES.

    No chemistry interpretation occurs here. Returns all validated rows for pipeline use.
    """
    if not content:
        raise ValueError("The uploaded CSV is empty")
    if len(content) > max_bytes:
        raise ValueError(f"CSV exceeds the {max_bytes // (1024 * 1024)} MB upload limit")
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("CSV must be encoded as UTF-8") from exc
    try:
        reader = csv.DictReader(io.StringIO(text, newline=""))
        headers = reader.fieldnames
        if not headers:
            raise ValueError("CSV must include a header row")
        if any(not name.strip() for name in headers):
            raise ValueError("CSV column names cannot be empty")
        if len(set(headers)) != len(headers):
            raise ValueError("CSV column names must be unique")
        missing = [name for name in (id_column, smiles_column) if name not in headers]
        if missing:
            raise ValueError("Missing required column(s): " + ", ".join(missing))
        rows = []
        ids = set()
        for line, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Row {line} has more fields than the header")
            clean = {key: (value or "").strip() for key, value in row.items() if key is not None}
            compound_id, smiles = clean[id_column], clean[smiles_column]
            if not compound_id:
                raise ValueError(f"Row {line}: {id_column} is required")
            if not smiles:
                raise ValueError(f"Row {line}: {smiles_column} is required")
            if compound_id in ids:
                raise ValueError(f"Row {line}: duplicate {id_column} '{compound_id}'")
            ids.add(compound_id)
            rows.append(clean)
        if not rows:
            raise ValueError("CSV has headers but no compound rows")
        return list(headers), rows
    except csv.Error as exc:
        raise ValueError(f"Could not parse CSV: {exc}") from exc


def preview_csv(content: bytes, *, max_bytes: int = MAX_UPLOAD_BYTES) -> tuple[list[str], list[dict[str, str]]]:
    """Legacy structural preview helper. Does not require a compound schema."""
    if not content:
        raise ValueError("The uploaded CSV is empty")
    if len(content) > max_bytes:
        raise ValueError(f"CSV exceeds the {max_bytes // (1024 * 1024)} MB upload limit")
    try:
        text = content.decode("utf-8-sig")
        reader = csv.DictReader(io.StringIO(text, newline=""))
        headers = reader.fieldnames
        if not headers:
            raise ValueError("CSV must include a header row")
        if any(not name.strip() for name in headers) or len(set(headers)) != len(headers):
            raise ValueError("CSV column names must be nonempty and unique")
        rows = []
        for line, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Row {line} has more fields than the header")
            rows.append({key: (value or "") for key, value in row.items() if key is not None})
        return list(headers), rows[:PREVIEW_ROWS]
    except UnicodeDecodeError as exc:
        raise ValueError("CSV must be encoded as UTF-8") from exc
    except csv.Error as exc:
        raise ValueError(f"Could not parse CSV: {exc}") from exc
