"""Dependency-light helpers for displaying real OmniPhyto outputs."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read stored CSV text rows without inventing or transforming results."""
    with Path(path).open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        return [{str(key): (value or "") for key, value in row.items() if key is not None} for row in reader]


def load_dashboard_data(project_root: Path) -> tuple[list[dict[str, str]], str, bool]:
    """Load only an actual pipeline manifest. There is deliberately no sample fallback."""
    path = Path(project_root) / "data" / "output" / "manifest.csv"
    if not path.is_file():
        return [], "", False
    return read_csv_rows(path), path.name, False


def guess_compound_column(rows: list[dict[str, str]]) -> str | None:
    if not rows:
        return None
    fields = list(rows[0])
    return next((name for name in ("compound_id", "compound_name", "name", "id") if name in fields), fields[0] if fields else None)


def filter_rows(rows: Iterable[dict[str, str]], *, query: str = "", field: str | None = None, value: str | None = None) -> list[dict[str, str]]:
    query = query.strip().casefold()
    kept = []
    for row in rows:
        if query and not any(query in str(item).casefold() for item in row.values()):
            continue
        if field and value is not None and str(row.get(field, "")) != value:
            continue
        kept.append(dict(row))
    return kept


def sort_rows(rows: Iterable[dict[str, str]], field: str | None, descending: bool = False) -> list[dict[str, str]]:
    result = [dict(row) for row in rows]
    if not field:
        return result
    populated = [row for row in result if str(row.get(field, "")).strip()]
    blank = [row for row in result if not str(row.get(field, "")).strip()]
    try:
        keyed = [(float(str(row[field]).strip()), row) for row in populated]
        populated = [row for _, row in sorted(keyed, key=lambda pair: pair[0], reverse=descending)]
    except (ValueError, TypeError):
        populated.sort(key=lambda row: str(row.get(field, "")).casefold(), reverse=descending)
    return populated + blank
