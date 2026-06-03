# AsmQC — shared parser helpers.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Small utilities shared by the tool parsers."""
from __future__ import annotations

import gzip
from pathlib import Path


def open_text(path: str | Path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_text(path: str | Path) -> str | None:
    """Return the file's text, or ``None`` if it cannot be read."""
    try:
        with open_text(path) as fh:
            return fh.read()
    except (OSError, UnicodeError):
        return None


def to_float(value: object) -> float | None:
    """Best-effort float conversion: strips ``%``, commas, surrounding space."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip().rstrip("%").replace(",", "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def to_int(value: object) -> int | None:
    f = to_float(value)
    return int(round(f)) if f is not None else None


def first_existing(*paths: str | Path | None) -> Path | None:
    for p in paths:
        if p and Path(p).exists():
            return Path(p)
    return None
