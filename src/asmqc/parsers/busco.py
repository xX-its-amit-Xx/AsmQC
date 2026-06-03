# AsmQC — BUSCO short_summary parser.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parse a BUSCO ``short_summary`` file (``.txt`` or ``.json``).

Format verified against the BUSCO source (gitlab.com/ezlab/busco) and real
example files spanning BUSCO 5.6 - 6.0:

* The canonical one-line notation::

      C:93.3%[S:92.0%,D:1.3%],F:3.1%,M:3.6%,n:954

* Six labelled count lines, e.g. ``890\\tComplete BUSCOs (C)``.
* The JSON variant exposes ``results`` with keys ``Complete percentage``,
  ``Single copy BUSCOs``, ``Multi copy percentage`` (note: a space, and
  "Multi copy" rather than "Duplicated"), ``n_markers`` and ``domain``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from asmqc.models import BuscoResult
from asmqc.parsers.common import read_text, to_float, to_int

# One-line BUSCO notation.  Order C,[S,D],F,M,n is fixed.
_ONE_LINE = re.compile(
    r"C:\s*([\d.]+)%\s*\[\s*S:\s*([\d.]+)%\s*,\s*D:\s*([\d.]+)%\s*\]\s*,\s*"
    r"F:\s*([\d.]+)%\s*,\s*M:\s*([\d.]+)%\s*,\s*n:\s*(\d+)"
)
_COUNT_LINE = re.compile(
    r"^\s*(\d+)\s+(Complete BUSCOs|Complete and single-copy BUSCOs|"
    r"Complete and duplicated BUSCOs|Fragmented BUSCOs|Missing BUSCOs|"
    r"Total BUSCO groups searched)\b"
)
_VERSION = re.compile(r"BUSCO version is:\s*([\w.]+)")
_LINEAGE = re.compile(r"lineage dataset is:\s*(\S+)")
_MODE = re.compile(r"run in mode:\s*(\w+)")


def parse_busco(path: str | Path) -> BuscoResult | None:
    text = read_text(path)
    if text is None:
        return None
    stripped = text.lstrip()
    # JSON files start with '{'; the official schema nests a 'results' object.
    if stripped.startswith("{"):
        result = _parse_json(text)
        if result is not None:
            return result
        # fall through: maybe a .json that is actually the txt notation
    return _parse_txt(text)


def _parse_json(text: str) -> BuscoResult | None:
    try:
        data = json.loads(text)
    except (ValueError, TypeError):
        return None
    if not isinstance(data, dict):
        return None
    results = data.get("results")
    if not isinstance(results, dict):
        return None

    r = BuscoResult()
    r.complete_pct = to_float(results.get("Complete percentage"))
    r.single_pct = to_float(results.get("Single copy percentage"))
    r.duplicated_pct = to_float(results.get("Multi copy percentage"))
    r.fragmented_pct = to_float(results.get("Fragmented percentage"))
    r.missing_pct = to_float(results.get("Missing percentage"))
    r.complete_n = to_int(results.get("Complete BUSCOs"))
    r.single_n = to_int(results.get("Single copy BUSCOs"))
    r.duplicated_n = to_int(results.get("Multi copy BUSCOs"))
    r.fragmented_n = to_int(results.get("Fragmented BUSCOs"))
    r.missing_n = to_int(results.get("Missing BUSCOs"))
    r.total_n = to_int(results.get("n_markers"))
    r.summary_line = results.get("one_line_summary")

    lineage = data.get("lineage_dataset")
    if isinstance(lineage, dict):
        r.lineage = lineage.get("name")
    params = data.get("parameters")
    if isinstance(params, dict):
        r.mode = params.get("mode") or r.mode
        if not r.lineage and isinstance(params.get("lineage_dataset"), str):
            r.lineage = Path(params["lineage_dataset"]).name
    versions = data.get("versions")
    if isinstance(versions, dict):
        r.busco_version = _stringify(versions.get("busco"))

    # If percentages were absent, derive them from counts.
    _fill_percentages(r)
    return r if _has_any(r) else None


def _parse_txt(text: str) -> BuscoResult | None:
    r = BuscoResult()
    m = _ONE_LINE.search(text)
    if m:
        r.complete_pct = to_float(m.group(1))
        r.single_pct = to_float(m.group(2))
        r.duplicated_pct = to_float(m.group(3))
        r.fragmented_pct = to_float(m.group(4))
        r.missing_pct = to_float(m.group(5))
        r.total_n = to_int(m.group(6))
        r.summary_line = m.group(0)

    label_to_field = {
        "Complete BUSCOs": "complete_n",
        "Complete and single-copy BUSCOs": "single_n",
        "Complete and duplicated BUSCOs": "duplicated_n",
        "Fragmented BUSCOs": "fragmented_n",
        "Missing BUSCOs": "missing_n",
        "Total BUSCO groups searched": "total_n",
    }
    for line in text.splitlines():
        cm = _COUNT_LINE.match(line)
        if cm:
            setattr(r, label_to_field[cm.group(2)], to_int(cm.group(1)))

    vm = _VERSION.search(text)
    if vm:
        r.busco_version = vm.group(1)
    lm = _LINEAGE.search(text)
    if lm:
        r.lineage = lm.group(1)
    mm = _MODE.search(text)
    if mm:
        r.mode = mm.group(1)

    _fill_percentages(r)
    return r if _has_any(r) else None


def _fill_percentages(r: BuscoResult) -> None:
    """Derive percentages from counts when only counts are present."""
    if not r.total_n:
        return
    total = r.total_n

    def pct(n: int | None) -> float | None:
        return round(100.0 * n / total, 4) if n is not None and total else None

    if r.complete_pct is None and r.complete_n is not None:
        r.complete_pct = pct(r.complete_n)
    if r.single_pct is None and r.single_n is not None:
        r.single_pct = pct(r.single_n)
    if r.duplicated_pct is None and r.duplicated_n is not None:
        r.duplicated_pct = pct(r.duplicated_n)
    if r.fragmented_pct is None and r.fragmented_n is not None:
        r.fragmented_pct = pct(r.fragmented_n)
    if r.missing_pct is None and r.missing_n is not None:
        r.missing_pct = pct(r.missing_n)
    # Complete = single + duplicated if we only have the parts.
    if r.complete_n is None and r.single_n is not None and r.duplicated_n is not None:
        r.complete_n = r.single_n + r.duplicated_n
        if r.complete_pct is None:
            r.complete_pct = pct(r.complete_n)


def _has_any(r: BuscoResult) -> bool:
    return any(
        v is not None
        for v in (r.complete_pct, r.complete_n, r.single_n, r.summary_line)
    )


def _stringify(value: object) -> str | None:
    return None if value is None else str(value)
