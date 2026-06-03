# AsmQC — gfastats summary parser.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parse gfastats default summary output.

Default mode prints ``<label>: <value>`` (colon-space).  ``--tabular`` mode uses a
single TAB and omits the ``+++Assembly summary+++`` header.  Labels/order verified
against gfastats source (``vgl-hub/gfalibs`` ``src/output.cpp``).

Caveat handled here: the ``Base composition (A:C:G:T)`` label *and* value both
contain colons, so we split on the **last** ``: `` (or the TAB in tabular mode).
"""
from __future__ import annotations

from pathlib import Path

from asmqc.models import AssemblyStats
from asmqc.parsers.common import open_text, to_float, to_int

_LABELS = {
    "# scaffolds": "num_scaffolds",
    "Total scaffold length": "total_length",
    "Average scaffold length": "mean_length",
    "Scaffold N50": "n50",
    "Scaffold auN": "auN",
    "Scaffold L50": "l50",
    "Largest scaffold": "largest",
    "Smallest scaffold": "smallest",
    "# contigs": "num_contigs",
    "Contig N50": "contig_n50",
    "# gaps in scaffolds": "gap_count",
    "Total gap length in scaffolds": "gap_bases",
    "GC content %": "gc_percent",
}
_INT_FIELDS = {
    "num_scaffolds",
    "total_length",
    "n50",
    "l50",
    "largest",
    "smallest",
    "num_contigs",
    "contig_n50",
    "gap_count",
    "gap_bases",
}


def _split_label_value(line: str) -> tuple[str, str] | None:
    line = line.rstrip("\n")
    if "\t" in line:  # tabular mode
        label, _, value = line.partition("\t")
        return label.strip(), value.strip()
    if ": " in line:
        label, _, value = line.rpartition(": ")
        return label.strip(), value.strip()
    return None


def parse_gfastats(path: str | Path) -> AssemblyStats | None:
    p = Path(path)
    if not p.exists():
        return None
    parsed: dict[str, str] = {}
    try:
        with open_text(p) as fh:
            for line in fh:
                if line.startswith("+++") or not line.strip():
                    continue
                kv = _split_label_value(line)
                if kv is None:
                    continue
                label, value = kv
                if label in _LABELS:
                    parsed[_LABELS[label]] = value
    except OSError:
        return None

    if not parsed:
        return None

    stats = AssemblyStats(source="gfastats")
    for field, raw in parsed.items():
        if field in _INT_FIELDS:
            setattr(stats, field, to_int(raw) or (0 if field in ("gap_count", "gap_bases") else None))
        else:
            setattr(stats, field, to_float(raw))
    stats.num_sequences = stats.num_scaffolds or 0
    if stats.total_length and stats.gap_bases:
        stats.n_per_100kbp = 100000.0 * stats.gap_bases / stats.total_length
    return stats
