# AsmQC — QUAST report.tsv parser.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parse a QUAST ``report.tsv``.

Tab-separated, ``<metric name>\\t<value>`` (one extra value column per assembly in
multi-assembly runs — we take the first).  Metric-name strings verified verbatim
against QUAST source (``quast_libs/reporting.py``).  Accepts the ``report.tsv``
file directly or a QUAST output directory containing it.
"""
from __future__ import annotations

from pathlib import Path

from asmqc.models import AssemblyStats
from asmqc.parsers.common import open_text, to_float, to_int


def _locate(path: Path) -> Path | None:
    if path.is_dir():
        for candidate in ("report.tsv", "transposed_report.tsv"):
            p = path / candidate
            if p.exists():
                return p
        return None
    return path if path.exists() else None


def parse_quast(path: str | Path) -> AssemblyStats | None:
    report = _locate(Path(path))
    if report is None:
        return None
    try:
        rows: dict[str, str] = {}
        assembly_label: str | None = None
        with open_text(report) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                metric = parts[0].strip()
                value = parts[1].strip()
                if metric == "Assembly":
                    assembly_label = value
                    continue
                rows[metric] = value
    except OSError:
        return None

    if not rows:
        return None

    stats = AssemblyStats(source="quast")
    stats.total_length = to_int(rows.get("Total length")) or 0
    # QUAST '# contigs' counts sequences >= min_contig (default 500), NOT
    # gap-split contigs; '# contigs (>= 0 bp)' counts every sequence. QUAST does
    # not model scaffold gaps, so it cannot supply a true gap-split contig count
    # — we report num_contigs == num_sequences rather than the filtered count
    # (which would be < num_scaffolds and is not what num_contigs means).
    num_all = to_int(rows.get("# contigs (>= 0 bp)"))
    stats.num_sequences = num_all if num_all is not None else (to_int(rows.get("# contigs")) or 0)
    stats.num_scaffolds = stats.num_sequences
    stats.num_contigs = stats.num_sequences
    stats.n50 = to_int(rows.get("N50"))
    stats.n90 = to_int(rows.get("N90"))
    stats.l50 = to_int(rows.get("L50"))
    stats.l90 = to_int(rows.get("L90"))
    stats.auN = to_float(rows.get("auN"))
    stats.largest = to_int(rows.get("Largest contig"))
    stats.gc_percent = to_float(rows.get("GC (%)"))
    stats.n_per_100kbp = to_float(rows.get("# N's per 100 kbp"))
    if stats.total_length and stats.n_per_100kbp is not None:
        stats.gap_bases = int(round(stats.n_per_100kbp * stats.total_length / 100000.0))
    if assembly_label:
        stats.per_seq = []
    return stats
