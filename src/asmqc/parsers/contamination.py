# AsmQC — contamination report parser (NCBI FCS-GX and Kraken2).
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parse contamination screens.

**FCS-GX action report** (``*.fcs_gx_report.txt``) — verified against
``ncbi/fcs-gx`` ``scripts/action_report.py``:

* line 1: ``##[["FCS genome report", 2, 1], {...}]`` (JSON metadata)
* line 2: ``#seq_id\\tstart_pos\\tend_pos\\tseq_len\\taction\\tdiv\\tagg_cont_cov\\ttop_tax_name``
* data rows: 8 TAB columns; ``action`` ∈ EXCLUDE / TRIM / FIX / REVIEW / REVIEW_RARE / INFO.

**Kraken2 report** — verified against the Kraken2 manual: TAB, **no header**, 6
columns ``percent clade_reads taxon_reads rank_code taxid name`` (the name is
indented by tree depth; ``--report-minimizer-data`` inserts 2 extra columns
between col 3 and col 4 — handled).
"""
from __future__ import annotations

import re
from pathlib import Path

from asmqc.models import ContaminationHit, ContaminationResult
from asmqc.parsers.common import open_text, to_float, to_int

_RANK = re.compile(r"^[URDKPCOFGS](\d+)?$")
_FCS_ACTIONS = {"EXCLUDE", "TRIM", "FIX", "REVIEW", "REVIEW_RARE", "INFO"}


def parse_contamination(path: str | Path) -> ContaminationResult | None:
    """Auto-detect FCS-GX vs Kraken2 and dispatch."""
    p = Path(path)
    if not p.exists():
        return None
    kind = _sniff(p)
    if kind == "fcs":
        return parse_fcs(p)
    if kind == "kraken":
        return parse_kraken(p)
    # Last resort: try both.
    return parse_fcs(p) or parse_kraken(p)


def _sniff(path: Path) -> str | None:
    try:
        with open_text(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                if line.startswith("##") or line.startswith("#seq_id"):
                    return "fcs"
                parts = line.rstrip("\n").split("\t")
                # FCS data row: 8 cols with a known action in col 5.
                if len(parts) >= 8 and parts[4].strip() in _FCS_ACTIONS:
                    return "fcs"
                # Kraken: >=6 cols, col 4 is a rank code.
                if len(parts) >= 6 and _RANK.match(parts[3].strip()):
                    return "kraken"
                if len(parts) >= 8 and _RANK.match(parts[5].strip()):
                    return "kraken"  # minimizer-data variant
                break
    except OSError:
        return None
    return None


def parse_fcs(path: str | Path) -> ContaminationResult | None:
    result = ContaminationResult(source="fcs")
    seqs: set[str] = set()
    try:
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                seq_id = parts[0].strip()
                start = to_int(parts[1])
                end = to_int(parts[2])
                seq_len = to_int(parts[3])
                action = parts[4].strip() if len(parts) > 4 else None
                div = parts[5].strip() if len(parts) > 5 else None
                cov = to_float(parts[6]) if len(parts) > 6 else None
                taxon = parts[7].strip() if len(parts) > 7 else None
                if action not in _FCS_ACTIONS:
                    continue
                span = 0
                if start is not None and end is not None and end >= start:
                    span = end - start + 1
                result.hits.append(
                    ContaminationHit(
                        seq_id=seq_id,
                        action=action,
                        start=start,
                        end=end,
                        seq_len=seq_len,
                        taxon=taxon,
                        coverage=cov,
                        note=div,
                    )
                )
                result.total_flagged_bases += span
                seqs.add(seq_id)
    except OSError:
        return None
    if not result.hits:
        return None
    result.n_sequences_flagged = len(seqs)
    return result


def parse_kraken(path: str | Path) -> ContaminationResult | None:
    result = ContaminationResult(source="kraken")
    domains: dict[str, float] = {}
    saw_any = False
    try:
        with open_text(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                # Handle the 8-column --report-minimizer-data variant.
                if len(parts) >= 8 and _RANK.match(parts[5].strip()):
                    percent, _clade, _taxon = parts[0], parts[1], parts[2]
                    rank, taxid, name = parts[5], parts[6], parts[7]
                else:
                    percent, _clade, _taxon, rank, taxid, name = parts[:6]
                rank = rank.strip()
                if not _RANK.match(rank):
                    continue
                saw_any = True
                pct = to_float(percent)
                clean_name = name.strip()
                if rank == "U":
                    result.unclassified_pct = pct
                    continue
                # Superkingdom / domain-level summary feeds the contamination flag.
                if rank == "D" and pct is not None:
                    domains[clean_name] = pct
    except OSError:
        return None
    if not saw_any:
        return None
    result.top_taxa = domains
    # Build informational hits for non-target superkingdoms.
    for name, pct in sorted(domains.items(), key=lambda kv: -kv[1]):
        result.hits.append(
            ContaminationHit(seq_id="(read-level)", taxon=name, coverage=pct, note="superkingdom")
        )
    return result
