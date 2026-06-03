# AsmQC — data models shared across parsers, the flag engine and report writers.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Typed data containers that form the contract between every AsmQC module.

Every parser returns one of the ``*Result`` dataclasses (or ``None`` when its
input is absent).  The flag engine reads them, the report writers serialise
them.  Keeping the schema here means the JSON output, the HTML/MD report and the
tests all agree on field names.

All containers expose :meth:`to_dict` producing a JSON-serialisable structure
(plain dict/list/str/int/float/bool/None — no numpy scalars, no enums).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


# --------------------------------------------------------------------------- #
# Severity
# --------------------------------------------------------------------------- #
class Severity(Enum):
    """Ordered severity levels for curation flags.

    The integer value encodes rank so the *overall* report status is simply the
    maximum severity across all flags.
    """

    PASS = 0   # metric meets the standard
    INFO = 1   # neutral context, no action
    NOTE = 2   # minor observation worth a curator's glance
    WARN = 3   # likely issue; curation recommended
    FLAG = 4   # strong signal of a problem; curation needed
    FAIL = 5   # hard failure / blocking issue

    @property
    def label(self) -> str:
        return self.name

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name

    @classmethod
    def from_str(cls, value: str) -> Severity:
        return cls[value.strip().upper()]


# --------------------------------------------------------------------------- #
# Assembly contiguity (from FASTA or QUAST/gfastats)
# --------------------------------------------------------------------------- #
@dataclass
class SeqStat:
    """Per-sequence statistics for one FASTA record."""

    name: str
    length: int
    gc_percent: float | None = None
    n_bases: int = 0           # number of N/n bases
    gap_count: int = 0         # number of runs of >= min_gap N's (assembly gaps)
    gap_bases: int = 0         # total bases inside those gap runs
    num_subcontigs: int = 1    # contigs obtained after splitting on gaps

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "length": int(self.length),
            "gc_percent": _round(self.gc_percent),
            "n_bases": int(self.n_bases),
            "gap_count": int(self.gap_count),
            "gap_bases": int(self.gap_bases),
            "num_subcontigs": int(self.num_subcontigs),
        }


@dataclass
class AssemblyStats:
    """Contiguity / size statistics for a whole assembly.

    ``source`` records where the numbers came from: ``"fasta"`` (computed by
    AsmQC), ``"quast"`` or ``"gfastats"``.
    """

    source: str = "fasta"
    total_length: int = 0
    num_sequences: int = 0          # records in the FASTA (== scaffolds usually)
    num_scaffolds: int = 0
    num_contigs: int = 0            # scaffolds split on gaps
    n50: int | None = None
    n90: int | None = None
    l50: int | None = None
    l90: int | None = None
    contig_n50: int | None = None
    auN: float | None = None     # area under the Nx curve (Lh-insensitive contiguity)
    largest: int | None = None
    smallest: int | None = None
    mean_length: float | None = None
    median_length: int | None = None
    gc_percent: float | None = None
    gap_count: int = 0
    gap_bases: int = 0
    n_per_100kbp: float | None = None
    num_short_seqs: int = 0         # sequences below the "short scaffold" threshold
    short_seq_threshold: int = 0
    # Sorted (descending) scaffold lengths — used for the cumulative-length plot.
    sorted_lengths: list[int] = field(default_factory=list)
    # Optional per-sequence detail (omitted from JSON when large unless requested).
    per_seq: list[SeqStat] = field(default_factory=list)

    def to_dict(self, include_per_seq: bool = False) -> dict[str, Any]:
        d: dict[str, Any] = {
            "source": self.source,
            "total_length": int(self.total_length),
            "num_sequences": int(self.num_sequences),
            "num_scaffolds": int(self.num_scaffolds),
            "num_contigs": int(self.num_contigs),
            "n50": _int_or_none(self.n50),
            "n90": _int_or_none(self.n90),
            "l50": _int_or_none(self.l50),
            "l90": _int_or_none(self.l90),
            "contig_n50": _int_or_none(self.contig_n50),
            "auN": _round(self.auN),
            "largest": _int_or_none(self.largest),
            "smallest": _int_or_none(self.smallest),
            "mean_length": _round(self.mean_length),
            "median_length": _int_or_none(self.median_length),
            "gc_percent": _round(self.gc_percent),
            "gap_count": int(self.gap_count),
            "gap_bases": int(self.gap_bases),
            "n_per_100kbp": _round(self.n_per_100kbp),
            "num_short_seqs": int(self.num_short_seqs),
            "short_seq_threshold": int(self.short_seq_threshold),
        }
        if include_per_seq:
            d["per_seq"] = [s.to_dict() for s in self.per_seq]
        return d


# --------------------------------------------------------------------------- #
# BUSCO
# --------------------------------------------------------------------------- #
@dataclass
class BuscoResult:
    complete_pct: float | None = None
    single_pct: float | None = None
    duplicated_pct: float | None = None
    fragmented_pct: float | None = None
    missing_pct: float | None = None
    complete_n: int | None = None
    single_n: int | None = None
    duplicated_n: int | None = None
    fragmented_n: int | None = None
    missing_n: int | None = None
    total_n: int | None = None
    lineage: str | None = None
    mode: str | None = None
    busco_version: str | None = None
    summary_line: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "complete_pct": _round(self.complete_pct),
            "single_pct": _round(self.single_pct),
            "duplicated_pct": _round(self.duplicated_pct),
            "fragmented_pct": _round(self.fragmented_pct),
            "missing_pct": _round(self.missing_pct),
            "complete_n": _int_or_none(self.complete_n),
            "single_n": _int_or_none(self.single_n),
            "duplicated_n": _int_or_none(self.duplicated_n),
            "fragmented_n": _int_or_none(self.fragmented_n),
            "missing_n": _int_or_none(self.missing_n),
            "total_n": _int_or_none(self.total_n),
            "lineage": self.lineage,
            "mode": self.mode,
            "busco_version": self.busco_version,
            "summary_line": self.summary_line,
        }


# --------------------------------------------------------------------------- #
# Merqury (k-mer based QV / completeness / spectra)
# --------------------------------------------------------------------------- #
@dataclass
class SpectraPoint:
    copy_class: str        # "read-only", "0"/"absent", "1", "2", "3", "4", ">4"
    multiplicity: int      # k-mer multiplicity (x-axis)
    count: int             # number of k-mers (y-axis)


@dataclass
class MerquryResult:
    qv: float | None = None
    error_rate: float | None = None
    completeness_pct: float | None = None
    kmers_asm_only: int | None = None
    kmers_total: int | None = None
    kmer_size: int | None = None
    assembly_label: str | None = None
    # Raw per-assembly rows from a multi-haplotype .qv (label -> qv).
    per_assembly_qv: dict[str, float] = field(default_factory=dict)
    # spectra-cn histogram points (for the k-mer copy-number plot).
    spectra: list[SpectraPoint] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "qv": _round(self.qv, 4),
            "error_rate": self.error_rate,
            "completeness_pct": _round(self.completeness_pct),
            "kmers_asm_only": _int_or_none(self.kmers_asm_only),
            "kmers_total": _int_or_none(self.kmers_total),
            "kmer_size": _int_or_none(self.kmer_size),
            "assembly_label": self.assembly_label,
            "per_assembly_qv": {k: _round(v, 4) for k, v in self.per_assembly_qv.items()},
            "has_spectra": bool(self.spectra),
            "n_spectra_points": len(self.spectra),
        }


# --------------------------------------------------------------------------- #
# Telomeres (tidk)
# --------------------------------------------------------------------------- #
@dataclass
class TelomereScaffold:
    name: str
    length: int | None = None
    start_telomere: bool = False
    end_telomere: bool = False
    start_repeats: int = 0
    end_repeats: int = 0
    total_repeats: int = 0

    @property
    def both_ends(self) -> bool:
        return self.start_telomere and self.end_telomere

    @property
    def any_end(self) -> bool:
        return self.start_telomere or self.end_telomere

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "length": _int_or_none(self.length),
            "start_telomere": bool(self.start_telomere),
            "end_telomere": bool(self.end_telomere),
            "both_ends": self.both_ends,
            "start_repeats": int(self.start_repeats),
            "end_repeats": int(self.end_repeats),
            "total_repeats": int(self.total_repeats),
        }


@dataclass
class TelomereResult:
    repeat_motif: str | None = None
    window_size: int | None = None
    scaffolds: list[TelomereScaffold] = field(default_factory=list)
    # Only consider "large" scaffolds (putative chromosomes) for completeness math.
    chromosome_min_length: int = 0

    @property
    def n_scaffolds(self) -> int:
        return len(self.scaffolds)

    @property
    def n_with_any(self) -> int:
        return sum(1 for s in self.scaffolds if s.any_end)

    @property
    def n_with_both(self) -> int:
        return sum(1 for s in self.scaffolds if s.both_ends)

    def to_dict(self) -> dict[str, Any]:
        return {
            "repeat_motif": self.repeat_motif,
            "window_size": _int_or_none(self.window_size),
            "n_scaffolds": self.n_scaffolds,
            "n_with_any_telomere": self.n_with_any,
            "n_with_both_telomeres": self.n_with_both,
            "chromosome_min_length": int(self.chromosome_min_length),
            "scaffolds": [s.to_dict() for s in self.scaffolds],
        }


# --------------------------------------------------------------------------- #
# Contamination (FCS-GX / Kraken2)
# --------------------------------------------------------------------------- #
@dataclass
class ContaminationHit:
    seq_id: str
    action: str | None = None         # EXCLUDE / TRIM / REVIEW / FIX (FCS) or None
    start: int | None = None
    end: int | None = None
    seq_len: int | None = None
    taxon: str | None = None
    tax_id: str | None = None
    coverage: float | None = None     # fraction / count depending on source
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "seq_id": self.seq_id,
            "action": self.action,
            "start": _int_or_none(self.start),
            "end": _int_or_none(self.end),
            "seq_len": _int_or_none(self.seq_len),
            "taxon": self.taxon,
            "tax_id": self.tax_id,
            "coverage": self.coverage,
            "note": self.note,
        }


@dataclass
class ContaminationResult:
    source: str = "unknown"              # "fcs" | "kraken"
    hits: list[ContaminationHit] = field(default_factory=list)
    total_flagged_bases: int = 0
    n_sequences_flagged: int = 0
    # For kraken: top non-host taxa percentages (name -> percent).
    top_taxa: dict[str, float] = field(default_factory=dict)
    unclassified_pct: float | None = None

    @property
    def n_hits(self) -> int:
        return len(self.hits)

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "n_hits": self.n_hits,
            "n_sequences_flagged": int(self.n_sequences_flagged),
            "total_flagged_bases": int(self.total_flagged_bases),
            "unclassified_pct": _round(self.unclassified_pct),
            "top_taxa": {k: _round(v) for k, v in self.top_taxa.items()},
            "hits": [h.to_dict() for h in self.hits],
        }


# --------------------------------------------------------------------------- #
# Curation flag
# --------------------------------------------------------------------------- #
@dataclass
class Flag:
    """One curation flag produced by the heuristics engine."""

    id: str
    severity: Severity
    title: str
    category: str                       # "contiguity" | "completeness" | ...
    message: str                        # what was observed (with numbers)
    explanation: str = ""               # plain-English likely cause
    action: str = ""                    # recommended curation action
    evidence: dict[str, Any] = field(default_factory=dict)
    metric: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "severity": self.severity.name,
            "severity_rank": self.severity.value,
            "title": self.title,
            "category": self.category,
            "message": self.message,
            "explanation": self.explanation,
            "action": self.action,
            "metric": self.metric,
            "evidence": _jsonable(self.evidence),
        }


# --------------------------------------------------------------------------- #
# Top-level report
# --------------------------------------------------------------------------- #
@dataclass
class QCReport:
    assembly_name: str
    stats: AssemblyStats | None = None
    busco: BuscoResult | None = None
    merqury: MerquryResult | None = None
    telomere: TelomereResult | None = None
    contamination: ContaminationResult | None = None
    flags: list[Flag] = field(default_factory=list)
    inputs: dict[str, Any] = field(default_factory=dict)     # which files were used
    config: dict[str, Any] = field(default_factory=dict)     # thresholds applied
    asmqc_version: str = ""
    generated_at: str | None = None                       # ISO timestamp (caller-set)
    warnings: list[str] = field(default_factory=list)        # non-fatal parse issues

    @property
    def overall_severity(self) -> Severity:
        if not self.flags:
            return Severity.PASS
        return max((f.severity for f in self.flags), key=lambda s: s.value)

    def counts_by_severity(self) -> dict[str, int]:
        counts = {s.name: 0 for s in Severity}
        for f in self.flags:
            counts[f.severity.name] += 1
        return counts

    def to_dict(self, include_per_seq: bool = False) -> dict[str, Any]:
        return {
            "asmqc_version": self.asmqc_version,
            "generated_at": self.generated_at,
            "assembly_name": self.assembly_name,
            "overall_status": self.overall_severity.name,
            "flag_counts": self.counts_by_severity(),
            "inputs": _jsonable(self.inputs),
            "config": _jsonable(self.config),
            "warnings": list(self.warnings),
            "assembly_stats": self.stats.to_dict(include_per_seq) if self.stats else None,
            "busco": self.busco.to_dict() if self.busco else None,
            "merqury": self.merqury.to_dict() if self.merqury else None,
            "telomere": self.telomere.to_dict() if self.telomere else None,
            "contamination": self.contamination.to_dict() if self.contamination else None,
            "flags": [f.to_dict() for f in self.flags],
        }


# --------------------------------------------------------------------------- #
# Serialisation helpers
# --------------------------------------------------------------------------- #
def _round(value: float | None, ndigits: int = 2) -> float | None:
    if value is None:
        return None
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return None


def _int_or_none(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _jsonable(obj: Any) -> Any:
    """Recursively coerce numpy scalars / enums / sets into JSON-safe types."""
    if obj is None or isinstance(obj, (str, bool)):
        return obj
    if isinstance(obj, Severity):
        return obj.name
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    # numpy / other numeric scalar
    if hasattr(obj, "item") and not isinstance(obj, type):
        try:
            return obj.item()
        except (ValueError, TypeError):
            pass
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return float(obj)
    return obj
