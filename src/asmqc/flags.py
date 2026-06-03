# AsmQC — curation flag engine.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Heuristics that turn parsed QC metrics into actionable curation flags.

The engine reads the populated result models on a :class:`~asmqc.models.QCReport`
and the thresholds in :class:`~asmqc.config.AsmQCConfig`, and returns a list of
:class:`~asmqc.models.Flag`.  Every flag carries a plain-English cause and a
recommended action (from :mod:`asmqc.explanations`).
"""
from __future__ import annotations

import re

from asmqc.config import AsmQCConfig
from asmqc.explanations import get as explain
from asmqc.models import (
    AssemblyStats,
    BuscoResult,
    ContaminationResult,
    Flag,
    MerquryResult,
    QCReport,
    Severity,
    TelomereResult,
)


def evaluate(report: QCReport, config: AsmQCConfig) -> list[Flag]:
    """Run every enabled heuristic family and collect the flags."""
    flags: list[Flag] = []
    if report.busco and config.enabled("busco"):
        flags.extend(_busco_flags(report.busco, config))
    if report.merqury and config.enabled("merqury"):
        flags.extend(_merqury_flags(report.merqury, config))
    if report.stats and config.enabled("contiguity"):
        flags.extend(_contiguity_flags(report.stats, config))
    if report.stats and config.enabled("gaps"):
        flags.extend(_gap_flags(report.stats, config))
    if report.telomere and config.enabled("telomere"):
        flags.extend(_telomere_flags(report.telomere, config))
    if report.contamination and config.enabled("contamination"):
        flags.extend(_contamination_flags(report.contamination, config))
    # Stable sort: most severe first, then by category.
    flags.sort(key=lambda f: (-f.severity.value, f.category, f.id))
    return flags


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _flag(id_, severity, title, category, message, key, metric=None, evidence=None) -> Flag:
    explanation, action = explain(key)
    if severity is Severity.PASS:
        explanation, action = explain("pass")
    return Flag(
        id=id_,
        severity=severity,
        title=title,
        category=category,
        message=message,
        explanation=explanation,
        action=action,
        evidence=evidence or {},
        metric=metric,
    )


def _fmt_bp(value: float | None) -> str:
    if value is None:
        return "n/a"
    value = float(value)
    if value >= 1e9:
        return f"{value / 1e9:.2f} Gb"
    if value >= 1e6:
        return f"{value / 1e6:.2f} Mb"
    if value >= 1e3:
        return f"{value / 1e3:.1f} kb"
    return f"{int(value)} bp"


def _span(hit) -> int:
    """1-based inclusive span of a contamination hit, or 0 if unknown."""
    if hit.start is None or hit.end is None or hit.end < hit.start:
        return 0
    return hit.end - hit.start + 1


def _clean_taxon(text: str | None) -> str:
    """Neutralise control chars and markup/table-breaking characters in a
    free-text taxon name before embedding it in a flag message."""
    if not text:
        return ""
    return re.sub(r"[<>|\x00-\x1f]", " ", str(text)).strip()


# --------------------------------------------------------------------------- #
# BUSCO
# --------------------------------------------------------------------------- #
def _busco_flags(b: BuscoResult, cfg: AsmQCConfig) -> list[Flag]:
    flags: list[Flag] = []
    lineage = f" ({b.lineage})" if b.lineage else ""

    dup_warn = cfg.get("busco", "duplicated_warn_pct", 5.0)
    dup_flag = cfg.get("busco", "duplicated_flag_pct", 10.0)
    if b.duplicated_pct is not None:
        ev = {"duplicated_pct": b.duplicated_pct, "threshold": dup_warn}
        if b.duplicated_pct >= dup_flag:
            flags.append(_flag(
                "busco_duplication", Severity.FLAG, "High BUSCO duplication", "completeness",
                f"BUSCO duplicated is {b.duplicated_pct:.1f}%{lineage}, at or above the "
                f"{dup_flag:.0f}% flag threshold — likely uncollapsed haplotypes.",
                "busco_duplication", "BUSCO duplicated %", ev))
        elif b.duplicated_pct >= dup_warn:
            flags.append(_flag(
                "busco_duplication", Severity.WARN, "Elevated BUSCO duplication", "completeness",
                f"BUSCO duplicated is {b.duplicated_pct:.1f}%{lineage}, above the "
                f"{dup_warn:.0f}% warning threshold — possible haplotype duplication.",
                "busco_duplication", "BUSCO duplicated %", ev))

    comp_warn = cfg.get("busco", "complete_warn_pct", 95.0)
    comp_flag = cfg.get("busco", "complete_flag_pct", 90.0)
    if b.complete_pct is not None:
        ev = {"complete_pct": b.complete_pct, "warn": comp_warn, "flag": comp_flag}
        if b.complete_pct < comp_flag:
            flags.append(_flag(
                "busco_completeness", Severity.FLAG, "Low BUSCO completeness", "completeness",
                f"BUSCO complete is {b.complete_pct:.1f}%{lineage}, below the "
                f"{comp_flag:.0f}% flag threshold — expected genes are missing or broken.",
                "busco_completeness", "BUSCO complete %", ev))
        elif b.complete_pct < comp_warn:
            flags.append(_flag(
                "busco_completeness", Severity.WARN, "Sub-standard BUSCO completeness",
                "completeness",
                f"BUSCO complete is {b.complete_pct:.1f}%{lineage}, below the "
                f"{comp_warn:.0f}% target.",
                "busco_completeness", "BUSCO complete %", ev))

    frag_warn = cfg.get("busco", "fragmented_warn_pct", 5.0)
    if b.fragmented_pct is not None and b.fragmented_pct >= frag_warn:
        flags.append(_flag(
            "busco_fragmented", Severity.WARN, "Many fragmented BUSCOs", "completeness",
            f"BUSCO fragmented is {b.fragmented_pct:.1f}%{lineage} (>= {frag_warn:.0f}%) — "
            "genes broken across contig boundaries; tracks low contiguity.",
            "busco_fragmented", "BUSCO fragmented %",
            {"fragmented_pct": b.fragmented_pct, "threshold": frag_warn}))

    miss_warn = cfg.get("busco", "missing_warn_pct", 5.0)
    if b.missing_pct is not None and b.missing_pct >= miss_warn:
        flags.append(_flag(
            "busco_missing", Severity.WARN, "Many missing BUSCOs", "completeness",
            f"BUSCO missing is {b.missing_pct:.1f}%{lineage} (>= {miss_warn:.0f}%).",
            "busco_missing", "BUSCO missing %",
            {"missing_pct": b.missing_pct, "threshold": miss_warn}))

    if not flags and b.complete_pct is not None:
        flags.append(_flag(
            "busco_ok", Severity.PASS, "BUSCO completeness OK", "completeness",
            f"BUSCO complete {b.complete_pct:.1f}% / duplicated "
            f"{(b.duplicated_pct or 0):.1f}%{lineage} meets the configured standard.",
            "pass", "BUSCO", {"complete_pct": b.complete_pct}))
    return flags


# --------------------------------------------------------------------------- #
# Merqury
# --------------------------------------------------------------------------- #
def _merqury_flags(m: MerquryResult, cfg: AsmQCConfig) -> list[Flag]:
    flags: list[Flag] = []
    qv_warn = cfg.get("merqury", "qv_warn", 40.0)
    qv_flag = cfg.get("merqury", "qv_flag", 30.0)
    if m.qv is not None:
        ev = {"qv": m.qv, "warn": qv_warn, "flag": qv_flag}
        if m.qv < qv_flag:
            flags.append(_flag(
                "merqury_qv", Severity.FLAG, "Low consensus accuracy (QV)", "correctness",
                f"Merqury QV is {m.qv:.1f}, below the Q{qv_flag:.0f} flag threshold "
                "— substantial residual base errors.",
                "merqury_qv", "Merqury QV", ev))
        elif m.qv < qv_warn:
            flags.append(_flag(
                "merqury_qv", Severity.WARN, "Consensus accuracy below standard", "correctness",
                f"Merqury QV is {m.qv:.1f}, below the Q{qv_warn:.0f} EBP standard "
                f"(Q40 = 99.99% / 1 error per 10 kb).",
                "merqury_qv", "Merqury QV", ev))
        else:
            flags.append(_flag(
                "merqury_qv", Severity.PASS, "Consensus accuracy OK", "correctness",
                f"Merqury QV {m.qv:.1f} meets the Q{qv_warn:.0f} standard.",
                "pass", "Merqury QV", ev))

    comp_warn = cfg.get("merqury", "completeness_warn_pct", 90.0)
    if m.completeness_pct is not None and m.completeness_pct < comp_warn:
        flags.append(_flag(
            "merqury_completeness", Severity.WARN, "Low k-mer completeness", "completeness",
            f"Merqury k-mer completeness is {m.completeness_pct:.1f}%, below "
            f"{comp_warn:.0f}% — read k-mers are missing from the assembly.",
            "merqury_completeness", "k-mer completeness",
            {"completeness_pct": m.completeness_pct, "threshold": comp_warn}))
    return flags


# --------------------------------------------------------------------------- #
# Contiguity
# --------------------------------------------------------------------------- #
def _contiguity_flags(s: AssemblyStats, cfg: AsmQCConfig) -> list[Flag]:
    flags: list[Flag] = []
    sn_warn = cfg.get("contiguity", "scaffold_n50_warn_bp", 1_000_000)
    sn_note = cfg.get("contiguity", "scaffold_n50_note_bp", 10_000_000)
    if s.n50 is not None:
        ev = {"scaffold_n50": s.n50, "warn": sn_warn, "note": sn_note}
        if s.n50 < sn_warn:
            flags.append(_flag(
                "contiguity_scaffold_n50", Severity.WARN, "Low scaffold N50 (fragmented)",
                "contiguity",
                f"Scaffold N50 is {_fmt_bp(s.n50)}, below {_fmt_bp(sn_warn)} — a fragmented "
                "draft, not chromosome-scale.",
                "contiguity_scaffold_n50", "Scaffold N50", ev))
        elif s.n50 < sn_note:
            flags.append(_flag(
                "contiguity_scaffold_n50", Severity.NOTE, "Scaffold N50 below chromosome scale",
                "contiguity",
                f"Scaffold N50 is {_fmt_bp(s.n50)}, below the {_fmt_bp(sn_note)} "
                "chromosome-scale target.",
                "contiguity_scaffold_n50", "Scaffold N50", ev))
        else:
            flags.append(_flag(
                "contiguity_scaffold_n50", Severity.PASS, "Scaffold N50 OK", "contiguity",
                f"Scaffold N50 is {_fmt_bp(s.n50)} (>= {_fmt_bp(sn_note)}).",
                "pass", "Scaffold N50", ev))

    cn_warn = cfg.get("contiguity", "contig_n50_warn_bp", 100_000)
    cn_note = cfg.get("contiguity", "contig_n50_note_bp", 1_000_000)
    if s.contig_n50 is not None:
        ev = {"contig_n50": s.contig_n50, "warn": cn_warn, "note": cn_note}
        if s.contig_n50 < cn_warn:
            flags.append(_flag(
                "contiguity_contig_n50", Severity.WARN, "Low contig N50", "contiguity",
                f"Contig N50 is {_fmt_bp(s.contig_n50)}, below {_fmt_bp(cn_warn)} — the "
                "underlying contigs are short.",
                "contiguity_contig_n50", "Contig N50", ev))
        elif s.contig_n50 < cn_note:
            flags.append(_flag(
                "contiguity_contig_n50", Severity.NOTE, "Contig N50 below EBP target",
                "contiguity",
                f"Contig N50 is {_fmt_bp(s.contig_n50)}, below the {_fmt_bp(cn_note)} "
                "EBP contig target.",
                "contiguity_contig_n50", "Contig N50", ev))

    # Fragmentation: fraction of assembly length in short scaffolds.
    short_len = cfg.get("contiguity", "short_scaffold_length_bp", 10_000)
    short_frac_warn = cfg.get("contiguity", "short_scaffold_fraction_warn", 0.10)
    if s.per_seq:
        # Use the per-seq-derived total as the denominator so numerator and
        # denominator are consistent even when total_length was overridden from
        # an external QUAST/gfastats report (which may filter short sequences).
        per_seq_total = sum(x.length for x in s.per_seq)
        short_bases = sum(x.length for x in s.per_seq if x.length < short_len)
        frac = (short_bases / per_seq_total) if per_seq_total else 0.0
        if frac >= short_frac_warn:
            flags.append(_flag(
                "contiguity_fragmentation", Severity.WARN, "Fragmented: short-scaffold debris",
                "contiguity",
                f"{frac * 100:.1f}% of the assembly is in scaffolds < {_fmt_bp(short_len)} "
                f"({s.num_short_seqs} sequences) — above the {short_frac_warn * 100:.0f}% "
                "threshold.",
                "contiguity_fragmentation", "short-scaffold fraction",
                {"short_fraction": round(frac, 4), "num_short_seqs": s.num_short_seqs,
                 "threshold": short_frac_warn}))

    many = cfg.get("contiguity", "num_sequences_note", 1000)
    if s.num_sequences and s.num_sequences >= many:
        flags.append(_flag(
            "contiguity_many_sequences", Severity.NOTE, "Very many sequences", "contiguity",
            f"The assembly has {s.num_sequences:,} sequences (>= {many:,}).",
            "contiguity_many_sequences", "sequence count",
            {"num_sequences": s.num_sequences, "threshold": many}))
    return flags


# --------------------------------------------------------------------------- #
# Gaps
# --------------------------------------------------------------------------- #
def _gap_flags(s: AssemblyStats, cfg: AsmQCConfig) -> list[Flag]:
    flags: list[Flag] = []
    if not s.total_length:
        return flags
    per100 = cfg.get("gaps", "gaps_per_100kbp_note", 5.0)
    gaps_per_100kbp = 100_000.0 * s.gap_count / s.total_length if s.total_length else 0.0
    gaps_per_gbp = 1e9 * s.gap_count / s.total_length if s.total_length else 0.0
    frac_warn = cfg.get("gaps", "gap_fraction_warn", 0.05)
    gap_fraction = (s.gap_bases / s.total_length) if s.total_length else 0.0

    if gap_fraction >= frac_warn:
        flags.append(_flag(
            "gaps_density", Severity.WARN, "High gap content", "scaffolding",
            f"{gap_fraction * 100:.1f}% of the assembly is gap (N) bases across "
            f"{s.gap_count:,} gaps — above the {frac_warn * 100:.0f}% threshold.",
            "gaps_density", "gap fraction",
            {"gap_fraction": round(gap_fraction, 4), "gap_count": s.gap_count,
             "gaps_per_gbp": round(gaps_per_gbp, 1)}))
    elif s.gap_count and gaps_per_100kbp >= per100:
        flags.append(_flag(
            "gaps_density", Severity.NOTE, "Notable gap density", "scaffolding",
            f"{s.gap_count:,} gaps ({gaps_per_100kbp:.1f} per 100 kb, "
            f"{gaps_per_gbp:,.0f} per Gbp). A T2T assembly has zero gaps.",
            "gaps_density", "gap density",
            {"gap_count": s.gap_count, "gaps_per_100kbp": round(gaps_per_100kbp, 2),
             "gaps_per_gbp": round(gaps_per_gbp, 1)}))
    return flags


# --------------------------------------------------------------------------- #
# Telomeres
# --------------------------------------------------------------------------- #
def _telomere_flags(t: TelomereResult, cfg: AsmQCConfig) -> list[Flag]:
    flags: list[Flag] = []
    chrom_min = cfg.get("telomere", "chromosome_min_length_bp", 1_000_000)
    both_note = cfg.get("telomere", "both_ends_fraction_note", 0.90)
    any_warn = cfg.get("telomere", "any_end_fraction_warn", 0.50)

    chroms = [s for s in t.scaffolds if (s.length or 0) >= chrom_min]
    if not chroms:
        chroms = t.scaffolds  # fall back to all if none reach the chromosome size
    n = len(chroms)
    if n == 0:
        return flags
    n_both = sum(1 for s in chroms if s.both_ends)
    n_any = sum(1 for s in chroms if s.any_end)
    both_frac = n_both / n
    any_frac = n_any / n
    motif = f" (motif {t.repeat_motif})" if t.repeat_motif else ""
    ev = {
        "n_chromosomes": n, "n_both_ends": n_both, "n_any_end": n_any,
        "both_fraction": round(both_frac, 3), "any_fraction": round(any_frac, 3),
    }

    if any_frac < any_warn:
        flags.append(_flag(
            "telomere_low", Severity.WARN, "Telomeres detected at few ends", "completeness",
            f"Only {n_any}/{n} chromosome-scale scaffolds have a telomere at any end "
            f"({any_frac * 100:.0f}%){motif} — below {any_warn * 100:.0f}%.",
            "telomere_low", "telomere presence", ev))
    elif both_frac < both_note:
        flags.append(_flag(
            "telomere_t2t", Severity.NOTE, "Not telomere-to-telomere", "completeness",
            f"{n_both}/{n} chromosome-scale scaffolds are telomere-capped at BOTH ends "
            f"({both_frac * 100:.0f}%){motif}; {both_note * 100:.0f}% needed for a T2T claim.",
            "telomere_incomplete", "telomere completeness", ev))
    else:
        flags.append(_flag(
            "telomere_t2t", Severity.PASS, "Telomeres at both ends", "completeness",
            f"{n_both}/{n} chromosome-scale scaffolds are capped at both ends "
            f"({both_frac * 100:.0f}%){motif}.",
            "pass", "telomere completeness", ev))
    return flags


# --------------------------------------------------------------------------- #
# Contamination
# --------------------------------------------------------------------------- #
def _contamination_flags(c: ContaminationResult, cfg: AsmQCConfig) -> list[Flag]:
    flags: list[Flag] = []
    if c.source == "fcs":
        flag_actions = set(cfg.get("contamination", "flag_actions", ["EXCLUDE", "TRIM"]))
        review_actions = set(cfg.get("contamination", "review_actions", ["REVIEW", "FIX"]))
        flagged = [h for h in c.hits if (h.action or "") in flag_actions]
        review = [h for h in c.hits if (h.action or "") in review_actions]
        big = cfg.get("contamination", "flagged_bases_flag", 100_000)
        if flagged:
            severity = Severity.FLAG
            taxa = sorted(_clean_taxon(h.taxon) for h in flagged if h.taxon)
            taxa_str = ", ".join(taxa[:4]) + ("…" if len(taxa) > 4 else "")
            # Totals over the EXCLUDE/TRIM subset only (the parser's
            # total_flagged_bases / n_sequences_flagged span every action incl.
            # REVIEW/INFO, which would overstate the removal/trim span).
            flagged_bases = sum(_span(h) for h in flagged)
            flagged_seqs = len({h.seq_id for h in flagged})
            flags.append(_flag(
                "contamination_fcs", severity, "FCS-GX flagged contamination", "contamination",
                f"FCS-GX flagged {len(flagged)} region(s) across "
                f"{flagged_seqs} sequence(s), {flagged_bases:,} bp "
                f"({_fmt_bp(flagged_bases)}) for removal/trim"
                + (f"; top taxa: {taxa_str}" if taxa_str else "") + ".",
                "contamination_fcs", "FCS-GX action",
                {"n_flagged": len(flagged), "total_flagged_bases": flagged_bases,
                 "n_sequences_flagged": flagged_seqs,
                 "actions": sorted({h.action for h in flagged if h.action}),
                 "large": flagged_bases >= big}))
        if review:
            flags.append(_flag(
                "contamination_fcs_review", Severity.NOTE, "FCS-GX review items", "contamination",
                f"FCS-GX marked {len(review)} region(s) for manual REVIEW.",
                "contamination_fcs", "FCS-GX review",
                {"n_review": len(review)}))
        if not flagged and not review:
            flags.append(_flag(
                "contamination_fcs", Severity.PASS, "No FCS-GX contamination", "contamination",
                "FCS-GX flagged no sequence for removal or review.",
                "pass", "FCS-GX", {}))
    elif c.source == "kraken":
        thr = cfg.get("contamination", "kraken_taxon_flag_pct", 1.0)
        # Treat non-Eukaryota superkingdoms as potential contaminants by default.
        suspects = {
            name: pct for name, pct in c.top_taxa.items()
            if name.lower() not in ("eukaryota", "root", "cellular organisms") and pct >= thr
        }
        if suspects:
            top = ", ".join(f"{_clean_taxon(n)} {p:.1f}%" for n, p in
                            sorted(suspects.items(), key=lambda kv: -kv[1])[:4])
            flags.append(_flag(
                "contamination_kraken", Severity.FLAG, "Kraken2 contaminant signal",
                "contamination",
                f"Kraken2 assigned {top} of the assembly to non-target superkingdoms "
                f"(>= {thr:.1f}%).",
                "contamination_kraken", "Kraken2 superkingdom",
                {"suspects": {k: round(v, 2) for k, v in suspects.items()},
                 "unclassified_pct": c.unclassified_pct}))
        else:
            flags.append(_flag(
                "contamination_kraken", Severity.PASS, "No Kraken2 contaminant signal",
                "contamination",
                "Kraken2 found no non-target superkingdom above the threshold.",
                "pass", "Kraken2", {"unclassified_pct": c.unclassified_pct}))
    return flags
