# AsmQC — orchestration: parse inputs, compute stats, run flags, write outputs.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Top-level pipeline that turns a set of QC input files into a :class:`QCReport`
and (optionally) the report.html / report.md / summary.json / browser tracks.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from asmqc import __version__
from asmqc.config import AsmQCConfig
from asmqc.fasta_stats import compute_assembly_stats
from asmqc.flags import evaluate
from asmqc.models import AssemblyStats, QCReport
from asmqc.parsers import (
    parse_busco,
    parse_contamination,
    parse_gfastats,
    parse_merqury,
    parse_quast,
    parse_tidk,
)


def build_report(
    fasta: str | Path | None = None,
    busco: str | Path | None = None,
    merqury: str | Path | None = None,
    quast: str | Path | None = None,
    gfastats: str | Path | None = None,
    tidk: str | Path | None = None,
    contamination: str | Path | None = None,
    config: AsmQCConfig | None = None,
    assembly_name: str | None = None,
    kmer_size: int | None = None,
) -> QCReport:
    """Parse whatever inputs are present and assemble a flagged report."""
    config = config or AsmQCConfig.default()
    warnings: list[str] = []
    inputs: dict[str, str] = {}

    name = assembly_name or (_infer_name(fasta) if fasta else None) or "assembly"
    report = QCReport(assembly_name=name, asmqc_version=__version__)
    report.generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    report.config = config.to_dict()

    # -- Contiguity stats ---------------------------------------------------
    min_gap = int(config.get("gaps", "min_gap_len", 10))
    short_len = int(config.get("contiguity", "short_scaffold_length_bp", 10000))
    stats = _resolve_stats(fasta, quast, gfastats, min_gap, short_len, warnings, inputs)
    report.stats = stats

    # -- BUSCO --------------------------------------------------------------
    if busco:
        report.busco = _try(parse_busco, busco, "BUSCO", warnings)
        if report.busco:
            inputs["busco"] = str(busco)

    # -- Merqury ------------------------------------------------------------
    if merqury:
        report.merqury = _try(
            lambda p: parse_merqury(p, kmer_size=kmer_size), merqury, "Merqury", warnings
        )
        if report.merqury:
            inputs["merqury"] = str(merqury)

    # -- Telomeres (tidk) ---------------------------------------------------
    if tidk:
        chrom_min = int(config.get("telomere", "chromosome_min_length_bp", 1_000_000))
        min_count = int(config.get("telomere", "tidk_min_window_count", 50))
        terminal = int(config.get("telomere", "tidk_terminal_windows", 2))
        report.telomere = _try(
            lambda p: parse_tidk(
                p, min_window_count=min_count, terminal_windows=terminal,
                chromosome_min_length=chrom_min,
            ),
            tidk, "tidk", warnings,
        )
        if report.telomere:
            inputs["tidk"] = str(tidk)

    # -- Contamination ------------------------------------------------------
    if contamination:
        report.contamination = _try(parse_contamination, contamination, "contamination", warnings)
        if report.contamination:
            inputs["contamination"] = str(contamination)

    if fasta:
        inputs.setdefault("fasta", str(fasta))

    report.inputs = inputs
    report.warnings = warnings
    report.flags = evaluate(report, config)
    return report


def _resolve_stats(
    fasta, quast, gfastats, min_gap, short_len, warnings, inputs
) -> AssemblyStats | None:
    """Prefer FASTA-computed stats (richest); refine headline numbers from a
    supplied QUAST/gfastats report when available; fall back to the report when
    no FASTA is given."""
    fasta_stats: AssemblyStats | None = None
    if fasta and Path(fasta).exists():
        fasta_stats = _try(
            lambda p: compute_assembly_stats(p, min_gap_len=min_gap,
                                              short_scaffold_length=short_len),
            fasta, "FASTA", warnings,
        )

    external: AssemblyStats | None = None
    if quast:
        external = _try(parse_quast, quast, "QUAST", warnings)
        if external:
            inputs["quast"] = str(quast)
    if external is None and gfastats:
        external = _try(parse_gfastats, gfastats, "gfastats", warnings)
        if external:
            inputs["gfastats"] = str(gfastats)

    if fasta_stats and external:
        return _merge_stats(fasta_stats, external)
    return fasta_stats or external


def _merge_stats(base: AssemblyStats, ext: AssemblyStats) -> AssemblyStats:
    """Keep FASTA per-seq richness; override headline N50-family numbers from the
    report.  Deliberately does NOT override ``num_contigs``/``contig_n50``: the
    FASTA computes those from real gap structure, whereas QUAST/gfastats either
    cannot (QUAST) or already agree (gfastats), so the FASTA values are kept to
    avoid the contig vs scaffold count being polluted by a filtered-sequence
    count."""
    for field_name in ("total_length", "n50", "n90", "l50", "l90",
                       "auN", "largest", "gc_percent", "n_per_100kbp"):
        value = getattr(ext, field_name, None)
        if value is not None:
            setattr(base, field_name, value)
    base.source = f"fasta+{ext.source}"
    return base


def _infer_name(path: str | Path) -> str:
    stem = Path(path).name
    for suffix in (".gz", ".fasta", ".fa", ".fna", ".fas"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
    return stem or "assembly"


def _try(func, arg, label, warnings):
    try:
        return func(arg)
    except Exception as exc:  # never let one bad input sink the report
        warnings.append(f"Failed to parse {label} input ({arg}): {exc}")
        return None


# --------------------------------------------------------------------------- #
# Output writing
# --------------------------------------------------------------------------- #
def write_outputs(
    report: QCReport,
    out_dir: str | Path,
    config: AsmQCConfig,
    fasta: str | Path | None = None,
    formats: tuple[str, ...] = ("html", "md", "json"),
    make_tracks: bool = False,
    make_jbrowse: bool = False,
    include_per_seq: bool = False,
) -> dict[str, Path]:
    """Render the requested outputs into *out_dir*; return written paths."""
    from asmqc import report as report_mod  # local import to avoid cycles

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    plots: dict[str, str] = {}
    if "html" in formats:
        from asmqc.plots import build_plots
        plots = build_plots(report, config)
        written["html"] = report_mod.write_html(report, out_dir / "report.html", plots)
    if "md" in formats:
        written["md"] = report_mod.write_markdown(report, out_dir / "report.md")
    if "json" in formats:
        written["json"] = report_mod.write_json(
            report, out_dir / "summary.json", include_per_seq=include_per_seq
        )

    if make_tracks or make_jbrowse:
        from asmqc.exporters import write_jbrowse_config, write_tracks
        min_gap = int(config.get("gaps", "min_gap_len", 10))
        tracks = write_tracks(report, out_dir, fasta_path=fasta, min_gap_len=min_gap)
        written["gff3"] = tracks["gff3"]
        written["bed"] = tracks["bed"]
        if make_jbrowse and fasta:
            written["jbrowse"] = write_jbrowse_config(
                out_dir, report.assembly_name, fasta,
                gff3_name="flags.gff3", bed_name="flags.bed",
            )
    return written


def run_report(
    out_dir: str | Path,
    config_path: str | Path | None = None,
    formats: tuple[str, ...] = ("html", "md", "json"),
    make_tracks: bool = False,
    make_jbrowse: bool = False,
    include_per_seq: bool = False,
    **inputs,
) -> tuple[QCReport, dict[str, Path]]:
    """Convenience: build the report and write outputs in one call."""
    config = AsmQCConfig.load(config_path)
    report = build_report(config=config, **inputs)
    written = write_outputs(
        report, out_dir, config,
        fasta=inputs.get("fasta"), formats=formats,
        make_tracks=make_tracks, make_jbrowse=make_jbrowse,
        include_per_seq=include_per_seq,
    )
    return report, written
