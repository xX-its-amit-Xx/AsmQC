# AsmQC — command-line interface.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""``asmqc`` command-line entry point (argparse subcommands).

Examples::

    asmqc run --fasta asm.fa --busco short_summary.txt --merqury merqury/ --out report/
    asmqc run --fasta asm.fa --quast quast/report.tsv --tidk windows.tsv --jbrowse --out report/
    asmqc init-config > thresholds.yaml
    asmqc version
"""
from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from asmqc import __version__
from asmqc.config import AsmQCConfig, default_config_text
from asmqc.core import build_report, write_outputs
from asmqc.models import Severity


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="asmqc",
        description="Aggregate genome-assembly QC outputs into one report with curation flags.",
    )
    parser.add_argument("--version", action="version", version=f"AsmQC {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>")

    # -- run ----------------------------------------------------------------
    run = sub.add_parser("run", help="build a consolidated QC report")
    run.add_argument("--fasta", help="assembly FASTA (.fa/.fasta[.gz]); used to compute "
                                     "contiguity stats when no QUAST/gfastats is given")
    run.add_argument("--busco", help="BUSCO short_summary .txt or .json")
    run.add_argument("--merqury", help="Merqury output dir, prefix, or .qv file")
    run.add_argument("--quast", help="QUAST report.tsv or output directory")
    run.add_argument("--gfastats", help="gfastats summary output file")
    run.add_argument("--tidk", help="tidk *_telomeric_repeat_windows.tsv")
    run.add_argument("--contamination", "--fcs", "--kraken", dest="contamination",
                     help="FCS-GX report or Kraken2 report (auto-detected)")
    run.add_argument("--out", "-o", required=True, help="output directory for the report")
    run.add_argument("--config", "-c", help="YAML thresholds (merged over defaults)")
    run.add_argument("--name", help="assembly name shown in the report")
    run.add_argument("--kmer-size", type=int, help="k-mer size used for Merqury (for context)")
    run.add_argument("--formats", default="html,md,json",
                     help="comma list of outputs: html,md,json (default: all)")
    run.add_argument("--tracks", action="store_true",
                     help="also export flags.gff3 and flags.bed")
    run.add_argument("--jbrowse", action="store_true",
                     help="also export a JBrowse2 config (implies --tracks; needs --fasta)")
    run.add_argument("--per-seq", action="store_true",
                     help="include per-sequence detail in summary.json")
    run.add_argument("--fail-on", choices=[s.name for s in Severity], default=None,
                     help="exit non-zero if the overall status reaches this severity "
                          "(e.g. FLAG) — for CI gating; default: always exit 0 on success")
    run.add_argument("--quiet", "-q", action="store_true", help="suppress the stdout summary")
    run.set_defaults(func=_cmd_run)

    # -- init-config --------------------------------------------------------
    init = sub.add_parser("init-config", help="print the default thresholds YAML")
    init.add_argument("--out", "-o", help="write to this file instead of stdout")
    init.set_defaults(func=_cmd_init_config)

    # -- version ------------------------------------------------------------
    ver = sub.add_parser("version", help="print the AsmQC version")
    ver.set_defaults(func=_cmd_version)

    return parser


def _cmd_run(args: argparse.Namespace) -> int:
    if not any([args.fasta, args.busco, args.merqury, args.quast,
                args.gfastats, args.tidk, args.contamination]):
        print("error: provide at least --fasta or one QC input", file=sys.stderr)
        return 2

    config = AsmQCConfig.load(args.config) if args.config else AsmQCConfig.default()
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())

    report = build_report(
        fasta=args.fasta,
        busco=args.busco,
        merqury=args.merqury,
        quast=args.quast,
        gfastats=args.gfastats,
        tidk=args.tidk,
        contamination=args.contamination,
        config=config,
        assembly_name=args.name,
        kmer_size=args.kmer_size,
    )
    written = write_outputs(
        report, args.out, config,
        fasta=args.fasta, formats=formats,
        make_tracks=args.tracks or args.jbrowse,
        make_jbrowse=args.jbrowse,
        include_per_seq=args.per_seq,
    )

    if not args.quiet:
        _print_summary(report, written)

    # A successful run exits 0; --fail-on lets CI gate on report severity.
    if args.fail_on:
        threshold = Severity.from_str(args.fail_on)
        if report.overall_severity.value >= threshold.value:
            return 1
    return 0


def _cmd_init_config(args: argparse.Namespace) -> int:
    text = default_config_text()
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Wrote default config to {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


def _cmd_version(args: argparse.Namespace) -> int:
    print(f"AsmQC {__version__}")
    return 0


def _print_summary(report, written: dict) -> None:
    line = "=" * 60
    print(line)
    print(f"AsmQC report: {report.assembly_name}")
    print(f"Overall status: {report.overall_severity.name}")
    counts = report.counts_by_severity()
    shown = " ".join(f"{k}={v}" for k, v in counts.items() if v)
    print(f"Flags: {shown or 'none'}")
    if report.flags:
        print("-" * 60)
        for f in report.flags:
            if f.severity.value >= Severity.NOTE.value:
                print(f"  [{f.severity.name:>4}] {f.title}: {f.message}")
    if report.warnings:
        print("-" * 60)
        for w in report.warnings:
            print(f"  ! {w}")
    print("-" * 60)
    for kind, path in written.items():
        print(f"  {kind:>8}: {path}")
    print(line)


def _make_stdio_safe() -> None:
    """Never crash on a legacy console: replace unencodable chars instead of
    raising UnicodeEncodeError (e.g. em-dashes on a Windows cp1252 console)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")  # type: ignore[union-attr]
        except (AttributeError, ValueError):
            pass


def main(argv: Sequence[str] | None = None) -> int:
    _make_stdio_safe()
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
