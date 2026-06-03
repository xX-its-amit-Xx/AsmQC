#!/usr/bin/env python3
# AsmQC — generate the illustrative example dataset and report.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Build a small, fully synthetic example assembly + QC inputs and run AsmQC on
it, writing the report into ``examples/synthetic/report/`` and the four annotated
figures into ``docs/img/``.

The data here is invented purely to exercise every code path and produce nice
figures for the README — it is NOT a real genome.  For worked examples on *real*
published assemblies, see ``cookbook/``.

Run from the repo root::

    python examples/make_example.py
"""
from __future__ import annotations

import base64
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "examples" / "synthetic"
INP = EX / "inputs"
IMG = ROOT / "docs" / "img"

random.seed(42)


def rand_seq(n: int, gc: float = 0.42) -> str:
    out = []
    for _ in range(n):
        out.append(random.choice("GC") if random.random() < gc else random.choice("AT"))
    return "".join(out)


def write_fasta() -> Path:
    """Five chromosome-scale scaffolds (one with two gaps) plus short debris and
    one bacterial-looking contaminant contig."""
    INP.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    sizes = [120_000, 100_000, 80_000, 60_000, 40_000]
    for i, size in enumerate(sizes, start=1):
        seq = rand_seq(size)
        if i == 2:  # chr2 gets two scaffolding gaps
            seq = seq[:30_000] + "N" * 200 + seq[30_000:60_000] + "N" * 500 + seq[60_000:]
        records.append((f"chr{i}", seq))
    # Short unplaced debris (fragmentation signal).
    for j in range(20):
        records.append((f"scaffold_{j:03d}", rand_seq(random.randint(800, 2500))))
    # A contaminant contig (will be flagged by the FCS input).
    records.append(("contig_bact01", rand_seq(6000, gc=0.62)))

    fa = INP / "assembly.fasta"
    with fa.open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for k in range(0, len(seq), 80):
                fh.write(seq[k:k + 80] + "\n")
    return fa


def write_busco() -> Path:
    p = INP / "short_summary.txt"
    p.write_text(
        "# BUSCO version is: 5.4.3\n"
        "# The lineage dataset is: aves_odb10 (Creation date: 2024-01-08, "
        "number of genomes: 50, number of BUSCOs: 8338)\n"
        "# BUSCO was run in mode: genome\n\n"
        "\t***** Results: *****\n\n"
        "\tC:96.4%[S:90.1%,D:6.3%],F:1.6%,M:2.0%,n:8338\n"
        "\t8038\tComplete BUSCOs (C)\n"
        "\t7513\tComplete and single-copy BUSCOs (S)\n"
        "\t525\tComplete and duplicated BUSCOs (D)\n"
        "\t133\tFragmented BUSCOs (F)\n"
        "\t167\tMissing BUSCOs (M)\n"
        "\t8338\tTotal BUSCO groups searched\n"
    )
    return p


def write_merqury() -> Path:
    mq = INP / "merqury"
    mq.mkdir(parents=True, exist_ok=True)
    (mq / "asm.qv").write_text("asm\t251000\t410000000\t38.6\t0.000138\n")
    (mq / "asm.completeness.stats").write_text("asm\tall\t96100000\t100000000\t96.10\n")
    # A spectra-cn histogram: a 1-copy peak, a small 2-copy (residual dup) peak,
    # a read-only (missing) tail and a high-copy tail.
    lines = ["Copies\tkmer_multiplicity\tCount"]
    for mult in range(1, 80):
        ro = max(0, int(40000 * 2.71 ** (-mult / 3)))
        lines.append(f"read-only\t{mult}\t{ro}")
    for mult in range(1, 80):
        # main heterozygous/haploid peak around 30x
        import math
        one = int(120000 * math.exp(-((mult - 30) ** 2) / (2 * 6 ** 2)))
        lines.append(f"1\t{mult}\t{one}")
    for mult in range(1, 80):
        import math
        two = int(20000 * math.exp(-((mult - 60) ** 2) / (2 * 7 ** 2)))
        lines.append(f"2\t{mult}\t{two}")
    for mult in range(1, 80):
        lines.append(f">4\t{mult}\t{max(0, int(1500 * 2.71 ** (-mult / 25)))}")
    (mq / "asm.spectra-cn.hist").write_text("\n".join(lines) + "\n")
    return mq


def write_tidk() -> Path:
    """Five chromosomes; chr1-3 capped both ends, chr4 one end, chr5 none."""
    p = INP / "asm_telomeric_repeat_windows.tsv"
    rows = ["id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat"]
    plan = {
        "chr1": (120_000, True, True),
        "chr2": (100_000, True, True),
        "chr3": (80_000, True, True),
        "chr4": (60_000, True, False),
        "chr5": (40_000, False, False),
    }
    for name, (length, start, end) in plan.items():
        win = 10_000
        ends = list(range(win, length, win)) + [length]
        for w in ends:
            fwd, rev = 0, 0
            if start and w == ends[0]:
                rev = 410
            if end and w == ends[-1]:
                fwd = 395
            rows.append(f"{name}\t{w}\t{fwd}\t{rev}\tTTAGGG")
    p.write_text("\n".join(rows) + "\n")
    return p


def write_fcs() -> Path:
    p = INP / "assembly.fcs_gx_report.txt"
    p.write_text(
        '##[["FCS genome report", 2, 1], {"git-rev": "v0.5.4"}]\n'
        "#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n"
        "contig_bact01\t1\t6000\t6000\tEXCLUDE\tprok:a-proteobacteria\t18\t"
        "Bradyrhizobium sp.\n"
    )
    return p


def write_example_config() -> Path:
    """Thresholds scaled for this miniature genome (so it reads like a
    chromosome-scale assembly despite small absolute sizes)."""
    p = EX / "example_thresholds.yaml"
    p.write_text(
        "# Thresholds scaled down for the miniature example genome.\n"
        "contiguity:\n"
        "  scaffold_n50_warn_bp: 20000\n"
        "  scaffold_n50_note_bp: 80000\n"
        "  contig_n50_warn_bp: 10000\n"
        "  contig_n50_note_bp: 50000\n"
        "  short_scaffold_length_bp: 5000\n"
        "telomere:\n"
        "  chromosome_min_length_bp: 30000\n"
    )
    return p


def main() -> None:
    from asmqc.config import AsmQCConfig
    from asmqc.core import build_report, write_outputs
    from asmqc.plots import build_plots

    fa = write_fasta()
    busco = write_busco()
    merqury = write_merqury()
    tidk = write_tidk()
    fcs = write_fcs()
    cfg_path = write_example_config()

    def rel(p):
        # repo-relative POSIX paths so the committed report shows clean paths
        return p.resolve().relative_to(ROOT).as_posix()

    cfg = AsmQCConfig.load(cfg_path)
    report = build_report(
        fasta=rel(fa), busco=rel(busco), merqury=rel(merqury), tidk=rel(tidk),
        contamination=rel(fcs), config=cfg, assembly_name="example_bird_v1",
    )
    out = EX / "report"
    write_outputs(report, out, cfg, fasta=rel(fa),
                  formats=("html", "md", "json"), make_tracks=True, make_jbrowse=True)

    # Save standalone PNGs for the README.
    IMG.mkdir(parents=True, exist_ok=True)
    plots = build_plots(report, cfg)
    for name, b64 in plots.items():
        (IMG / f"example_{name}.png").write_bytes(base64.b64decode(b64))

    print(f"Overall status: {report.overall_severity.name}")
    print("Flags:", {k: v for k, v in report.counts_by_severity().items() if v})
    print(f"Report written to {out}")
    print(f"Figures written to {IMG}")


if __name__ == "__main__":
    main()
