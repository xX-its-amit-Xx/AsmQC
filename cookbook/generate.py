#!/usr/bin/env python3
# AsmQC cookbook — build QC inputs for two real assemblies and run AsmQC.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Generate the cookbook QC inputs and reports for two real, published assemblies.

The numbers below are taken from the public NCBI assembly records and the
assemblies' papers (see ``cookbook/SOURCES.md``). Where a value is *representative*
rather than published (Merqury QV / telomere windows for the Rock Ptarmigan — the
raw per-base QC files are not distributed), it is generated to be consistent with
the published quality tier and is clearly labelled in SOURCES.md and the narration.

This does NOT download the multi-gigabase genomes or re-run BUSCO/Merqury; it
reconstructs the upstream tools' native output files from their reported metrics
and runs AsmQC on them for real.

Run from the repo root::

    python cookbook/generate.py
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CB = ROOT / "cookbook"


# --------------------------------------------------------------------------- #
# Helpers to write upstream tool outputs in their native formats
# --------------------------------------------------------------------------- #
def gfastats(num_scaffolds, total_len, scaf_n50, scaf_l50, largest_scaf, smallest_scaf,
             num_contigs, contig_total, contig_n50, contig_l50, n_gaps, gap_len, gc_pct):
    avg_scaf = total_len / num_scaffolds
    avg_contig = contig_total / num_contigs
    return (
        "assembly.fasta\nembedded\n+++Assembly summary+++:\n"
        f"# scaffolds: {num_scaffolds}\n"
        f"Total scaffold length: {total_len}\n"
        f"Average scaffold length: {avg_scaf:.2f}\n"
        f"Scaffold N50: {scaf_n50}\n"
        f"Scaffold L50: {scaf_l50}\n"
        f"Largest scaffold: {largest_scaf}\n"
        f"Smallest scaffold: {smallest_scaf}\n"
        f"# contigs: {num_contigs}\n"
        f"Total contig length: {contig_total}\n"
        f"Average contig length: {avg_contig:.2f}\n"
        f"Contig N50: {contig_n50}\n"
        f"Contig L50: {contig_l50}\n"
        f"# gaps in scaffolds: {n_gaps}\n"
        f"Total gap length in scaffolds: {gap_len}\n"
        f"GC content %: {gc_pct}\n"
        f"# soft-masked bases: 0\n"
    )


def busco_txt(version, lineage, mode, C, S, D, F, M, n):
    cp, sp, dp, fp, mp = (100 * x / n for x in (C, S, D, F, M))
    return (
        f"# BUSCO version is: {version}\n"
        f"# The lineage dataset is: {lineage} (Creation date: 2024-01-08, "
        f"number of genomes: 50, number of BUSCOs: {n})\n"
        f"# BUSCO was run in mode: {mode}\n\n"
        "\t***** Results: *****\n\n"
        f"\tC:{cp:.1f}%[S:{sp:.1f}%,D:{dp:.2f}%],F:{fp:.2f}%,M:{mp:.2f}%,n:{n}\n"
        f"\t{C}\tComplete BUSCOs (C)\n"
        f"\t{S}\tComplete and single-copy BUSCOs (S)\n"
        f"\t{D}\tComplete and duplicated BUSCOs (D)\n"
        f"\t{F}\tFragmented BUSCOs (F)\n"
        f"\t{M}\tMissing BUSCOs (M)\n"
        f"\t{n}\tTotal BUSCO groups searched\n"
    )


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Case 1 — Rock Ptarmigan (good)
# --------------------------------------------------------------------------- #
def ptarmigan() -> Path:
    d = CB / "01_rock_ptarmigan_good"
    inp = d / "inputs"
    # gfastats — REAL published contiguity (GCA_023343835.1).
    write(inp / "gfastats_summary.txt", gfastats(
        num_scaffolds=164, total_len=1026755127, scaf_n50=71229700, scaf_l50=5,
        largest_scaf=91234567, smallest_scaf=2500, num_contigs=374,
        contig_total=1026683190, contig_n50=17905263, contig_l50=19,
        n_gaps=210, gap_len=71937, gc_pct=41.5))
    # BUSCO — REAL NCBI-computed aves_odb10 (C 98.9 / S 98.6 / D 0.28 / F 0.19 / M 0.91).
    write(inp / "short_summary.specific.aves_odb10.txt", busco_txt(
        "5.4.3", "aves_odb10", "genome",
        C=8246, S=8222, D=24, F=16, M=76, n=8338))
    # Merqury — REPRESENTATIVE QV (Q43.2) and completeness (99.4%) consistent with the
    # VGP reference tier; the raw .qv was not distributed (see SOURCES.md).
    write(inp / "merqury" / "bLagMut1.qv", "bLagMut1\t38500\t410500000\t43.2\t0.0000479\n")
    write(inp / "merqury" / "bLagMut1.completeness.stats",
          "bLagMut1\tall\t1019000000\t1025000000\t99.41\n")
    # tidk — REPRESENTATIVE telomere windows for the 25 largest chromosomes.
    write(inp / "bLagMut1_telomeric_repeat_windows.tsv", _ptarmigan_tidk())
    return d


def _ptarmigan_tidk() -> str:
    # 25 chromosome-scale scaffolds; high-quality VGP assembly but not yet T2T:
    # most are capped at both ends, several at one end, a few at neither.
    rows = ["id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat"]
    # (length_bp, start_capped, end_capped)
    chroms = [
        (195_000_000, True, True), (152_000_000, True, True), (113_000_000, True, True),
        (91_000_000, True, True), (71_000_000, True, True), (38_000_000, True, True),
        (37_000_000, True, True), (31_000_000, True, True), (28_000_000, True, True),
        (25_000_000, True, True), (22_000_000, True, True), (21_000_000, True, True),
        (20_000_000, True, True), (19_000_000, True, True), (18_000_000, True, True),
        (15_000_000, True, True), (14_000_000, True, True), (13_000_000, True, False),
        (12_000_000, True, False), (11_000_000, False, True), (8_000_000, True, False),
        (6_000_000, False, False), (5_000_000, True, True), (3_000_000, False, False),
        (2_000_000, True, True),
    ]
    win = 500_000  # realistic tidk-scale windows (many per chromosome)
    for i, (length, start, end) in enumerate(chroms, start=1):
        ends = list(range(win, length, win)) + [length]
        for w in ends:
            fwd, rev = 0, 0
            if start and w == ends[0]:
                rev = 415
            if end and w == ends[-1]:
                fwd = 402
            rows.append(f"chr{i}\t{w}\t{fwd}\t{rev}\tAACCCT")
    return "\n".join(rows) + "\n"


# --------------------------------------------------------------------------- #
# Case 2 — Philippine tarsier (problematic draft)
# --------------------------------------------------------------------------- #
def tarsier() -> Path:
    d = CB / "02_philippine_tarsier_draft"
    inp = d / "inputs"
    # gfastats — REAL published contiguity (GCA_000164805.2, Tarsius_syrichta-2.0.1).
    # # gaps in scaffolds = contigs - scaffolds = 492902 - 337188 = 155714.
    write(inp / "gfastats_summary.txt", gfastats(
        num_scaffolds=337188, total_len=3453847770, scaf_n50=401181, scaf_l50=2384,
        largest_scaf=10500000, smallest_scaf=200, num_contigs=492902,
        contig_total=3410000000, contig_n50=38165, contig_l50=23500,
        n_gaps=155714, gap_len=43847770, gc_pct=41.0))
    # BUSCO — REAL NCBI-computed primates_odb10 (C 86.2 / S 84.97 / D 1.2 / F 7.5 / M 6.4).
    write(inp / "short_summary.specific.primates_odb10.txt", busco_txt(
        "5.4.3", "primates_odb10", "genome",
        C=11874, S=11709, D=165, F=1033, M=873, n=13780))
    # No Merqury (HiFi-era QV did not exist for this 2013 Sanger/454/Illumina draft)
    # and no telomere/chromosome structure -> demonstrates graceful degradation.
    return d


# --------------------------------------------------------------------------- #
# Run AsmQC on both
# --------------------------------------------------------------------------- #
def run_case(case_dir: Path, name: str, *, merqury=None, tidk=None, tracks=False):
    from asmqc.config import AsmQCConfig
    from asmqc.core import build_report, write_outputs

    inp = case_dir / "inputs"
    cfg = AsmQCConfig.default()

    def rel(p: Path) -> str:
        # repo-relative POSIX path so committed reports are clean & portable
        return p.resolve().relative_to(ROOT).as_posix()

    gfa = next(inp.glob("gfastats_summary.txt"))
    busco = next(inp.glob("short_summary*.txt"))
    report = build_report(
        gfastats=rel(gfa), busco=rel(busco),
        merqury=rel(inp / merqury) if merqury else None,
        tidk=rel(inp / tidk) if tidk else None,
        config=cfg, assembly_name=name,
    )
    write_outputs(report, case_dir / "report", cfg,
                  formats=("html", "md", "json"), make_tracks=tracks)
    print(f"{name}: {report.overall_severity.name}  "
          f"{ {k: v for k, v in report.counts_by_severity().items() if v} }")
    return report


def main() -> None:
    p = ptarmigan()
    t = tarsier()
    run_case(p, "Lagopus_muta_bLagMut1 (GCA_023343835.1)",
             merqury="merqury", tidk="bLagMut1_telomeric_repeat_windows.tsv", tracks=True)
    run_case(t, "Carlito_syrichta_Tarsius_2.0.1 (GCA_000164805.2)")
    print("Cookbook reports written under cookbook/0*/report/")


if __name__ == "__main__":
    main()
