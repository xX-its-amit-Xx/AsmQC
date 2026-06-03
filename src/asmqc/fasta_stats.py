# AsmQC — contiguity statistics computed directly from an assembly FASTA.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Compute N50/N90/L50/L90, auN, gap and contig statistics from a FASTA.

This is the fallback used when no QUAST/gfastats report is supplied.  The reader
streams one record at a time and transparently handles ``.gz`` input, so it
copes with multi-gigabase genomes without loading the whole file into memory.

A *scaffold* is a FASTA record.  A *contig* is a maximal run of non-gap sequence
within a scaffold, where a *gap* is a run of at least ``min_gap_len`` ``N`` bases
(short internal N runs stay inside a contig, matching gfastats/QUAST behaviour).
"""
from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from pathlib import Path

from asmqc.models import AssemblyStats, SeqStat


def open_text(path: str | Path):
    """Open *path* for text reading, transparently decompressing ``.gz``."""
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(record_id, sequence)`` for each record.

    ``record_id`` is the first whitespace-delimited token of the header (the
    convention every downstream tool — samtools, BUSCO, tidk — uses).
    """
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def compute_seq_stat(name: str, seq: str, min_gap_len: int) -> tuple[SeqStat, list[int]]:
    """Per-record statistics plus the list of contig lengths in that record."""
    length = len(seq)
    upper = seq.upper()
    n_bases = upper.count("N")
    gc = upper.count("G") + upper.count("C")
    at = upper.count("A") + upper.count("T")
    denom = gc + at
    gc_percent = (100.0 * gc / denom) if denom else None

    gap_re = re.compile(r"N{%d,}" % max(1, int(min_gap_len)))
    gap_count = 0
    gap_bases = 0
    for m in gap_re.finditer(upper):
        gap_count += 1
        gap_bases += m.end() - m.start()

    # Contigs = non-empty pieces after splitting on gap runs.
    contig_lengths = [len(piece) for piece in gap_re.split(upper) if piece]
    if not contig_lengths and length and gap_count == 0:
        # A non-empty scaffold with no gap run (e.g. only short sub-threshold N
        # runs) is one real contig. A scaffold that is entirely gap (one long
        # N run) is counted *only* as a gap, never as a contig — otherwise the
        # same N bases would be double-counted in gap_bases and contig_n50.
        contig_lengths = [length]

    stat = SeqStat(
        name=name,
        length=length,
        gc_percent=gc_percent,
        n_bases=n_bases,
        gap_count=gap_count,
        gap_bases=gap_bases,
        num_subcontigs=len(contig_lengths) if length else 0,
    )
    return stat, contig_lengths


def nx_lx(lengths_desc: list[int], total: int, fraction: float) -> tuple[int | None, int | None]:
    """Return ``(Nx, Lx)`` for the given fraction (0.5 -> N50/L50)."""
    if not lengths_desc or total <= 0:
        return None, None
    target = fraction * total
    cumulative = 0
    for i, length in enumerate(lengths_desc, start=1):
        cumulative += length
        if cumulative >= target:
            return length, i
    return lengths_desc[-1], len(lengths_desc)


def auN(lengths: list[int], total: int) -> float | None:
    """Area under the Nx curve: sum(L_i^2) / total.  Length-threshold free."""
    if not lengths or total <= 0:
        return None
    return sum(length * length for length in lengths) / total


def compute_assembly_stats(
    fasta_path: str | Path,
    min_gap_len: int = 10,
    short_scaffold_length: int = 10000,
    keep_per_seq: bool = True,
) -> AssemblyStats:
    """Compute a full :class:`AssemblyStats` from *fasta_path*."""
    per_seq: list[SeqStat] = []
    scaffold_lengths: list[int] = []
    contig_lengths: list[int] = []
    total_length = 0
    total_gc = 0
    total_at = 0
    total_n = 0
    gap_count = 0
    gap_bases = 0

    for name, seq in iter_fasta(fasta_path):
        stat, contigs = compute_seq_stat(name, seq, min_gap_len)
        if keep_per_seq:
            per_seq.append(stat)
        scaffold_lengths.append(stat.length)
        contig_lengths.extend(contigs)
        total_length += stat.length
        total_n += stat.n_bases
        gap_count += stat.gap_count
        gap_bases += stat.gap_bases
        # Recompute GC/AT totals from the record (cheap, avoids storing per-base).
        upper = seq.upper()
        total_gc += upper.count("G") + upper.count("C")
        total_at += upper.count("A") + upper.count("T")

    num_sequences = len(scaffold_lengths)
    stats = AssemblyStats(source="fasta")
    stats.total_length = total_length
    stats.num_sequences = num_sequences
    stats.num_scaffolds = num_sequences
    stats.num_contigs = len(contig_lengths)
    stats.gap_count = gap_count
    stats.gap_bases = gap_bases
    stats.short_seq_threshold = short_scaffold_length
    stats.per_seq = per_seq if keep_per_seq else []

    if num_sequences == 0:
        return stats

    scaffold_desc = sorted(scaffold_lengths, reverse=True)
    contig_desc = sorted(contig_lengths, reverse=True)
    stats.sorted_lengths = scaffold_desc

    stats.n50, stats.l50 = nx_lx(scaffold_desc, total_length, 0.5)
    stats.n90, stats.l90 = nx_lx(scaffold_desc, total_length, 0.9)
    contig_total = sum(contig_lengths)
    stats.contig_n50, _ = nx_lx(contig_desc, contig_total, 0.5)
    stats.auN = auN(scaffold_desc, total_length)
    stats.largest = scaffold_desc[0]
    stats.smallest = scaffold_desc[-1]
    stats.mean_length = total_length / num_sequences
    stats.median_length = _median(scaffold_desc)
    denom = total_gc + total_at
    stats.gc_percent = (100.0 * total_gc / denom) if denom else None
    stats.n_per_100kbp = (100000.0 * total_n / total_length) if total_length else None
    stats.num_short_seqs = sum(1 for length in scaffold_lengths if length < short_scaffold_length)
    return stats


def _median(values_desc: list[int]) -> int | None:
    if not values_desc:
        return None
    n = len(values_desc)
    mid = n // 2
    if n % 2 == 1:
        return values_desc[mid]
    return (values_desc[mid - 1] + values_desc[mid]) // 2
