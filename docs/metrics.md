# Metrics reference

A deeper companion to the README's [field guide](../README.md#reading-an-assembly-report-a-field-guide).
For each metric: the definition, how AsmQC computes/obtains it, and the
good/warn/flag bars (all configurable in `thresholds.yaml`).

## EBP metric code

The Earth BioGenome Project summarises an assembly as **`x.y.QV`**:

- `x` = ⌊log₁₀(contig N50)⌋ — `6` = ≥ 1 Mb, `7` = ≥ 10 Mb, `5` = ≥ 100 kb
- `y` = ⌊log₁₀(scaffold N50)⌋, or the letter **`C`** for chromosome-scale
- `QV` = Merqury consensus QV

Reference standard ≈ **`6.7.Q40`** (contig N50 1 Mb, scaffold N50 10 Mb, QV 40).
T2T tier = **`C.C.Q60`** (gapless, telomeres both ends, QV 60).

---

## Contiguity

| Metric | Definition | AsmQC source |
| --- | --- | --- |
| Total length | Sum of all sequence lengths | FASTA / QUAST / gfastats |
| N50 / N90 | Length where cumulative size first reaches 50% / 90% (sequences sorted longest-first) | computed |
| L50 / L90 | Number of sequences needed to reach 50% / 90% | computed |
| auN | Area under the Nx curve = Σ Lᵢ² / total — a length-threshold-free contiguity score | computed |
| Contig N50 | N50 of contigs (scaffolds split on gaps ≥ `min_gap_len`) | computed |
| Gap count / bases | Runs of ≥ `min_gap_len` N's, and total N's in them | computed / gfastats |

**Defaults:** scaffold N50 warns < 1 Mb, notes < 10 Mb; contig N50 warns < 100 kb,
notes < 1 Mb. Fragmentation warns when > 10% of total length is in scaffolds shorter
than `short_scaffold_length_bp` (10 kb). These suit vertebrate-scale genomes —
**scale them to your organism.**

> N50 vs NG50: NG50 normalises to a known genome size and is what EBP uses for
> cross-species comparison. AsmQC reports N50 (no genome size assumed) and lets you
> set absolute bp thresholds.

## Completeness

**BUSCO** (Benchmarking Universal Single-Copy Orthologs): fraction of lineage genes
found complete (C = single S + duplicated D), fragmented (F) or missing (M).

- Complete: good ≥ 95%, warn < 95%, flag < 90%.
- Duplicated: good < 5%, warn ≥ 5%, flag ≥ 10% → uncollapsed haplotigs (run purge_dups),
  unless explained by real whole-genome duplication.
- Fragmented / Missing ≥ 5% → low contiguity or wrong lineage dataset.

**Merqury k-mer completeness:** fraction of reliable read k-mers present in the
assembly (captures non-genic sequence too). Warn < 90%.

**Telomeres (tidk):** a chromosome is T2T-complete when telomere arrays cap both
ends. AsmQC counts, among scaffolds ≥ `chromosome_min_length_bp`, the fraction capped
at both ends (note if < 90%) or at any end (warn if < 50%).

## Correctness

**Merqury QV** (Phred): `QV = −10·log₁₀(E)`, where `E` is the per-base error
probability inferred from assembly-only k-mers. Q30 = 99.9%, **Q40 = 99.99% (EBP
bar)**, Q50 = 99.999%, Q60 = T2T-grade. Warn < 40, flag < 30.

**Contamination:** FCS-GX `EXCLUDE`/`TRIM` spans → flag (coordinates exported);
Kraken2 non-target superkingdom ≥ 1% → flag.

## What "good" looks like (worked numbers)

| | Reference (good) | Fragmented draft (bad) |
| --- | --- | --- |
| Scaffold N50 | 71 Mb | 0.40 Mb |
| Contig N50 | 18 Mb | 38 kb |
| # scaffolds | 164 | 337,188 |
| BUSCO complete | 98.9% | 86.2% |
| BUSCO duplicated | 0.28% | 1.2% (F 7.5%, M 6.4%) |
| Merqury QV | ≥ 40 | not measured |

These are the two real assemblies in the [cookbook](../cookbook/) (Rock Ptarmigan vs
Philippine tarsier).

## References

Lawniczak/Lewin et al. 2022 (EBP standards, PNAS); Rhie et al. 2021 (VGP, Nature);
Rhie et al. 2020 (Merqury, Genome Biology); Manni et al. 2021 (BUSCO, Current
Protocols); ERGA — Mc Cartney et al. 2024 (npj Biodiversity).
