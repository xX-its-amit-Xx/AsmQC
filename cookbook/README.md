# AsmQC cookbook — worked examples on real assemblies

Two real, published, non-model vertebrate genomes from NCBI, run through AsmQC to
show the difference between **what passing QC looks like** and **how the flags catch
a problematic draft**.

| | [Rock Ptarmigan](01_rock_ptarmigan_good/) | [Philippine tarsier](02_philippine_tarsier_draft/) |
| --- | --- | --- |
| Species | *Lagopus muta* | *Carlito syrichta* |
| Accession | `GCA_023343835.1` | `GCA_000164805.2` |
| Class | **Good** (VGP chromosome-level, 2023) | **Problematic** (short-read draft, 2013) |
| Assembly level | Chromosome (40 chromosomes) | Scaffold (no chromosomes) |
| Scaffold N50 | **71.2 Mb** | **0.40 Mb** |
| Contig N50 | **17.9 Mb** | **38 kb** |
| # scaffolds | 164 | 337,188 |
| BUSCO complete | 98.9% | 86.2% |
| BUSCO dup / frag / missing | 0.28% / 0.19% / 0.91% | 1.2% / 7.5% / 6.4% |
| **AsmQC overall** | **NOTE** (passes 3 C's; not yet T2T) | **FLAG** (low contiguity + completeness) |

Open each case's `report/report.html` in a browser, or read the `NARRATION.md` for a
guided tour of what the flags caught.

## Reproduce

```bash
python cookbook/generate.py
```

This reconstructs each tool's native output files from the assemblies' **published
metrics** and runs AsmQC on them. See [`SOURCES.md`](SOURCES.md) for every number's
provenance — which values are taken verbatim from NCBI/the paper, and which (the Rock
Ptarmigan's Merqury QV and telomere windows) are *representative* of the published
quality tier because the raw per-base QC files are not distributed.

> **Why not download the genomes and re-run BUSCO/Merqury here?** Those pipelines need
> the multi-gigabyte read sets and reference databases and hours of compute — out of
> scope for a repository. AsmQC's job starts *after* those tools run, so the cookbook
> feeds it the tools' real reported outputs. The contiguity and BUSCO numbers in both
> reports are the genuine published values.

## What to take away

- **Rock Ptarmigan** passes Contiguity (71 Mb scaffold N50), Completeness (98.9% BUSCO,
  negligible duplication) and Correctness (QV ≥ 40). The single `NOTE` is that not every
  chromosome is telomere-capped — i.e. excellent, but not yet telomere-to-telomere. This
  is exactly the profile of a modern reference genome.
- **Philippine tarsier** trips six flags: a `FLAG` for sub-90% BUSCO completeness and
  `WARN`s for fragmented/missing genes, sub-100 kb contig N50 and sub-1 Mb scaffold N50,
  plus a `NOTE` for its 337,188 sequences. The flags trace a single story: a short-read
  draft never reached chromosome scale, and that fragmentation degraded the gene space.
