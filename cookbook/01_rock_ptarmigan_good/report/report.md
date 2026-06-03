# AsmQC report — Lagopus_muta_bLagMut1 (GCA_023343835.1)

**Overall status:** NOTE  

_Generated 2026-06-03 06:42:22 UTC by AsmQC 0.1.0_


**Flag summary:** PASS: 3 · NOTE: 1


## Curation flags

### [NOTE] Not telomere-to-telomere

- **Observation:** 19/25 chromosome-scale scaffolds are telomere-capped at BOTH ends (76%) (motif AACCCT); 90% needed for a T2T claim.
- **Likely cause:** Telomeric repeat arrays were not detected at both ends of most chromosome-scale scaffolds. That means the assembly is not yet telomere-to-telomere: some chromosome ends are missing or unresolved (telomeres are highly repetitive and hard to assemble). This is expected for most non-T2T assemblies and is a 'note', not a failure.
- **Suggested action:** If T2T completeness is a goal, target the missing ends with ultra-long reads and verify with the telomere map. Otherwise record which chromosomes are arm-to-arm complete.

### [PASS] BUSCO completeness OK

- **Observation:** BUSCO complete 98.9% / duplicated 0.3% (aves_odb10) meets the configured standard.
- **Likely cause:** This metric meets the configured quality standard.
- **Suggested action:** No action required.

### [PASS] Scaffold N50 OK

- **Observation:** Scaffold N50 is 71.23 Mb (>= 10.00 Mb).
- **Likely cause:** This metric meets the configured quality standard.
- **Suggested action:** No action required.

### [PASS] Consensus accuracy OK

- **Observation:** Merqury QV 43.2 meets the Q40 standard.
- **Likely cause:** This metric meets the configured quality standard.
- **Suggested action:** No action required.


## Metrics

| Group | Metric | Value |
| --- | --- | --- |
| Contiguity | Total length | 1.03 Gb |
| Contiguity | Sequences (scaffolds) | 164 |
| Contiguity | Contigs | 374 |
| Contiguity | Scaffold N50 | 71.23 Mb |
| Contiguity | Scaffold N90 | — |
| Contiguity | L50 | 5 |
| Contiguity | Contig N50 | 17.91 Mb |
| Contiguity | Largest | 91.23 Mb |
| Contiguity | GC content | 41.5% |
| Contiguity | Gaps | 210 |
| Contiguity | Gap bases | 71.9 kb |
| Contiguity | N per 100 kbp | 7.01 |
| Completeness (BUSCO) | Complete | 98.9% |
| Completeness (BUSCO) | Single-copy | 98.6% |
| Completeness (BUSCO) | Duplicated | 0.3% |
| Completeness (BUSCO) | Fragmented | 0.2% |
| Completeness (BUSCO) | Missing | 0.9% |
| Completeness (BUSCO) | Lineage | aves_odb10 |
| Correctness (Merqury) | Consensus QV | 43.20 |
| Correctness (Merqury) | k-mer completeness | 99.4% |
| Telomeres | Scaffolds w/ both ends | 19 |
| Telomeres | Scaffolds w/ any end | 23 |
| Telomeres | Repeat motif | AACCCT |

## Inputs

- **gfastats:** cookbook/01_rock_ptarmigan_good/inputs/gfastats_summary.txt
- **busco:** cookbook/01_rock_ptarmigan_good/inputs/short_summary.specific.aves_odb10.txt
- **merqury:** cookbook/01_rock_ptarmigan_good/inputs/merqury
- **tidk:** cookbook/01_rock_ptarmigan_good/inputs/bLagMut1_telomeric_repeat_windows.tsv

---

Thresholds follow the Earth BioGenome Project / VGP 3C framework (Contiguity, Completeness, Correctness). See the AsmQC README for what each metric means and what good vs bad looks like.

