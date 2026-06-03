# AsmQC report — Carlito_syrichta_Tarsius_2.0.1 (GCA_000164805.2)

**Overall status:** FLAG  

_Generated 2026-06-03 06:42:23 UTC by AsmQC 0.1.0_


**Flag summary:** NOTE: 1 · WARN: 4 · FLAG: 1


## Curation flags

### [FLAG] Low BUSCO completeness

- **Observation:** BUSCO complete is 86.2% (primates_odb10), below the 90% flag threshold — expected genes are missing or broken.
- **Likely cause:** A low Complete percentage means expected single-copy genes are not found intact in the assembly. Causes: genuinely missing sequence (low coverage, collapsed repeats), genes broken across contig boundaries, or — for distant lineages — real biological gene loss. High Missing+Fragmented usually points to a fragmented or incomplete assembly rather than biology.
- **Suggested action:** Check coverage and k-mer completeness (Merqury) to see if sequence is missing; consider additional long-read data or a re-assembly. Confirm the BUSCO lineage dataset is appropriate for the species before concluding genes are truly absent.

### [WARN] Many fragmented BUSCOs

- **Observation:** BUSCO fragmented is 7.5% (primates_odb10) (>= 5%) — genes broken across contig boundaries; tracks low contiguity.
- **Likely cause:** Fragmented BUSCOs are found only partially — typically a gene split across two contigs/scaffolds because the assembly breaks inside it. This tracks low contiguity: short contigs interrupt gene models.
- **Suggested action:** Improve contiguity (long-read scaffolding, gap filling). Fragmentation usually falls as contig N50 rises; re-evaluate after scaffolding.

### [WARN] Many missing BUSCOs

- **Observation:** BUSCO missing is 6.3% (primates_odb10) (>= 5%).
- **Likely cause:** Missing BUSCOs are not detected at all. Combined with low completeness this indicates absent sequence; alone it can also reflect an over-stringent or mismatched lineage dataset.
- **Suggested action:** Verify the lineage dataset, then investigate coverage gaps. Persistent high Missing after good coverage suggests real dropout regions to target with additional data.

### [WARN] Low contig N50

- **Observation:** Contig N50 is 38.2 kb, below 100.0 kb — the underlying contigs are short.
- **Likely cause:** Contig N50 measures the underlying contiguity before scaffolding (gaps removed). Short contigs mean the primary assembly is broken — usually short or noisy reads, or unresolved repeats. Scaffolding can join these but cannot recover the missing sequence in the gaps.
- **Suggested action:** Prefer long high-accuracy reads (PacBio HiFi / ONT) and re-assemble; contig N50 is set at the contigging stage, not by scaffolding.

### [WARN] Low scaffold N50 (fragmented)

- **Observation:** Scaffold N50 is 401.2 kb, below 1.00 Mb — a fragmented draft, not chromosome-scale.
- **Likely cause:** Scaffold N50 summarises how long the scaffolds are: half the genome sits in scaffolds at least this long. A low value means the assembly has not reached chromosome scale — it is a fragmented draft rather than a reference. The EBP target is chromosome-scale scaffolds (often 10 Mb+ N50).
- **Suggested action:** Add long-range information (Hi-C, optical maps, or linkage) to scaffold contigs into chromosomes, then re-assess. Confirm against the expected karyotype.

### [NOTE] Very many sequences

- **Observation:** The assembly has 337,188 sequences (>= 1,000).
- **Likely cause:** The assembly contains a very large number of sequences for its size. On its own this is a smell rather than a verdict, but it usually accompanies fragmentation, retained haplotigs, or contaminant contigs.
- **Suggested action:** Cross-check with the fragmentation, duplication and contamination flags; reduce the unplaced/short fraction during curation.


## Metrics

| Group | Metric | Value |
| --- | --- | --- |
| Contiguity | Total length | 3.45 Gb |
| Contiguity | Sequences (scaffolds) | 337,188 |
| Contiguity | Contigs | 492,902 |
| Contiguity | Scaffold N50 | 401.2 kb |
| Contiguity | Scaffold N90 | — |
| Contiguity | L50 | 2,384 |
| Contiguity | Contig N50 | 38.2 kb |
| Contiguity | Largest | 10.50 Mb |
| Contiguity | GC content | 41.0% |
| Contiguity | Gaps | 155,714 |
| Contiguity | Gap bases | 43.85 Mb |
| Contiguity | N per 100 kbp | 1,269.53 |
| Completeness (BUSCO) | Complete | 86.2% |
| Completeness (BUSCO) | Single-copy | 85.0% |
| Completeness (BUSCO) | Duplicated | 1.2% |
| Completeness (BUSCO) | Fragmented | 7.5% |
| Completeness (BUSCO) | Missing | 6.3% |
| Completeness (BUSCO) | Lineage | primates_odb10 |

## Inputs

- **gfastats:** cookbook/02_philippine_tarsier_draft/inputs/gfastats_summary.txt
- **busco:** cookbook/02_philippine_tarsier_draft/inputs/short_summary.specific.primates_odb10.txt

---

Thresholds follow the Earth BioGenome Project / VGP 3C framework (Contiguity, Completeness, Correctness). See the AsmQC README for what each metric means and what good vs bad looks like.

