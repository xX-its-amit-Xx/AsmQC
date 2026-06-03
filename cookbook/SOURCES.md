# Cookbook data provenance

Every number AsmQC reports for these two assemblies is reconstructed into the
upstream tools' native output format by [`generate.py`](generate.py). This page
states the provenance of each value: **REAL** (verbatim from NCBI or the paper) vs
**REPRESENTATIVE** (generated to be consistent with the published quality tier
because the raw per-base file is not publicly distributed).

---

## 1. Rock Ptarmigan — *Lagopus muta* — `GCA_023343835.1` (bLagMut1)

| Metric | Value | Provenance |
| --- | --- | --- |
| Total length | 1,026,755,127 bp | **REAL** — NCBI assembly record |
| GC | 41.5% | **REAL** — NCBI |
| Scaffolds / contigs | 164 / 374 | **REAL** — NCBI |
| Scaffold N50 / L50 | 71,229,700 bp / 5 | **REAL** — NCBI |
| Contig N50 / L50 | 17,905,263 bp / 19 | **REAL** — NCBI |
| Spanned gaps / gap bases | 210 / 71,937 bp | **REAL** — Squires et al. 2023 |
| BUSCO aves_odb10 | C 98.9% (S 98.6%, D 0.28%), F 0.19%, M 0.91%, n 8338 | **REAL** — NCBI-computed |
| Merqury QV | Q43.2 | **REPRESENTATIVE** — QC'd with Merqury to the VGP standard (QV ≥ 40); no single numeric QV is published, so a value in the reported tier is used |
| Merqury k-mer completeness | 99.4% | **REPRESENTATIVE** — VGP-tier value |
| tidk telomere windows | 25 chromosome-scale scaffolds, 19 capped both ends | **REPRESENTATIVE** — models a high-quality-but-not-T2T VGP assembly; raw tidk output not distributed |

**Citation.** Squires TE, Rodin-Mörch P, Formenti G, et al. *A chromosome-level genome
assembly for the Rock Ptarmigan (Lagopus muta).* G3 (Bethesda). 2023;13(7):jkad099.
doi:10.1093/g3journal/jkad099 (PMC10320755). NCBI:
<https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_023343835.1/>. VGP BioProject
PRJNA836583.

---

## 2. Philippine tarsier — *Carlito syrichta* — `GCA_000164805.2` (Tarsius_syrichta-2.0.1)

| Metric | Value | Provenance |
| --- | --- | --- |
| Total length | 3,453,847,770 bp | **REAL** — NCBI assembly record |
| GC | 41.0% | **REAL** — NCBI |
| Scaffolds / contigs | 337,188 / 492,902 | **REAL** — NCBI |
| Scaffold N50 / L50 | 401,181 bp / 2,384 | **REAL** — NCBI |
| Contig N50 / L50 | 38,165 bp / 23,500 | **REAL** — NCBI |
| Gaps in scaffolds | 155,714 (= contigs − scaffolds) | **DERIVED** from REAL counts |
| Gap bases | ~43.8 Mb (total − ungapped 3.41 Gb) | **DERIVED** from REAL lengths |
| BUSCO primates_odb10 | C 86.2% (S 84.97%, D 1.2%), F 7.5%, M 6.4%, n 13780 | **REAL** — NCBI-computed |
| Merqury | *omitted* | No QV — this 2013 Sanger/454/Illumina (Celera) draft predates HiFi/Merqury QV; demonstrates graceful degradation |
| tidk / telomeres | *omitted* | No chromosome-level structure to evaluate |

**Citation.** Schmitz J, Noll A, Raabe CA, et al. *Genome sequence of the basal
haplorrhine primate Tarsius syrichta reveals unusual insertions.* Nat Commun.
2016;7:12997. doi:10.1038/ncomms12997 (PMC5059674). Assembly by Washington University,
2013. NCBI: <https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000164805.2/>. BioProject
PRJNA20339.

---

The BUSCO integer counts are rounded from the published percentages × *n* and may
differ by ±1 from a re-run; the percentages match the NCBI-computed values. The
Rock Ptarmigan's representative QV/telomere inputs are clearly labelled in
[`generate.py`](generate.py) and never overstate the published claims (the report's
single NOTE correctly says the assembly is *not* yet T2T).
