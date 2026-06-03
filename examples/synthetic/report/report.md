# AsmQC report — example_bird_v1

**Overall status:** FLAG  

_Generated 2026-06-03 06:42:21 UTC by AsmQC 0.1.0_


**Flag summary:** PASS: 1 · NOTE: 1 · WARN: 2 · FLAG: 1


## Curation flags

### [FLAG] FCS-GX flagged contamination

- **Observation:** FCS-GX flagged 1 region(s) across 1 sequence(s), 6,000 bp (6.0 kb) for removal/trim; top taxa: Bradyrhizobium sp..
- **Likely cause:** NCBI FCS-GX flagged sequence whose best taxonomic match is a foreign organism (bacteria, other host, adaptor, organelle in the wrong place). EXCLUDE/TRIM actions are high-confidence foreign sequence that NCBI would remove before release; it often hides in small, low-coverage contigs.
- **Suggested action:** Remove EXCLUDE sequences and TRIM the flagged spans (the GFF/BED export marks exact coordinates). Re-run FCS after cleaning and confirm no on-target sequence was lost. REVIEW calls need a manual look.

### [WARN] Elevated BUSCO duplication

- **Observation:** BUSCO duplicated is 6.3% (aves_odb10), above the 5% warning threshold — possible haplotype duplication.
- **Likely cause:** BUSCO looks for genes that should occur exactly once in this lineage. When many of them appear twice, the assembly is probably carrying both haplotypes of a diploid genome in the same (supposedly haploid) set of sequences — i.e. uncollapsed 'haplotigs' / false duplications. A truly duplicated region (recent whole-genome duplication) is possible but rarer.
- **Suggested action:** Run purge_dups (or purge_haplotigs) to collapse retained haplotigs into a primary haplotype, then re-check BUSCO duplication and the Merqury copy-number spectrum (the read-only/2-copy peaks). If duplication is biological, document the expected ploidy/WGD.

### [WARN] Consensus accuracy below standard

- **Observation:** Merqury QV is 38.6, below the Q40 EBP standard (Q40 = 99.99% / 1 error per 10 kb).
- **Likely cause:** Merqury's QV is a reference-free estimate of per-base consensus accuracy from k-mers that appear in the assembly but not the reads (assumed errors). A low QV means residual base-level errors — common with noisy long reads that were not polished, or under-polished consensus. Q40 ≈ 1 error / 10 kb (the EBP standard); below Q30 most genes will carry a base error.
- **Suggested action:** Polish the consensus (e.g. with the HiFi reads via DeepVariant/racon/medaka as appropriate) and re-run Merqury. Check the k-mer spectrum for a residual error peak at low multiplicity.

### [NOTE] Not telomere-to-telomere

- **Observation:** 3/5 chromosome-scale scaffolds are telomere-capped at BOTH ends (60%) (motif TTAGGG); 90% needed for a T2T claim.
- **Likely cause:** Telomeric repeat arrays were not detected at both ends of most chromosome-scale scaffolds. That means the assembly is not yet telomere-to-telomere: some chromosome ends are missing or unresolved (telomeres are highly repetitive and hard to assemble). This is expected for most non-T2T assemblies and is a 'note', not a failure.
- **Suggested action:** If T2T completeness is a goal, target the missing ends with ultra-long reads and verify with the telomere map. Otherwise record which chromosomes are arm-to-arm complete.

### [PASS] Scaffold N50 OK

- **Observation:** Scaffold N50 is 100.7 kb (>= 80.0 kb).
- **Likely cause:** This metric meets the configured quality standard.
- **Suggested action:** No action required.


## Metrics

| Group | Metric | Value |
| --- | --- | --- |
| Contiguity | Total length | 438.3 kb |
| Contiguity | Sequences (scaffolds) | 26 |
| Contiguity | Contigs | 28 |
| Contiguity | Scaffold N50 | 100.7 kb |
| Contiguity | Scaffold N90 | 40.0 kb |
| Contiguity | L50 | 2 |
| Contiguity | Contig N50 | 60.0 kb |
| Contiguity | Largest | 120.0 kb |
| Contiguity | GC content | 42.4% |
| Contiguity | Gaps | 2 |
| Contiguity | Gap bases | 700 bp |
| Contiguity | N per 100 kbp | 159.70 |
| Completeness (BUSCO) | Complete | 96.4% |
| Completeness (BUSCO) | Single-copy | 90.1% |
| Completeness (BUSCO) | Duplicated | 6.3% |
| Completeness (BUSCO) | Fragmented | 1.6% |
| Completeness (BUSCO) | Missing | 2.0% |
| Completeness (BUSCO) | Lineage | aves_odb10 |
| Correctness (Merqury) | Consensus QV | 38.60 |
| Correctness (Merqury) | k-mer completeness | 96.1% |
| Telomeres | Scaffolds w/ both ends | 3 |
| Telomeres | Scaffolds w/ any end | 4 |
| Telomeres | Repeat motif | TTAGGG |
| Contamination | Source | FCS |
| Contamination | Sequences flagged | 1 |
| Contamination | Flagged bases | 6.0 kb |

## Inputs

- **busco:** examples/synthetic/inputs/short_summary.txt
- **merqury:** examples/synthetic/inputs/merqury
- **tidk:** examples/synthetic/inputs/asm_telomeric_repeat_windows.tsv
- **contamination:** examples/synthetic/inputs/assembly.fcs_gx_report.txt
- **fasta:** examples/synthetic/inputs/assembly.fasta

---

Thresholds follow the Earth BioGenome Project / VGP 3C framework (Contiguity, Completeness, Correctness). See the AsmQC README for what each metric means and what good vs bad looks like.

