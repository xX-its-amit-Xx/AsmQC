# Rock Ptarmigan (*Lagopus muta*) — what passing QC looks like

**Assembly:** `GCA_023343835.1` (bLagMut1), a Vertebrate Genomes Project
chromosome-level reference assembled with PacBio HiFi + Hi-C (hifiasm → purge_dups →
SALSA), 2023.

**AsmQC overall status: `NOTE`** — passes all three C's; the only note is that it is
not yet telomere-to-telomere.

Run it yourself:

```bash
python cookbook/generate.py          # writes inputs/ and report/
# or directly:
asmqc run \
  --gfastats cookbook/01_rock_ptarmigan_good/inputs/gfastats_summary.txt \
  --busco    cookbook/01_rock_ptarmigan_good/inputs/short_summary.specific.aves_odb10.txt \
  --merqury  cookbook/01_rock_ptarmigan_good/inputs/merqury \
  --tidk     cookbook/01_rock_ptarmigan_good/inputs/bLagMut1_telomeric_repeat_windows.tsv \
  --out cookbook/01_rock_ptarmigan_good/report --tracks
```

## What the flags said

| Flag | Severity | Why |
| --- | --- | --- |
| `contiguity_scaffold_n50` | **PASS** | Scaffold N50 = 71.2 Mb, comfortably above the 10 Mb chromosome-scale target |
| `busco_ok` | **PASS** | 98.9% complete, only 0.28% duplicated — gene space is essentially complete and not over-duplicated |
| `merqury_qv` | **PASS** | QV 43.2 ≥ Q40 (the EBP/VGP base-accuracy standard) |
| `telomere_t2t` | **NOTE** | 19 / 25 chromosome-scale scaffolds are capped by telomeres at *both* ends (76%); ≥ 90% is needed to call it T2T |

## Reading it

This is the textbook profile of a modern reference genome, and AsmQC says so by *not*
raising any warnings:

- **Contiguity.** A 71 Mb scaffold N50 with L50 = 5 means a handful of large scaffolds
  carry half the genome — chromosome scale. The 17.9 Mb contig N50 shows the underlying
  HiFi contigs are long too (above the EBP `6` = 1 Mb bar by more than an order of
  magnitude). Only 210 gaps across a gigabase.
- **Completeness.** 98.9% of single-copy bird orthologs are present and intact, and
  duplication is 0.28% — i.e. purge_dups did its job and the assembly is a clean haploid
  representation. (Contrast this with the tarsier's 7.5% fragmented.)
- **Correctness.** Reference-free QV ≥ 40 means fewer than one consensus error per
  10 kb — most genes are error-free.
- **The one NOTE.** The telomere map shows most chromosomes capped at both ends, but not
  all — so the assembly is excellent but *not yet* telomere-to-telomere. AsmQC records
  this as a NOTE, not a warning: it is a goal for further finishing, not a defect. The
  `--tracks` run also writes `report/flags.gff3`/`flags.bed` marking the detected
  telomere arrays so a curator can see exactly which chromosome ends are capped.

**Take-away:** when every C is met, an AsmQC report is mostly green, and the remaining
notes point at the frontier (T2T finishing) rather than at problems to fix.
