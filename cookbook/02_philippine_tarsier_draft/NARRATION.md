# Philippine tarsier (*Carlito syrichta*) — how the flags catch a draft

**Assembly:** `GCA_000164805.2` (Tarsius_syrichta-2.0.1), a 2013 Washington University
draft built with the Celera assembler from Sanger + 454 + Illumina reads — a
pre-long-read, pre-HiFi assembly.

**AsmQC overall status: `FLAG`.**

Run it yourself:

```bash
python cookbook/generate.py          # writes inputs/ and report/
# or directly:
asmqc run \
  --gfastats cookbook/02_philippine_tarsier_draft/inputs/gfastats_summary.txt \
  --busco    cookbook/02_philippine_tarsier_draft/inputs/short_summary.specific.primates_odb10.txt \
  --out cookbook/02_philippine_tarsier_draft/report
```

## What the flags said

| Flag | Severity | Why |
| --- | --- | --- |
| `busco_completeness` | **FLAG** | 86.2% complete — below the 90% floor; expected primate genes are missing or broken |
| `busco_fragmented` | **WARN** | 7.5% of BUSCOs are fragmented (split across contig ends) |
| `busco_missing` | **WARN** | 6.3% of BUSCOs are missing entirely |
| `contiguity_contig_n50` | **WARN** | Contig N50 = 38 kb — far below the 100 kb floor |
| `contiguity_scaffold_n50` | **WARN** | Scaffold N50 = 0.40 Mb — a fragmented draft, not chromosome-scale |
| `contiguity_many_sequences` | **NOTE** | 337,188 sequences |

## Reading it

Notice how the flags tell **one coherent story**, not six unrelated problems:

1. **Low contiguity is the root cause.** Short reads + the Celera assembler produced a
   38 kb contig N50 and a 0.40 Mb scaffold N50, scattered across **337,188** sequences
   (L50 = 2,384 — it takes ~2,384 scaffolds to cover half the genome). The
   `many_sequences` NOTE and both contiguity WARNs all point here.
2. **Fragmentation degrades the gene space.** When contigs are only ~38 kb, genes get
   cut at contig boundaries — so 7.5% of BUSCOs come back *fragmented* and another 6.3%
   *missing*, dragging Complete down to 86.2% and tripping the completeness **FLAG**.
   This is the causal link AsmQC's explanations spell out: *fragmented BUSCOs track low
   contiguity.* Improving contiguity is the lever that would move completeness.
3. **Graceful degradation.** No Merqury QV exists for this 2013 assembly (it predates
   reference-free k-mer QV), and there is no chromosome structure for a telomere
   analysis — so AsmQC simply omits the Correctness and telomere sections rather than
   inventing them. The report shows what *can* be assessed and says nothing about what
   can't.

**Take-away:** AsmQC doesn't just list metrics — by attaching a likely cause and action
to each flag, it lets you see that the contiguity and completeness flags are the same
problem viewed two ways, and that the fix (long-read re-assembly / scaffolding) targets
the root, not the symptoms.

> The original submission was screened for non-tarsier contaminants, but small contigs
> are exactly where residual contamination hides — running FCS-GX/Kraken2 and feeding
> the report to `asmqc run --contamination …` would add a contamination panel and (for
> FCS) a browser track of flagged regions. See [`examples/`](../../examples/) for that
> workflow on a dataset where contamination is present.
