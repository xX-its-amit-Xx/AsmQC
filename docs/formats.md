# Input formats AsmQC parses

AsmQC reads each tool's **native** output. The exact layouts below were verified
against each tool's source code and real example files. Parsers are defensive: an
unrecognised or partial file yields a parse warning, not a crash.

---

## BUSCO — `--busco`

Accepts either the text or JSON `short_summary`. Detected by content (`{` → JSON).

**`short_summary.specific.<lineage>.<label>.txt`** — plain text, TAB-indented. AsmQC reads:

- the one-line notation (order is fixed): `C:<C>%[S:<S>%,D:<D>%],F:<F>%,M:<M>%,n:<n>`
- the six labelled count lines, e.g. `7338\tComplete BUSCOs (C)` …
  `Complete and single-copy BUSCOs (S)`, `Complete and duplicated BUSCOs (D)`,
  `Fragmented BUSCOs (F)`, `Missing BUSCOs (M)`, `Total BUSCO groups searched`
- the header comments for version, `lineage dataset is: <name>`, and `run in mode:`

**`short_summary.*.json`** — reads `results.{Complete percentage, Complete BUSCOs,
Single copy percentage, Single copy BUSCOs, Multi copy percentage, Multi copy BUSCOs,
Fragmented percentage, Fragmented BUSCOs, Missing percentage, Missing BUSCOs,
n_markers}` and `lineage_dataset.name`. Note the JSON spells duplicated as
**"Multi copy"** (with a space), which AsmQC maps to `duplicated`.

If only counts are present, percentages are derived from them.

---

## Merqury — `--merqury`

Point at a **directory**, a **prefix**, or any one of the three files (siblings are
auto-discovered). All are TAB-delimited.

| File | Header? | Columns AsmQC reads |
| --- | --- | --- |
| `<p>.qv` | no | `name`, `asm-only k-mers`, `total k-mers`, **`QV`**, `error rate` |
| `<p>.completeness.stats` | no | `name`, `set(all)`, `found`, `total`, **`completeness %`** |
| `<p>.spectra-cn.hist` | yes | `Copies`, `kmer_multiplicity`, `Count` |

In multi-assembly `.qv` files a `Both` row, if present, is used for the headline QV.
`spectra-cn.hist` `Copies` categories are `read-only` (in reads, missing from
assembly), `1`/`2`/`3`/`4` (copies in the assembly), and `>4`. A `spectra-asm.hist`
(header column `Assembly`) is parsed the same way if that's all you have.

---

## QUAST — `--quast`  ·  gfastats — `--gfastats`

Both produce an assembly-stats record. If a FASTA is **also** given, AsmQC keeps the
FASTA's per-sequence richness (for gap export, fragmentation, telomere sizing) and
overrides the headline numbers (N50/N90/L50/total/GC/…) from the report.

**QUAST `report.tsv`** — TAB, `<metric name>\t<value>`. AsmQC reads the verbatim
metric strings `Total length`, `# contigs`, `# contigs (>= 0 bp)`, `Largest contig`,
`N50`, `N90`, `L50`, `L90`, `auN`, `GC (%)`, `# N's per 100 kbp`. (Multi-assembly
reports have one value column per assembly; AsmQC takes the first.)

**gfastats** default summary — `<label>: <value>` (colon-space; `--tabular` uses a
TAB). AsmQC reads `# scaffolds`, `Total scaffold length`, `Scaffold N50/L50`,
`Largest/Smallest scaffold`, `Average scaffold length`, `# contigs`, `Contig N50`,
`# gaps in scaffolds`, `Total gap length in scaffolds`, `GC content %`. The first two
stdout lines (input path, `embedded`) and `+++Assembly summary+++` are skipped; lines
are split on the **last** `: ` so `Base composition (A:C:G:T)` parses correctly.

---

## tidk — `--tidk`

**`*_telomeric_repeat_windows.tsv`** — TAB, with header
`id  window  forward_repeat_number  reverse_repeat_number  telomeric_repeat`. The
`.bedgraph` variant (no header, `id start end count`) is also accepted.

`window` is the window **end coordinate** (`window_size`, 2×, …, last clamped to the
scaffold length), so `max(window)` per id ≈ scaffold length. AsmQC infers a telomere
at a scaffold **end** when the motif count in the terminal window(s) reaches
`tidk_min_window_count` (default 50). By convention the start telomere shows as the
reverse-complement motif and the end telomere as the forward motif; AsmQC takes the
max of the two to stay orientation-robust.

---

## Contamination — `--contamination` (FCS-GX or Kraken2, auto-detected)

**FCS-GX `*.fcs_gx_report.txt`** — line 1 is a `##[[…]]` JSON header, line 2 is
`#seq_id\t…`, then 8 TAB columns: `seq_id, start_pos, end_pos, seq_len, action, div,
agg_cont_cov, top_tax_name`. `action` ∈ `EXCLUDE/TRIM/FIX/REVIEW/REVIEW_RARE/INFO`.
Coordinates are 1-based inclusive and feed the GFF/BED export.

**Kraken2 standard report** — TAB, no header, 6 columns: `percent, clade_reads,
taxon_reads, rank_code, taxid, name` (name indented by tree depth). The 8-column
`--report-minimizer-data` variant is handled. AsmQC summarises superkingdom (rank
`D`) percentages and the unclassified fraction; non-target superkingdoms above
`kraken_taxon_flag_pct` raise a flag.

---

## Sources

- BUSCO user guide & source (gitlab.com/ezlab/busco); Manni et al. 2021.
- Merqury wiki & source (github.com/marbl/merqury, `eval/qv.sh`, `eval/spectra-cn.sh`); Rhie et al. 2020.
- QUAST manual & `quast_libs/reporting.py`; gfastats `vgl-hub/gfalibs` `src/output.cpp`.
- tidk source (github.com/tolkit/telomeric-identifier, `src/search.rs`/`src/finder.rs`).
- NCBI FCS-GX `scripts/action_report.py`; Kraken2 `docs/MANUAL.markdown`.
