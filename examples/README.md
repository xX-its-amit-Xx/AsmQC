# Example dataset

A small, **fully synthetic** assembly + QC inputs used to exercise every AsmQC code
path and to produce the figures in the main README. **This is not a real genome** —
for worked examples on real published assemblies, see [`../cookbook/`](../cookbook/).

## Regenerate

```bash
python examples/make_example.py
```

This writes:

- `synthetic/inputs/` — a 450 kB FASTA (5 chromosome-scale scaffolds with gaps + short
  debris + a contaminant contig), a BUSCO `short_summary.txt`, a Merqury directory, a
  tidk windows file, and an FCS-GX report.
- `synthetic/example_thresholds.yaml` — thresholds scaled down for the miniature genome.
- `synthetic/report/` — the generated `report.html`, `report.md`, `summary.json`, and
  the `flags.gff3` / `flags.bed` / `jbrowse2_config.json` tracks.
- `../docs/img/example_*.png` — the four figures used in the README.

## Run AsmQC on it directly

```bash
asmqc run \
  --fasta examples/synthetic/inputs/assembly.fasta \
  --busco examples/synthetic/inputs/short_summary.txt \
  --merqury examples/synthetic/inputs/merqury \
  --tidk examples/synthetic/inputs/asm_telomeric_repeat_windows.tsv \
  --contamination examples/synthetic/inputs/assembly.fcs_gx_report.txt \
  --config examples/synthetic/example_thresholds.yaml \
  --out examples/synthetic/report --jbrowse
```

The example is built to land on **overall status `FLAG`** with a spread of flags
(a contamination `FLAG`, a duplication/contiguity `WARN`, a telomere `NOTE`, and
passing checks) so the report shows the full range of severities.
