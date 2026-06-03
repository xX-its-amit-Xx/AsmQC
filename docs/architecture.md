# Architecture

AsmQC is a small, layered pipeline. Data flows one way; every layer depends only on
the ones above it.

```
 inputs (files)                         CLI / Python API
      │                                   (cli.py, core.run_report)
      ▼                                          │
 parsers/*  ──►  models.*  ◄── fasta_stats   core.build_report
 (busco,         (typed                            │
  merqury,        dataclasses:           ┌──────────┼───────────┐
  quast,          AssemblyStats,         ▼          ▼           ▼
  gfastats,       BuscoResult,        flags.py    plots.py   report.py
  tidk,           MerquryResult, …)   (heuristics  (matplotlib (html/md/json)
  contamination)        │              + config)    base64)        │
                        │                  │                       ▼
                        └──────────────────┴──────────► exporters/ (gff3/bed/jbrowse2)
```

## Layers

1. **`models.py`** — the contract. Plain dataclasses (`AssemblyStats`, `BuscoResult`,
   `MerquryResult`, `TelomereResult`, `ContaminationResult`, `Flag`, `QCReport`) with
   `to_dict()` for JSON. Every other module agrees on these field names, so the JSON
   schema, the report, and the tests can't drift apart.

2. **`parsers/`** and **`fasta_stats.py`** — turn files into models. Each `parse_*`
   returns a model or `None`; none of them raise on malformed content. `fasta_stats`
   streams the FASTA (handles `.gz`) so it scales to multi-gigabase genomes.

3. **`config.py`** — loads `data/default_config.yaml` and deep-merges a user YAML over
   it. The defaults encode EBP/VGP thresholds and are fully commented.

4. **`flags.py`** + **`explanations.py`** — the curation heuristics. `flags.evaluate`
   reads the populated models + config and emits `Flag`s; `explanations` supplies the
   plain-English cause + action for each.

5. **`plots.py`** — matplotlib figures (Agg backend) returned as base64 PNGs. Every
   builder is defensive and returns `None` on failure, so a plotting hiccup never
   sinks the report.

6. **`report.py`** — renders a `QCReport` to a self-contained HTML page (inline CSS +
   base64 images), Markdown, and JSON.

7. **`exporters/`** — `tracks.py` (GFF3/BED of located flags; gaps recovered from the
   FASTA) and `jbrowse.py` (a JBrowse2 `config.json`).

8. **`core.py`** — orchestration: `build_report` parses everything and runs flags;
   `write_outputs` renders. `cli.py` is a thin argparse wrapper.

## Design choices

- **Degrade gracefully.** Any subset of inputs works; a missing/garbled input becomes
  a warning on the report, never an exception. `core._try` wraps every parser.
- **FASTA is the richest stats source.** When both a FASTA and a QUAST/gfastats report
  are supplied, AsmQC keeps the FASTA's per-sequence detail (needed for gap export,
  fragmentation, telomere sizing) and overrides the headline numbers from the report.
- **Stable, machine-readable JSON.** `summary.json` is the integration surface — its
  keys are fixed in `models.py`.
- **Few dependencies.** Core needs only PyYAML, matplotlib, numpy. pyfaidx/biopython
  are optional accelerators.

## Adding a parser

1. Add a model (or reuse one) in `models.py`.
2. Write `parsers/<tool>.py` with `parse_<tool>(path) -> Model | None`; export it from
   `parsers/__init__.py`.
3. Add a fixture under `tests/fixtures/` and a test.
4. Wire it into `core.build_report` and, if it produces flags, `flags.evaluate`.
