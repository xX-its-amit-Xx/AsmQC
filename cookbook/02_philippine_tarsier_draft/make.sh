#!/usr/bin/env bash
# Build the Philippine tarsier (GCA_000164805.2) AsmQC report.
# Inputs are reconstructed from published metrics by cookbook/generate.py
# (see cookbook/SOURCES.md). Run from the repository root.
set -euo pipefail
HERE="cookbook/02_philippine_tarsier_draft"

asmqc run \
  --gfastats "$HERE/inputs/gfastats_summary.txt" \
  --busco    "$HERE/inputs/short_summary.specific.primates_odb10.txt" \
  --name "Carlito_syrichta_Tarsius_2.0.1 (GCA_000164805.2)" \
  --out "$HERE/report"

echo "Open $HERE/report/report.html"
