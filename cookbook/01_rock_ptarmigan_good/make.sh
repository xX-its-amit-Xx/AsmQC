#!/usr/bin/env bash
# Build the Rock Ptarmigan (GCA_023343835.1) AsmQC report.
# Inputs are reconstructed from published metrics by cookbook/generate.py
# (see cookbook/SOURCES.md). Run from the repository root.
set -euo pipefail
HERE="cookbook/01_rock_ptarmigan_good"

asmqc run \
  --gfastats "$HERE/inputs/gfastats_summary.txt" \
  --busco    "$HERE/inputs/short_summary.specific.aves_odb10.txt" \
  --merqury  "$HERE/inputs/merqury" \
  --tidk     "$HERE/inputs/bLagMut1_telomeric_repeat_windows.tsv" \
  --name "Lagopus_muta_bLagMut1 (GCA_023343835.1)" \
  --out "$HERE/report" \
  --tracks

echo "Open $HERE/report/report.html"
