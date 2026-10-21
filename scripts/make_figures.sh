#!/usr/bin/env bash
# Regenerate every figure used in the README from the shipped results and
# benchmark records. Needs matplotlib and numpy (see README, "Plots and CSV").
#
# Usage (from the repo root):
#     scripts/make_figures.sh            # uses .venv/bin/python
#     PYTHON=python3 scripts/make_figures.sh
set -euo pipefail

PY=${PYTHON:-.venv/bin/python}
RUN=${RUN:-1000_20_0}

for k in 4 5 6 7 8 9 10 11 12; do
    "$PY" scripts/plot_heatmap.py "results/qml_bi_$k/$RUN" -o "docs/figures/heatmaps/qml_bi_$k.png"
done
"$PY" scripts/plot_time.py
"$PY" scripts/plot_memory.py
