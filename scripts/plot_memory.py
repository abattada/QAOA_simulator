"""Plot the memory footprint comparison against QuEST as a dot plot.

Values are peak resident memory read from the system monitor while each
simulator ran the same grid search (screenshots in
benchmarks/memory_screenshots/). QAOA_simulator evaluates ~12 grid points
concurrently — one state vector per OpenMP thread — so its footprint is
roughly 12 state vectors versus QuEST's single one: memory traded for
wall-clock speed.

Usage:
    python3 scripts/plot_memory.py [-o docs/figures/memory_benchmark.png]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from plot_style import (C_MYSIM, C_QUEST, INK_MUTED, INK_PRIMARY,
                        INK_SECONDARY, style_axes, style_figure)

QUBITS = [18, 20, 22]
MYSIM_MB = [49.6, 201.2, 805.4]
QUEST_MB = [4.1, 33.0, 66.8]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default="docs/figures/memory_benchmark.png")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(8, 3.2))
    style_figure(fig)
    style_axes(ax, grid_axis="x")

    rows = range(len(QUBITS))
    ax.scatter(QUEST_MB, rows, s=70, color=C_QUEST, zorder=3, label="QuEST")
    ax.scatter(MYSIM_MB, rows, s=70, color=C_MYSIM, zorder=3,
               label="QAOA_simulator (12 concurrent grid points)")
    for y, (mine, quest) in zip(rows, zip(MYSIM_MB, QUEST_MB)):
        ax.text(mine, y + 0.22, f"{mine:g} MB", ha="center",
                color=INK_SECONDARY, fontsize=9)
        ax.text(quest, y + 0.22, f"{quest:g} MB", ha="center",
                color=INK_SECONDARY, fontsize=9)

    ax.set_xscale("log")
    ax.set_yticks(list(rows), [f"{q} qubits" for q in QUBITS])
    ax.set_ylim(-0.6, len(QUBITS) + 0.4)
    ax.set_xlabel("peak memory, MB (log scale)", color=INK_SECONDARY, fontsize=10)
    ax.set_title("Memory footprint during the grid search",
                 color=INK_PRIMARY, fontsize=12, loc="left", pad=10)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc="upper left", frameon=False,
              fontsize=9, labelcolor=INK_SECONDARY)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
