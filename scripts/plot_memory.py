"""Plot the memory footprint comparison against QuEST as a dot plot.

QAOA_simulator values are peak resident set sizes measured by
scripts/run_benchmarks.sh (benchmarks/my_memory_record.txt): one state
vector plus the 2-byte-per-bitstring cut table. QuEST values were read from
the system monitor during an earlier run (screenshots in
benchmarks/memory_screenshots/quest_*.png).

Usage:
    python3 scripts/plot_memory.py [-o docs/figures/memory_benchmark.png]
"""
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt

from plot_style import (C_MYSIM, C_QUEST, INK_MUTED, INK_PRIMARY,
                        INK_SECONDARY, style_axes, style_figure)

QUBITS = [18, 20, 22]
QUEST_MB = [4.1, 33.0, 66.8]
MEMORY_RECORD = Path("benchmarks/my_memory_record.txt")


def read_mysim_mb():
    """{k: maxrss_mb} from the record; key k = graph qml_bi_k = 2k qubits."""
    entries = re.findall(r"(\d+):\s*\{\s*maxrss_mb:\s*([0-9.]+)\s*\}",
                         MEMORY_RECORD.read_text())
    by_qubits = {2 * int(k): float(mb) for k, mb in entries}
    return [by_qubits[q] for q in QUBITS]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default="docs/figures/memory_benchmark.png")
    args = parser.parse_args()
    mysim_mb = read_mysim_mb()

    fig, ax = plt.subplots(figsize=(8, 3.2))
    style_figure(fig)
    style_axes(ax, grid_axis="x")

    rows = range(len(QUBITS))
    ax.scatter(QUEST_MB, rows, s=70, color=C_QUEST, zorder=3, label="QuEST")
    ax.scatter(mysim_mb, rows, s=70, color=C_MYSIM, zorder=3,
               label="QAOA_simulator (single state + cut table)")
    for y, (mine, quest) in zip(rows, zip(mysim_mb, QUEST_MB)):
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
