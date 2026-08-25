"""Plot the runtime benchmark against QuEST.

Both simulators ran the identical workload: a full 20x20 (gamma, beta)
grid search of the p=1 QAOA MaxCut circuit on the HamLib bipartite graphs,
measured with `time` on the same 12-thread machine. Raw records live in
benchmarks/my_time_record.txt and benchmarks/quest_time_record.txt.

Usage:
    python3 scripts/plot_time.py [-o docs/figures/runtime_benchmark.png]
"""
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt

from plot_style import (C_MYSIM, C_QUEST, INK_MUTED, INK_PRIMARY,
                        INK_SECONDARY, style_axes, style_figure)

# Measured records (benchmarks/*_time_record.txt); key k = graph qml_bi_k = 2k qubits
DATA = {
    "my_time": {
        2: {"real": "0m0.005s", "user": "0m0.015s"},
        3: {"real": "0m0.006s", "user": "0m0.024s"},
        4: {"real": "0m0.009s", "user": "0m0.068s"},
        5: {"real": "0m0.031s", "user": "0m0.309s"},
        6: {"real": "0m0.148s", "user": "0m1.710s"},
        7: {"real": "0m0.776s", "user": "0m9.128s"},
        8: {"real": "0m4.062s", "user": "0m48.059s"},
        9: {"real": "0m21.329s", "user": "4m12.455s"},
        10: {"real": "2m12.202s", "user": "25m43.121s"},
        11: {"real": "11m1.137s", "user": "129m8.744s"},
        12: {"real": "51m52.128s", "user": "607m29.901s"},
    },
    "quest_time": {
        2: {"real": "0m0.050s", "user": "0m0.452s"},
        3: {"real": "0m0.074s", "user": "0m0.836s"},
        4: {"real": "0m0.123s", "user": "0m1.410s"},
        5: {"real": "0m0.181s", "user": "0m2.096s"},
        6: {"real": "0m0.316s", "user": "0m3.710s"},
        7: {"real": "0m1.038s", "user": "0m11.895s"},
        8: {"real": "0m3.668s", "user": "0m41.874s"},
        9: {"real": "0m15.379s", "user": "3m0.428s"},
        10: {"real": "1m8.750s", "user": "13m36.383s"},
        11: {"real": "11m30.237s", "user": "135m21.026s"},
        12: {"real": "66m38.876s", "user": "709m23.835s"},
    },
}


def parse_time(time_str):
    """Convert a time string like '2m12.202s' into seconds."""
    match = re.match(r"(?:(\d+)m)?([\d.]+)s", time_str)
    minutes = int(match.group(1)) if match.group(1) else 0
    return minutes * 60 + float(match.group(2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default="docs/figures/runtime_benchmark.png")
    args = parser.parse_args()

    keys = sorted(DATA["my_time"])
    qubits = [2 * k for k in keys]
    series = {
        which: {kind: [parse_time(DATA[which][k][kind]) for k in keys]
                for kind in ("real", "user")}
        for which in ("my_time", "quest_time")
    }

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    style_figure(fig)

    panels = [("real", "Wall-clock time"), ("user", "CPU time")]
    for ax, (kind, title) in zip(axes, panels):
        style_axes(ax)
        ax.plot(qubits, series["my_time"][kind], color=C_MYSIM, linewidth=2,
                marker="o", markersize=6, label="QAOA_simulator (this work)")
        ax.plot(qubits, series["quest_time"][kind], color=C_QUEST, linewidth=2,
                marker="o", markersize=6, label="QuEST")
        ax.set_yscale("log")
        ax.set_xticks(qubits)
        ax.set_xlabel("qubits", color=INK_SECONDARY, fontsize=10)
        ax.set_ylabel("seconds (log scale)", color=INK_SECONDARY, fontsize=10)
        ax.set_title(title, color=INK_PRIMARY, fontsize=12, loc="left", pad=10)

        my_end = series["my_time"][kind][-1]
        quest_end = series["quest_time"][kind][-1]
        ax.text(0.97, 0.28, f"{quest_end / my_end:.2f}× faster\nat 24 qubits",
                transform=ax.transAxes, ha="right",
                color=INK_SECONDARY, fontsize=9)

    axes[0].legend(loc="upper left", frameon=False, fontsize=9,
                   labelcolor=INK_SECONDARY)
    fig.suptitle("p=1 QAOA MaxCut, full 20×20 (γ, β) grid search — 12 hardware threads",
                 color=INK_MUTED, fontsize=9, x=0.99, y=0.99, ha="right")
    fig.tight_layout()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
