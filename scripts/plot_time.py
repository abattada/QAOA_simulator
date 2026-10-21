"""Plot the runtime benchmark against QuEST.

Both simulators ran a full 20x20 (gamma, beta) grid search of p=1 QAOA
MaxCut on the HamLib bipartite graphs, timed with `time` on the same
machine. QuEST applied the circuit gate by gate (12 threads, recorded
earlier); this simulator uses its fused kernels (6 threads, one per core).
Raw records live in benchmarks/my_time_record.txt (written by
scripts/run_benchmarks.sh) and benchmarks/quest_time_record.txt.

Usage:
    python3 scripts/plot_time.py [-o docs/figures/runtime_benchmark.png]
"""
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt

from plot_style import (C_MYSIM, C_QUEST, INK_MUTED, INK_PRIMARY,
                        INK_SECONDARY, style_axes, style_figure)

RECORDS = {
    "my_time": Path("benchmarks/my_time_record.txt"),
    "quest_time": Path("benchmarks/quest_time_record.txt"),
}
# key k in the records = graph qml_bi_k = 2k qubits
ENTRY_RE = re.compile(r"(\d+):\s*\{\s*real:\s*([0-9ms.]+),\s*user:\s*([0-9ms.]+)\s*\}")


def parse_time(time_str):
    """Convert a time string like '2m12.202s' into seconds."""
    match = re.match(r"(?:(\d+)m)?([\d.]+)s", time_str)
    minutes = int(match.group(1)) if match.group(1) else 0
    return minutes * 60 + float(match.group(2))


def parse_record(path):
    """Read a {k: {real: XmY.YYYs, user: XmY.YYYs}} record into seconds."""
    return {int(k): {"real": parse_time(r), "user": parse_time(u)}
            for k, r, u in ENTRY_RE.findall(path.read_text())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default="docs/figures/runtime_benchmark.png")
    args = parser.parse_args()

    data = {which: parse_record(path) for which, path in RECORDS.items()}
    keys = sorted(set(data["my_time"]) & set(data["quest_time"]))
    qubits = [2 * k for k in keys]
    series = {
        which: {kind: [data[which][k][kind] for k in keys]
                for kind in ("real", "user")}
        for which in data
    }

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    style_figure(fig)

    panels = [("real", "Wall-clock time"), ("user", "CPU time")]
    for ax, (kind, title) in zip(axes, panels):
        style_axes(ax)
        ax.plot(qubits, series["my_time"][kind], color=C_MYSIM, linewidth=2,
                marker="o", markersize=6, label="QAOA_simulator (fused kernels, 6 threads)")
        ax.plot(qubits, series["quest_time"][kind], color=C_QUEST, linewidth=2,
                marker="o", markersize=6, label="QuEST (gate by gate, 12 threads, earlier record)")
        ax.set_yscale("log")
        ax.set_xticks(qubits)
        ax.set_xlabel("qubits", color=INK_SECONDARY, fontsize=10)
        ax.set_ylabel("seconds (log scale)", color=INK_SECONDARY, fontsize=10)
        ax.set_title(title, color=INK_PRIMARY, fontsize=12, loc="left", pad=10)

        my_end = series["my_time"][kind][-1]
        quest_end = series["quest_time"][kind][-1]
        ratio = quest_end / my_end
        verdict = f"{ratio:.0f}× faster" if ratio >= 1 else f"{1 / ratio:.1f}× slower"
        ax.text(0.97, 0.28, f"{verdict}\nat {qubits[-1]} qubits",
                transform=ax.transAxes, ha="right",
                color=INK_SECONDARY, fontsize=9)

    axes[0].legend(loc="upper left", frameon=False, fontsize=9,
                   labelcolor=INK_SECONDARY)
    fig.suptitle("p=1 QAOA MaxCut, full 20×20 (γ, β) grid search, Ryzen 5 3600X",
                 color=INK_MUTED, fontsize=9, x=0.99, y=0.99, ha="right")
    fig.tight_layout()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
