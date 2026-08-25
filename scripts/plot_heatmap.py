"""Render a (gamma, beta) grid-search landscape as a heatmap.

Reads a raw simulator output file (results/<graph>/<shots>_<steps>_<run>),
pivots the grid into a matrix and draws the cut value found at every
(gamma, beta) sample, with the best parameter point ringed.

By default the image is written next to the result file
(results/qml_bi_10/1_20_0 -> results/qml_bi_10/1_20_0.png); pass -o to
regenerate the README gallery figures in docs/figures/heatmaps/.

Usage:
    python3 scripts/plot_heatmap.py results/qml_bi_10/1_20_0 [-o out.png]
"""
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_style import (INK_MUTED, INK_PRIMARY, INK_SECONDARY, SEQ_BLUE,
                        SURFACE, style_figure)

LINE_RE = re.compile(
    r"Gamma:\s*([0-9.eE+-]+),\s*Beta:\s*([0-9.eE+-]+),\s*Max Cut:\s*([0-9.eE+-]+)")

PI_TICKS = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
PI_LABELS = ["0", "π/4", "π/2", "3π/4", "π"]


def parse_result(path):
    text = Path(path).read_text()
    samples = [(float(g), float(b), float(c)) for g, b, c in LINE_RE.findall(text)]
    meta = {}
    for key, pattern in [
        ("qubits", r"QAOA \(p=1\) with (\d+) qubits"),
        ("edges", r"edges=(\d+)"),
        ("best_cut", r"Best found cut\s*=\s*([0-9.]+)"),
        ("best_gamma", r"Best gamma\s*=\s*([0-9.]+)"),
        ("best_beta", r"Best beta\s*=\s*([0-9.]+)"),
    ]:
        m = re.search(pattern, text)
        meta[key] = float(m.group(1)) if m else None
    m = re.search(r"Input graph:\s*(\S+)", text)
    meta["graph"] = m.group(1) if m else Path(path).parent.name
    return samples, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_file", help="raw simulator output file")
    parser.add_argument("-o", "--output", default=None,
                        help="output image path (default: <result_file>.png)")
    args = parser.parse_args()

    samples, meta = parse_result(args.result_file)
    gammas = sorted({s[0] for s in samples})
    betas = sorted({s[1] for s in samples})
    grid = np.full((len(betas), len(gammas)), np.nan)
    g_idx = {g: i for i, g in enumerate(gammas)}
    b_idx = {b: i for i, b in enumerate(betas)}
    for g, b, c in samples:
        grid[b_idx[b], g_idx[g]] = c

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    style_figure(fig)
    ax.set_facecolor(SURFACE)

    mesh = ax.pcolormesh(gammas, betas, grid, cmap=SEQ_BLUE,
                         edgecolors=SURFACE, linewidth=0.6, vmin=0)

    # Ring the best (gamma, beta) sample: white halo + ink ring, never colour alone
    ax.scatter([meta["best_gamma"]], [meta["best_beta"]], s=170,
               facecolors="none", edgecolors=SURFACE, linewidths=3.4, zorder=3)
    ax.scatter([meta["best_gamma"]], [meta["best_beta"]], s=170,
               facecolors="none", edgecolors=INK_PRIMARY, linewidths=1.6, zorder=4)

    ax.set_xticks(PI_TICKS, PI_LABELS)
    ax.set_yticks(PI_TICKS, PI_LABELS)
    ax.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, labelsize=9)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel("γ", color=INK_SECONDARY, fontsize=12)
    ax.set_ylabel("β", color=INK_SECONDARY, fontsize=12, rotation=0, labelpad=10)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("cut value", color=INK_SECONDARY, fontsize=10)
    cbar.ax.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, labelsize=9)
    cbar.outline.set_visible(False)

    n = int(meta["qubits"])
    ax.set_title(f"{meta['graph']} — {n} qubits, {int(meta['edges'])} edges",
                 color=INK_PRIMARY, fontsize=13, pad=28, loc="left")
    ax.text(0, 1.02, f"best cut {int(meta['best_cut'])} at "
            f"γ={meta['best_gamma']:.2f}, β={meta['best_beta']:.2f}"
            f"  (○)",
            transform=ax.transAxes, color=INK_MUTED, fontsize=9)

    out = args.output or f"{args.result_file}.png"
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
