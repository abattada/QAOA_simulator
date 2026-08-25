"""Shared chart style for all figures in this repo.

Palette validated for colour-vision deficiency (worst-pair CVD deltaE 24.7,
normal-vision deltaE 33.6, contrast >= 3:1 on the light surface).
"""
from matplotlib.colors import LinearSegmentedColormap

# Surfaces & ink
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Series colours (fixed assignment: this work = blue, QuEST = orange)
C_MYSIM = "#2a78d6"
C_QUEST = "#eb6834"

# Sequential blue ramp (light -> dark) for magnitude encodings (heatmaps)
SEQ_BLUE_STEPS = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE_STEPS)


def style_axes(ax, grid_axis="both"):
    """Recessive grid and axes: hairline grid, muted ticks, no top/right spines."""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, labelsize=9)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)


def style_figure(fig):
    fig.patch.set_facecolor(SURFACE)
