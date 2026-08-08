import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

config = ["Full", "Skip 5m", "Slice in y"]

# remaining events without SBT Veto
LL_He = np.array([173.8, 173.8, 173.8])
LL_SBT = np.array([20997.5, 20997.5])
LL_SBT_f =  [843.27]


LX_He = np.array([1123.43, 1123.43, 1123.43])
LX_SBT_f = np.array([5902.19, 4842.65, 807.09])


# remaining events when SBT veto 90 MeV is applied
LL_He_Veto = np.array([173.8, 173.8, 173.8])
LL_SBT_Veto_f = np.array([20.34, 17.78, 1.3])

LX_He_Veto = np.array([1120.66, 1120.99, 1120.66])
LX_SBT_Veto_f = np.array([5.72, 4.1, 1.24])

# dataviz reference palette, same He/SBT colors as plot_cut_efficiency.py
COLOR_HE = "#f6c5db"
COLOR_SBT = "#b8c6de"
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
BASELINE = "#c3c2b7"
FACTORIZED_EDGE = "#b3401f"  # warning-orange ring on factorized points, regardless of He/SBT color

DIRECT_MARKER = 'o'
FACTORIZED_MARKER = 'D'
DODGE = 0.08  # horizontal offset so a no-veto/veto pair stays visible even when the values coincide

def direct_series(values):
    """All points in `values` are direct MC counts."""
    values = np.asarray(values, dtype=float)
    return values, np.zeros(len(values), dtype=bool)


def factorized_series(values):
    """All points in `values` are factorized (backfilled) estimates."""
    values = np.asarray(values, dtype=float)
    return values, np.ones(len(values), dtype=bool)


def combine(direct, factorized):
    """Concatenate a direct-count array with a factorized array into one
    per-config series, remembering which points are factorized so each point
    (not the whole series) can be marked correctly."""
    direct = np.asarray(direct, dtype=float)
    factorized = np.asarray(factorized, dtype=float)
    values = np.concatenate([direct, factorized])
    mask = np.concatenate([np.zeros(len(direct), dtype=bool), np.ones(len(factorized), dtype=bool)])
    assert len(values) == len(config), "series length must match config"
    return values, mask


# one panel per (channel, region): no-veto series vs. with-veto series.
# each series is (values, factorized_mask) -- factorized_mask flags, per point
# (not per whole array), which values are backfilled rather than directly counted.
# "ylim_split" is a list of (lo, hi) tiers, top to bottom, for panels whose values
# span multiple orders of magnitude. LL SBT needs 3 tiers (~21000 / ~843 / ~1-20):
# with only 2 tiers the 843 point sits so close to the tier boundary that its
# marker visually clips/bleeds into the break line.
panels = [
    dict(title="LL He", no_veto=direct_series(LL_He),
         veto=direct_series(LL_He_Veto), color=COLOR_HE, ylim_split=None),
    dict(title="LL SBT", no_veto=combine(LL_SBT, LL_SBT_f),
         veto=factorized_series(LL_SBT_Veto_f), color=COLOR_SBT,
         ylim_split=[(15000, 22000), (700, 1000), (-1, 25)]),
    dict(title="LX He", no_veto=direct_series(LX_He),
         veto=direct_series(LX_He_Veto), color=COLOR_HE, ylim_split=None),
    dict(title="LX SBT", no_veto=factorized_series(LX_SBT_f),
         veto=factorized_series(LX_SBT_Veto_f), color=COLOR_SBT,
         ylim_split=[(600, 6300), (-0.3, 7)]),
]


def draw_pair(ax, no_veto, veto, color):
    """Draw the no-veto (solid) and with-veto (open ring) series onto ax.
    Each series is split into its direct and factorized points, so a single
    series can mix both (e.g. two direct configs + one backfilled config)."""
    no_veto_vals, no_veto_mask = no_veto
    veto_vals, veto_mask = veto
    x = np.arange(len(config))

    direct = ~no_veto_mask
    if direct.any():
        ax.scatter(
            x[direct] - DODGE, no_veto_vals[direct], label="no veto",
            marker=DIRECT_MARKER, s=90, facecolor=color, edgecolor=INK_PRIMARY,
            linewidth=0.8, zorder=3,
        )
    if no_veto_mask.any():
        ax.scatter(
            x[no_veto_mask] - DODGE, no_veto_vals[no_veto_mask], label="no veto (factorized)",
            marker=FACTORIZED_MARKER, s=90, facecolor=color, edgecolor=FACTORIZED_EDGE,
            linewidth=2.4, zorder=3,
        )

    direct = ~veto_mask
    if direct.any():
        ax.scatter(
            x[direct] + DODGE, veto_vals[direct], label="SBT veto 90 MeV",
            marker=DIRECT_MARKER, s=150, facecolor='none', edgecolor=color,
            linewidth=2.2, zorder=4,
        )
    if veto_mask.any():
        # factorized veto points additionally get a dashed ring to flag they are backfilled, not directly counted
        ax.scatter(
            x[veto_mask] + DODGE, veto_vals[veto_mask], label="SBT veto 90 MeV (factorized)",
            marker=FACTORIZED_MARKER, s=150, facecolor='none', edgecolor=FACTORIZED_EDGE,
            linewidth=2.6, linestyle="--", zorder=4,
        )
    return x


def style_axis(ax, show_xticklabels, prune=None):
    ax.set_facecolor(SURFACE)
    ax.set_xticks(np.arange(len(config)))
    ax.set_xticklabels(config if show_xticklabels else [], color=INK_PRIMARY, fontsize=22)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune=prune))
    ax.tick_params(axis='y', colors=INK_MUTED, labelsize=20)
    ax.tick_params(axis='x', colors=INK_MUTED, length=0)
    for spine_name in ("top", "right"):
        ax.spines[spine_name].set_visible(False)
    ax.spines["left"].set_color(BASELINE)
    ax.spines["bottom"].set_color(BASELINE)


def add_break_marks(ax_top, ax_bottom):
    """Diagonal slashes at the split, marking the skipped middle range."""
    ax_top.spines["bottom"].set_visible(False)
    ax_bottom.spines["top"].set_visible(False)
    ax_top.tick_params(bottom=False, labelbottom=False)

    d = 0.6
    kwargs = dict(
        marker=[(-1, -d), (1, d)], markersize=14, linestyle="none",
        color=INK_MUTED, mec=INK_MUTED, mew=1.4, clip_on=False,
    )
    ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
    ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **kwargs)


# one shared legend for the whole figure: marker shape/fill encodes veto and
# factorized status the same way in every panel; only the fill color (He vs
# SBT) differs, and that's already named in each panel's title.
LEGEND_HANDLES = [
    Line2D([0], [0], marker=DIRECT_MARKER, linestyle='None', markersize=16,
           markerfacecolor=INK_MUTED, markeredgecolor=INK_PRIMARY, markeredgewidth=1.2,
           label="no veto"),
    Line2D([0], [0], marker=FACTORIZED_MARKER, linestyle='None', markersize=16,
           markerfacecolor=INK_MUTED, markeredgecolor=FACTORIZED_EDGE, markeredgewidth=2.6,
           label="no veto (factorized)"),
    Line2D([0], [0], marker=DIRECT_MARKER, linestyle='None', markersize=18,
           markerfacecolor='none', markeredgecolor=INK_MUTED, markeredgewidth=2.2,
           label="SBT veto 90 MeV"),
    Line2D([0], [0], marker=FACTORIZED_MARKER, linestyle='None', markersize=18,
           markerfacecolor='none', markeredgecolor=FACTORIZED_EDGE, markeredgewidth=2.8,
           label="SBT veto 90 MeV (factorized)"),
]

fig = plt.figure(figsize=(20, 13))
fig.patch.set_facecolor(SURFACE)
outer = fig.add_gridspec(
    2, 2, hspace=1.0, wspace=0.55,
    left=0.09, right=0.97, top=0.89, bottom=0.08,
)

for cell, panel in zip(outer, panels):
    tiers = panel["ylim_split"]
    if tiers is None:
        ax = fig.add_subplot(cell)
        draw_pair(ax, panel["no_veto"], panel["veto"], panel["color"])
        style_axis(ax, show_xticklabels=True)
        ax.set_title(panel["title"], color=INK_PRIMARY, fontsize=26)
        ax.set_ylabel("Remaining Events", color=INK_PRIMARY, fontsize=22, labelpad=18)
    else:
        n = len(tiers)
        inner = GridSpecFromSubplotSpec(n, 1, subplot_spec=cell, height_ratios=[1] * n, hspace=0.08)

        tier_axes = []
        for i, ylim in enumerate(tiers):
            shared = tier_axes[0] if tier_axes else None
            ax = fig.add_subplot(inner[i], sharex=shared)
            draw_pair(ax, panel["no_veto"], panel["veto"], panel["color"])
            ax.set_ylim(*ylim)
            # drop the tick right at each break so it can't collide with the
            # neighboring tier's boundary tick
            if n == 1:
                prune = None
            elif i == 0:
                prune = 'lower'
            elif i == n - 1:
                prune = 'upper'
            else:
                prune = 'both'
            style_axis(ax, show_xticklabels=(i == n - 1), prune=prune)
            tier_axes.append(ax)

        for ax_top, ax_bottom in zip(tier_axes, tier_axes[1:]):
            add_break_marks(ax_top, ax_bottom)

        tier_axes[0].set_title(panel["title"], color=INK_PRIMARY, fontsize=26)
        tier_axes[-1].set_ylabel("Remaining Events", color=INK_PRIMARY, fontsize=22)
        tier_axes[-1].yaxis.set_label_coords(-0.22, n / 2)

fig.suptitle("Muon Background - Remaining Events", color=INK_PRIMARY, fontsize=30, y=0.97)
fig.legend(
    handles=LEGEND_HANDLES, loc="center", bbox_to_anchor=(0.53, 0.485),
    frameon=True, fontsize=20, labelcolor=INK_SECONDARY,
    facecolor="#ececea", edgecolor=BASELINE, framealpha=0.95,
).get_frame().set_linewidth(1.2)
fig.savefig("MuonBackground_remaining_events.png", facecolor=SURFACE)
plt.show()
