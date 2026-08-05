#!/usr/bin/env python3
"""Plot He vs. SBT efficiency (%) of a single cut, per PID category.

Reads an `EffCuts.txt` produced by `selectionsteps.persist_cut_efficiencies()`
(one `tabulate` grid table per PID category: ee, mumu, emu, ex, mux, ll, lx,
all) and plots a chosen cut's He%/SBT% as grouped bars across categories.

Usage:
    python plot_cut_efficiency.py EffCuts.txt --cut "IP<10" [-o output.png]
    python plot_cut_efficiency.py EffCuts.txt --list-cuts
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# dataviz reference palette, categorical slots 1-2 (light mode)
COLOR_HE = "#f6c5db"   # blue
COLOR_SBT = "#b8c6de"  # orange
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
BASELINE = "#c3c2b7"

BLOCK_HEADER_RE = re.compile(r"^===\s*(.+?)\s*===$")


def parse_effcuts(path):
    """Return {category: {cut_name: (he_pct, sbt_pct)}}, categories in file order."""
    lines = Path(path).read_text().splitlines()

    blocks = []  # [(category, [lines...])]
    current = None
    for line in lines:
        m = BLOCK_HEADER_RE.match(line.strip())
        if m:
            current = (m.group(1), [])
            blocks.append(current)
        elif current is not None:
            current[1].append(line)

    result = {}
    for category, block_lines in blocks:
        rows = [l for l in block_lines if l.strip().startswith("|")]
        if not rows:
            continue
        header_cells = [c.strip() for c in rows[0].strip().strip("|").split("|")]
        cut_idx = header_cells.index("Cut")
        he_idx = header_cells.index("He %")
        sbt_idx = header_cells.index("SBT %")

        cuts = {}
        for row in rows[1:]:
            cells = [c.strip() for c in row.strip().strip("|").split("|")]
            if len(cells) <= max(cut_idx, he_idx, sbt_idx):
                continue
            cut_name = cells[cut_idx]
            try:
                he_pct = float(cells[he_idx])
                sbt_pct = float(cells[sbt_idx])
            except ValueError:
                continue
            cuts[cut_name] = (he_pct, sbt_pct)
        result[category] = cuts

    return result


def list_cuts(data):
    seen = []
    for cuts in data.values():
        for name in cuts:
            if name not in seen:
                seen.append(name)
    return seen


def sanitize_filename(name):
    return re.sub(r"[^A-Za-z0-9.+_-]+", "_", name)


def plot(data, cut, output):
    categories = []
    he_vals = []
    sbt_vals = []
    missing = []
    for category, cuts in data.items():
        if cut not in cuts:
            missing.append(category)
            continue
        he_pct, sbt_pct = cuts[cut]
        categories.append(category)
        he_vals.append(he_pct)
        sbt_vals.append(sbt_pct)

    if not categories:
        raise ValueError(f"cut '{cut}' not found in any category")
    if missing:
        print(f"Warning: cut '{cut}' missing in categories {missing}, skipped")

    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    x = list(range(len(categories)))
    width = 0.35
    he_bars = ax.bar(
        [xi - width / 2 for xi in x], he_vals, width=width,
        color=COLOR_HE, edgecolor=SURFACE, linewidth=1.0, label="He", zorder=3,
    )
    sbt_bars = ax.bar(
        [xi + width / 2 for xi in x], sbt_vals, width=width,
        color=COLOR_SBT, edgecolor=SURFACE, linewidth=1.0, label="SBT", zorder=3,
    )

    y_max = max(he_vals + sbt_vals + [0])
    headroom = max(y_max * 0.03, 1.0)
    for bars in (he_bars, sbt_bars):
        for b in bars:
            h = b.get_height()
            ax.text(
                b.get_x() + b.get_width() / 2, h + headroom, f"{h:.2f}",
                ha="center", va="bottom", fontsize=7.5, color=INK_SECONDARY,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(categories, color=INK_PRIMARY, fontsize=9)
    ax.set_ylabel("Efficiency (%)", color=INK_PRIMARY, fontsize=10)
    ax.set_ylim(0, y_max * 1.15)
    ax.set_title(f"Cut efficiency: {cut}", color=INK_PRIMARY, fontsize=12, pad=14)

    ax.set_facecolor(SURFACE)
    ax.tick_params(axis="y", colors=INK_MUTED, labelsize=9)
    ax.tick_params(axis="x", colors=INK_MUTED, length=0)
    for spine_name in ("top", "right", "left"):
        ax.spines[spine_name].set_visible(False)
    ax.spines["bottom"].set_color(BASELINE)

    ax.legend(loc="upper right", frameon=False, fontsize=9, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", help="EffCuts.txt to parse")
    parser.add_argument("--cut", default=None, help="exact cut name to plot (see --list-cuts)")
    parser.add_argument("--list-cuts", action="store_true", help="print available cut names and exit")
    parser.add_argument("-o", "--output", default=None, help="output PNG path")
    args = parser.parse_args()

    data = parse_effcuts(args.file)

    if args.list_cuts:
        for name in list_cuts(data):
            print(name)
        return

    if not args.cut:
        print("Error: --cut is required (or use --list-cuts to see available names)", file=sys.stderr)
        sys.exit(1)

    if args.cut not in list_cuts(data):
        print(f"Error: cut '{args.cut}' not found. Run with --list-cuts to see available names.", file=sys.stderr)
        sys.exit(1)

    input_path = Path(args.file)
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(f"eff_{sanitize_filename(args.cut)}.png")

    plot(data, args.cut, output_path)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
