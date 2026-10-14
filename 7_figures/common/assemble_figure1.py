#!/usr/bin/env python3
"""Assemble two audited source plots as the manuscript's Figure 1."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from audit_panel_alignment import require_matplotlib_panel_alignment


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-a", type=Path, required=True)
    parser.add_argument("--panel-b", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args()

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 9,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })

    images = [plt.imread(args.panel_a), plt.imread(args.panel_b)]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 9.45))
    fig.subplots_adjust(left=.018, right=.982, bottom=.012, top=.988, hspace=.025)
    for axis, image, letter in zip(axes, images, ("a", "b")):
        axis.imshow(image)
        axis.set_axis_off()
        axis.text(.005, .995, letter, transform=axis.transAxes, fontsize=11,
                  fontweight="bold", ha="left", va="top", clip_on=False)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=("a", "b"),
        column_groups=(("a", "b"),),
        json_out=Path(str(args.out_prefix) + ".alignment.json"),
        overlay_svg=Path(str(args.out_prefix) + ".alignment.svg"),
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )
    fig.savefig(Path(str(args.out_prefix) + ".png"), dpi=600)
    fig.savefig(Path(str(args.out_prefix) + ".pdf"))
    fig.savefig(Path(str(args.out_prefix) + ".svg"))
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
