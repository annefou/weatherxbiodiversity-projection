"""Concordance heatmap for the four-variant cross-substrate comparison.

Reads results/variant_comparison_<horizon>.json and produces:
  - figures/variant_concordance_<horizon>.png — Spearman ρ heatmap, variant × n_cells filter
  - figures/variant_pairs_<horizon>.png       — per-variant scatter of η_64 vs η_128

Usage:
    python scripts/plot_variant_concordance.py --horizon 2030_2039
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path("/Users/annef/Documents/ScienceLive/weatherxbiodiversity-projection")


def main(horizon: str):
    summary = json.loads((REPO / f"results/variant_comparison_{horizon}.json").read_text())
    df = pd.read_csv(REPO / f"results/variant_comparison_{horizon}.csv")

    plt.style.use("seaborn-v0_8-whitegrid")

    variants = ["a", "b", "c", "d", "d2"]
    var_labels = [
        "(a) Full GLMM η\n(within-substrate)",
        "(b) Main-effects η\n(within-substrate)",
        "(c) Full η,\nshared CEA ref",
        "(d) Mean future TEI\n(substrate-invariant)",
        "(d2) Frac TEI_fut>0.5\n(substrate-invariant)",
    ]
    n_filters = list(summary["by_min_n"].keys())

    rho_grid = np.full((len(variants), len(n_filters)), np.nan)
    n_species_per_filter = []
    for j, n in enumerate(n_filters):
        n_species_per_filter.append(summary["by_min_n"][n]["n_species"])
        for i, v in enumerate(variants):
            r = summary["by_min_n"][n]["variants"][v]["spearman_rho"]
            if r is not None:
                rho_grid[i, j] = r

    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=150)
    im = ax.imshow(rho_grid, vmin=0.0, vmax=1.0, cmap="RdYlGn", aspect="auto")
    ax.set_xticks(range(len(n_filters)))
    ax.set_xticklabels(
        [f"≥{n}\n(n={ns})" for n, ns in zip(n_filters, n_species_per_filter)]
    )
    ax.set_yticks(range(len(variants)))
    ax.set_yticklabels(var_labels)
    ax.set_xlabel("Per-species cell-count filter (n_cells at both substrates)")
    ax.set_title(
        f"Cross-substrate concordance (Spearman ρ): nside=64 vs nside=128\n"
        f"horizon {horizon.replace('_', '–')}"
    )
    for i in range(len(variants)):
        for j in range(len(n_filters)):
            v = rho_grid[i, j]
            if not np.isnan(v):
                color = "white" if v < 0.5 else "black"
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", color=color, fontsize=9)
    cbar = fig.colorbar(im, ax=ax, label="Spearman ρ (rank concordance)")
    fig.tight_layout()
    out1 = REPO / "figures" / f"variant_concordance_{horizon}.png"
    out1.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out1)
    plt.close(fig)
    print(f"  wrote {out1}")

    # Per-variant scatter — colour by n_cells using the analytical thresholds
    # we actually report on (the same filters used in the concordance heatmap).
    # Discrete bins so the visual signal is "red = unreliable, green = trust".
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), dpi=150)
    axes = axes.flatten()
    pairs = [
        ("a_eta_64", "a_eta_128", "(a) Full GLMM η"),
        ("b_eta_64", "b_eta_128", "(b) Main-effects only η"),
        ("c_eta_64", "c_eta_128", "(c) Shared CEA ref η"),
        ("d_frac_64", "d_frac_128", "(d) Mean future TEI"),
        ("d2_frac_64", "d2_frac_128", "(d2) Frac TEI>0.5"),
    ]
    n_min = np.minimum(df["n_cells_64"], df["n_cells_128"]).astype(float)
    n_top = float(max(n_min.max(), 20)) + 1
    bin_edges = [1, 5, 10, 20, n_top]
    bin_colors = ["#d73027", "#fdae61", "#a6d96a", "#1a9850"]
    bin_labels = ["n<5\n(unreliable)", "n=5–9\n(borderline)", "n=10–19\n(reliable)", "n≥20\n(very reliable)"]
    cmap = plt.matplotlib.colors.ListedColormap(bin_colors)
    norm = plt.matplotlib.colors.BoundaryNorm(bin_edges, cmap.N)
    for ax, (xc, yc, title) in zip(axes, pairs):
        sc = ax.scatter(df[xc], df[yc], c=n_min, cmap=cmap, norm=norm,
                        s=55, alpha=0.9, edgecolor="k", linewidth=0.5)
        lo = float(np.nanmin([df[xc].min(), df[yc].min()]))
        hi = float(np.nanmax([df[xc].max(), df[yc].max()]))
        ax.plot([lo, hi], [lo, hi], "k--", alpha=0.4, linewidth=1)
        ax.set_xlabel("nside=64")
        ax.set_ylabel("nside=128")
        ax.set_title(title)
    axes[-1].axis("off")
    cbar = fig.colorbar(
        sc, ax=axes[-1],
        ticks=[3, 7, 14.5, (20 + n_top) / 2],  # midpoint of each bin
        shrink=0.85, fraction=0.5,
    )
    cbar.ax.set_yticklabels(bin_labels, fontsize=8)
    cbar.set_label("min(n_cells_64, n_cells_128)", fontsize=9)
    fig.suptitle(
        f"Per-species substrate concordance, horizon {horizon.replace('_', '–')} — "
        f"each point is one Bombus species; diagonal = perfect agreement",
        fontsize=11,
    )
    fig.tight_layout()
    out2 = REPO / "figures" / f"variant_pairs_{horizon}.png"
    fig.savefig(out2)
    plt.close(fig)
    print(f"  wrote {out2}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--horizon", default="2030_2039", choices=("2020_2029", "2030_2039"))
    args = ap.parse_args()
    main(args.horizon)
