# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 08 — Tier-2 projection figures (HEALPix nside=64)
#
# Produces:
#
#   * `figures/projection_species_rank.png` — two-panel ranked bar chart
#     (near-term left, mid-term right) of per-species posterior-mean
#     extirpation probability with 95 % HDI error bars; top-3 species
#     per horizon highlighted in gold.
#   * `figures/projection_risk_map_2020_2029.png` — Iberia HEALPix
#     nside=64 risk map of community-mean (over species) extirpation
#     probability (near-term).
#   * `figures/projection_risk_map_2030_2039.png` — same, mid-term.
#   * `figures/projection_summary.png` — combined panel for the
#     Jupyter Book / nanopub Outcome draft (rank chart + the more
#     impactful map).
#
# Per `DOMAIN.md`: HEALPix is always NESTED, and we use **healpix-geo**
# (NOT healpy) for the cell→lat/lon mapping, since healpy is not
# geo-aware (cosmology-first, no CRS handling) and accumulates small
# biases over decadal timeseries.

# %%
import json
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from healpix_geo.nested import vertices as hp_vertices
from matplotlib.collections import PolyCollection

# %%
plt.style.use("seaborn-v0_8-whitegrid")

ROOT = Path("..").resolve()
HPORT = ROOT / "healpix_port"
OUT_DIR = HPORT / "outputs_iberia"
RESULTS_DIR = ROOT / "results"
PRECOMP = ROOT / "data" / "precomputed"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Iberia nside=64 cells (depth=6, NESTED).
DEPTH = 6
IBERIA_PIX_64 = np.load(PRECOMP / "iberia_pix_nside64_nested.npy").astype(np.uint64)
N_64 = len(IBERIA_PIX_64)
print(f"Iberia HEALPix nside=64 cells: {N_64}")

# %% [markdown]
# ## Load projection summary + per-cell rasters

# %%
with open(RESULTS_DIR / "projection_headline.json") as f:
    summary = json.load(f)

HORIZONS = ["2020_2029", "2030_2039"]
HORIZON_TITLES = {
    "2020_2029": "Near-term (2020-2029)",
    "2030_2039": "Mid-term (2030-2039)",
}

per_cell = {}
for h in HORIZONS:
    p = RESULTS_DIR / f"projection_{h}.nc"
    if p.exists():
        ds_p = xr.open_dataset(p)
        per_cell[h] = ds_p["community_mean_p_extirpation"].values.astype(float)
        print(f"  loaded {p.name}: shape {per_cell[h].shape}")
    else:
        print(f"  [missing] {p}")

# %% [markdown]
# ## Cell vertices via healpix-geo (sphere, NESTED)
#
# `vertices(ipix, depth)` returns (N, 4) lon/lat arrays — the four
# corners of each HEALPix cell on the sphere. We wrap longitudes from
# [0, 360) to [-180, 180] for plotting. Cells that straddle the
# antimeridian are not an issue at Iberia latitudes — IBERIA is fully
# in the western hemisphere of the (-180, 180) frame.

# %%
lon_v, lat_v = hp_vertices(IBERIA_PIX_64, DEPTH, ellipsoid="WGS84")
# (N_64, 4) for both
lon_v = np.where(lon_v > 180.0, lon_v - 360.0, lon_v)
print(f"vertex array shapes: lon={lon_v.shape}  lat={lat_v.shape}")
print(
    f"lon vertex range: {lon_v.min():.2f} .. {lon_v.max():.2f}  "
    f"lat vertex range: {lat_v.min():.2f} .. {lat_v.max():.2f}"
)

# Build (N_64, 4, 2) polygon array for matplotlib.
poly_xy = np.stack([lon_v, lat_v], axis=-1)         # (N_64, 4, 2)

# %% [markdown]
# ## Helper: ranked bar chart on a given matplotlib axis

# %%
GOLD = "#d4a017"
DARK_GOLD = "#8a6a0c"
TEAL = "#2c7bb6"
ORANGE = "#d7191c"
DATA_FOOTER = (
    "Source: DestinE Climate DT SSP3-7.0 "
    "(licence-restricted; access via DestinE Data Lake)"
)


def _plot_rank(ax, records, title, color, order_species=None,
               highlight_species=None):
    """Horizontal bar chart of per-species posterior-mean p_extirp.

    If `order_species` is given, the species are reordered to match
    (any species not present in `records` are dropped); otherwise
    the input order of `records` is preserved.

    `highlight_species` (set or list) controls which bars are drawn
    in dark gold; defaults to the top-3 species *in the displayed
    order*.
    """
    by_name = {r["species"]: r for r in records}
    if order_species is not None:
        species = [sp for sp in order_species if sp in by_name]
    else:
        species = [r["species"] for r in records]

    means = np.array([by_name[sp]["post_mean_p_extirpation"] for sp in species])
    los = np.array([by_name[sp]["hdi95_low"] for sp in species])
    his = np.array([by_name[sp]["hdi95_high"] for sp in species])

    y = np.arange(len(species))
    err = np.vstack([means - los, his - means])

    if highlight_species is None:
        # Highlight the top-3 by posterior-mean within this panel
        top3 = set(sorted(species, key=lambda s: -by_name[s]["post_mean_p_extirpation"])[:3])
    else:
        top3 = set(highlight_species)
    bar_colors = [DARK_GOLD if sp in top3 else color for sp in species]

    ax.barh(y, means, color=bar_colors, alpha=0.85, edgecolor="white")
    ax.errorbar(means, y, xerr=err, fmt="none", ecolor="black",
                elinewidth=0.8, capsize=2, alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"B. {sp}" for sp in species], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Posterior-mean extirpation probability\n(95% HDI)")
    ax.set_title(title, fontsize=11)
    ax.set_xlim(0, max(0.05, means.max() * 1.25))
    ax.grid(axis="x", linewidth=0.3, alpha=0.5)


# %% [markdown]
# ## Helper: HEALPix-cell polygon map on a given cartopy axis

# %%
def _draw_healpix_map(ax, raster_per_cell, title):
    """Draw the 110-cell Iberia raster as colour-filled polygons."""
    proj = ccrs.PlateCarree()
    ax.set_extent([-10.5, 4.5, 35.0, 44.5], crs=proj)
    ax.add_feature(cfeature.LAND, facecolor="#f5f5f5", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="#e8f0fb", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, zorder=2)
    ax.add_feature(cfeature.BORDERS, linewidth=0.4, linestyle=":", zorder=2)

    valid = np.isfinite(raster_per_cell)
    if valid.sum() == 0:
        ax.set_title(title + " (no data)", fontsize=11)
        return None

    vmax = max(0.01, np.nanpercentile(raster_per_cell, 98))
    cmap = plt.get_cmap("YlOrRd")
    norm = plt.Normalize(vmin=0, vmax=vmax)

    polys = poly_xy[valid]
    vals = raster_per_cell[valid]
    pc = PolyCollection(
        polys,
        array=vals,
        cmap=cmap,
        norm=norm,
        edgecolors="black",
        linewidths=0.2,
        transform=proj,
        zorder=1,
    )
    ax.add_collection(pc)
    ax.set_title(title, fontsize=11)
    return pc


# %% [markdown]
# ## Figure 1 — species-rank chart (two panels)

# %%
fig, axes = plt.subplots(1, 2, figsize=(13, 9), sharey=False)
records_near = summary["horizons"]["2020_2029"]["species_ranked"]
records_mid = summary["horizons"]["2030_2039"]["species_ranked"]

# Sort *both* panels by mid-term posterior-mean, so the eye can track
# rank shifts across horizons. Highlight the top-3 most-vulnerable in
# each horizon (per its own ranking) in dark gold.
def _top3_species(records):
    return [r["species"] for r in
            sorted(records, key=lambda r: -r["post_mean_p_extirpation"])[:3]]


mid_order = [
    r["species"] for r in
    sorted(records_mid, key=lambda r: -r["post_mean_p_extirpation"])
]
top3_near = set(_top3_species(records_near))
top3_mid = set(_top3_species(records_mid))

_plot_rank(axes[0], records_near, HORIZON_TITLES["2020_2029"], TEAL,
           order_species=mid_order, highlight_species=top3_near)
_plot_rank(axes[1], records_mid, HORIZON_TITLES["2030_2039"], ORANGE,
           order_species=mid_order, highlight_species=top3_mid)
fig.suptitle(
    "Iberian Bombus extirpation risk under DestinE Climate DT SSP3-7.0\n"
    f"Top-3 most-vulnerable per horizon highlighted in gold "
    f"(N = {summary['method']['n_posterior_draws']} posterior draws, "
    f"HEALPix nside=64 NESTED)",
    fontsize=12,
)
fig.text(
    0.5, 0.005, DATA_FOOTER,
    ha="center", va="bottom", fontsize=8, color="dimgray", style="italic",
)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
out = FIG_DIR / "projection_species_rank.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"Saved {out}")

# %% [markdown]
# ## Figures 2 + 3 — per-cell community-mean risk maps

# %%

def _plot_map(horizon: str, raster: np.ndarray) -> Path:
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(7.5, 6))
    ax = plt.axes(projection=proj)
    pc = _draw_healpix_map(
        ax, raster,
        f"Iberian Bombus extirpation risk -- {HORIZON_TITLES[horizon]}",
    )
    if pc is not None:
        cbar = plt.colorbar(pc, ax=ax, orientation="vertical",
                            fraction=0.04, pad=0.03)
        cbar.set_label("Community-mean extirpation probability")
    fig.text(
        0.5, 0.01, DATA_FOOTER,
        ha="center", va="bottom", fontsize=8, color="dimgray", style="italic",
    )
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    out = FIG_DIR / f"projection_risk_map_{horizon}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.show()
    return out


for h in HORIZONS:
    if h not in per_cell:
        print(f"  [skip] {h}")
        continue
    out = _plot_map(h, per_cell[h])
    print(f"Saved {out}")


# %% [markdown]
# ## Figure 4 — combined summary panel for the Jupyter Book

# %%
# The combined panel is the headline figure for the FORRT Outcome.
# Select the horizon with the higher posterior-mean community risk
# (median across cells) — this is more representative of the broad
# shift than the per-cell max (which is sensitive to the rare species
# with very few historical cells).
def _median_risk(h):
    if h not in per_cell:
        return -np.inf
    arr = per_cell[h]
    finite = arr[np.isfinite(arr)]
    return float(np.median(finite)) if finite.size else -np.inf


impactful = "2030_2039"
if per_cell:
    impactful = max(per_cell.keys(), key=_median_risk)
print(f"Combined panel uses {impactful} for the map "
      f"(median community risk = {_median_risk(impactful):.3f}).")

fig = plt.figure(figsize=(15, 7.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1.0])

ax_rank = fig.add_subplot(gs[0, 0])
records = summary["horizons"][impactful]["species_ranked"]
_plot_rank(
    ax_rank, records,
    f"Species rank -- {HORIZON_TITLES[impactful]}",
    ORANGE if impactful == "2030_2039" else TEAL,
)

if impactful in per_cell:
    ax_map = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
    pc = _draw_healpix_map(
        ax_map, per_cell[impactful],
        f"Risk map -- {HORIZON_TITLES[impactful]}",
    )
    if pc is not None:
        plt.colorbar(
            pc, ax=ax_map, orientation="vertical",
            fraction=0.04, pad=0.03,
            label="Community-mean extirpation probability",
        )

fig.suptitle(
    "Iberian Bombus extirpation projection -- DestinE Climate DT SSP3-7.0",
    fontsize=13,
)
fig.text(
    0.5, 0.005, DATA_FOOTER,
    ha="center", va="bottom", fontsize=8, color="dimgray", style="italic",
)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
out = FIG_DIR / "projection_summary.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"Saved {out}")
