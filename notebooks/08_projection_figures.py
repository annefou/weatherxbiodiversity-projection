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

# %%
import _tier2_guard; _tier2_guard.ensure_destine_or_skip()

# %% [markdown]
# # 08 — Tier-2 projection figures
#
# Produces:
#
#   * `figures/projection_species_rank.png` — two-panel ranked bar chart
#     (near-term left, mid-term right) of per-species posterior
#     mean extirpation probability with 95 % HDI error bars.
#   * `figures/projection_risk_map_2020_2029.png` — Iberia map of
#     community-mean (over species) extirpation probability (near-term).
#   * `figures/projection_risk_map_2030_2039.png` — same, mid-term.
#   * `figures/projection_summary.png` — combined panel for the Jupyter
#     Book / nanopub Outcome draft (rank chart + the more impactful
#     map).

# %%
import json
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

# %%
plt.style.use("seaborn-v0_8-whitegrid")

ROOT = Path("..").resolve()
PORT = ROOT / "soroye_port"
OUT_DIR = PORT / "outputs_iberia"
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Load projection summary + per-cell rasters

# %%
with open(RESULTS_DIR / "projection_headline.json") as f:
    summary = json.load(f)

HORIZONS = ["2020_2029", "2030_2039"]
HORIZON_TITLES = {
    "2020_2029": "Near-term (2020–2029)",
    "2030_2039": "Mid-term (2030–2039)",
}

per_cell = {}
for h in HORIZONS:
    p = RESULTS_DIR / f"projection_per_cell_{h}.npy"
    if p.exists():
        per_cell[h] = np.load(p)
        print(f"  loaded {p.name}: shape {per_cell[h].shape}")
    else:
        print(f"  [missing] {p}")

# %% [markdown]
# ## Recover CEA cell centre lat/lon for cartopy plotting
#
# Same construction as Tier 1 — we keep this self-contained so the
# figure module does not import from `soroye_port/`.

# %%
RES_M = 100_000
X_MIN, X_MAX = -20_037_507, 20_062_493
Y_MIN, Y_MAX = -5_263_885, 6_336_115

cea_to_ll = Transformer.from_crs(
    "+proj=cea +lat_ts=0 +lon_0=0 +ellps=WGS84",
    "EPSG:4326",
    always_xy=True,
)
n_x = (X_MAX - X_MIN) // RES_M
n_y = (Y_MAX - Y_MIN) // RES_M
x_centers = X_MIN + (np.arange(n_x) + 0.5) * RES_M
y_centers = Y_MAX - (np.arange(n_y) + 0.5) * RES_M
xx, yy = np.meshgrid(x_centers, y_centers)
LON_CELL, LAT_CELL = cea_to_ll.transform(xx.ravel(), yy.ravel())
LON_CELL = LON_CELL.reshape((n_y, n_x))
LAT_CELL = LAT_CELL.reshape((n_y, n_x))

# Cell-edge meshes (one extra row + col) for pcolormesh.
xe = X_MIN + np.arange(n_x + 1) * RES_M
ye = Y_MAX - np.arange(n_y + 1) * RES_M
xxe, yye = np.meshgrid(xe, ye)
LON_E, LAT_E = cea_to_ll.transform(xxe.ravel(), yye.ravel())
LON_E = LON_E.reshape((n_y + 1, n_x + 1))
LAT_E = LAT_E.reshape((n_y + 1, n_x + 1))

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


def _plot_rank(ax, records, title, color):
    species = [r["species"] for r in records]
    means = np.array([r["post_mean_p_extirpation"] for r in records])
    los = np.array([r["hdi95_low"] for r in records])
    his = np.array([r["hdi95_high"] for r in records])

    y = np.arange(len(species))
    err = np.vstack([means - los, his - means])
    bar_colors = [DARK_GOLD if i < 3 else color for i in range(len(species))]
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
# ## Figure 1 — species-rank chart (two panels)

# %%
fig, axes = plt.subplots(1, 2, figsize=(13, 9), sharey=False)
records_mid = summary["horizons"]["2020_2029"]["species_ranked"]
records_end = summary["horizons"]["2030_2039"]["species_ranked"]

_plot_rank(axes[0], records_mid, HORIZON_TITLES["2020_2029"], TEAL)
_plot_rank(axes[1], records_end, HORIZON_TITLES["2030_2039"], ORANGE)
fig.suptitle(
    "Iberian Bombus extirpation risk under DestinE Climate DT SSP3-7.0\n"
    f"Top-3 most-vulnerable per horizon highlighted in gold "
    f"(N = {summary['method']['n_posterior_draws']} posterior draws)",
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
    # Iberia + adjacent extent.
    ax.set_extent([-10.5, 4.5, 35.0, 44.5], crs=proj)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linewidth=0.4, linestyle=":")
    ax.add_feature(cfeature.LAND, facecolor="#f5f5f5")
    ax.add_feature(cfeature.OCEAN, facecolor="#e8f0fb")

    # Mask cells outside the extent (NaN where no data).
    pcm = ax.pcolormesh(
        LON_E, LAT_E, raster,
        cmap="YlOrRd", vmin=0, vmax=max(0.01, np.nanpercentile(raster, 98)),
        transform=proj, shading="flat",
    )
    cbar = plt.colorbar(pcm, ax=ax, orientation="vertical",
                        fraction=0.04, pad=0.03)
    cbar.set_label("Community-mean extirpation probability")

    ax.set_title(
        f"Iberian Bombus extirpation risk — {HORIZON_TITLES[horizon]}",
        fontsize=11,
    )
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
        continue
    out = _plot_map(h, per_cell[h])
    print(f"Saved {out}")


# %% [markdown]
# ## Figure 4 — combined summary panel for the Jupyter Book

# %%
# Pick the "more impactful" horizon by max community-mean p.
if per_cell:
    impactful = max(
        per_cell.keys(),
        key=lambda h: float(np.nanmax(per_cell[h])) if np.isfinite(np.nanmax(per_cell[h])) else 0.0,
    )
else:
    impactful = "2030_2039"
print(f"Combined panel uses {impactful} for the map.")

fig = plt.figure(figsize=(15, 7.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1.0])

ax_rank = fig.add_subplot(gs[0, 0])
records = summary["horizons"][impactful]["species_ranked"]
_plot_rank(ax_rank, records,
           f"Species rank — {HORIZON_TITLES[impactful]}",
           ORANGE if impactful == "2030_2039" else TEAL)

if impactful in per_cell:
    ax_map = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
    ax_map.set_extent([-10.5, 4.5, 35.0, 44.5], crs=ccrs.PlateCarree())
    ax_map.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax_map.add_feature(cfeature.BORDERS, linewidth=0.4, linestyle=":")
    ax_map.add_feature(cfeature.LAND, facecolor="#f5f5f5")
    ax_map.add_feature(cfeature.OCEAN, facecolor="#e8f0fb")
    raster = per_cell[impactful]
    pcm = ax_map.pcolormesh(
        LON_E, LAT_E, raster,
        cmap="YlOrRd", vmin=0,
        vmax=max(0.01, np.nanpercentile(raster, 98)),
        transform=ccrs.PlateCarree(), shading="flat",
    )
    plt.colorbar(pcm, ax=ax_map, orientation="vertical",
                 fraction=0.04, pad=0.03,
                 label="Community-mean extirpation probability")
    ax_map.set_title(
        f"Risk map — {HORIZON_TITLES[impactful]}", fontsize=11,
    )

fig.suptitle(
    "Iberian Bombus extirpation projection — DestinE Climate DT SSP3-7.0",
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
