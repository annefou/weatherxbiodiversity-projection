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
# # 07 — Future-climate extirpation projection (Tier 2)
#
# Combines:
#
#   * `results/posterior_bambi.nc` — joint posterior of the GLMM fixed
#     effects from Tier 1 (2 chains × 2000 draws = 4 000 samples).
#   * `soroye_port/outputs_iberia/climate_tei_pei_future_<horizon>.npz`
#     — DestinE-derived future TEI / PEI on the Tier-1 CEA grid.
#   * `soroye_port/outputs_iberia/sampling_continent.npz` — recent-period
#     mean sampling effort per cell (held fixed for the projection).
#
# For each (species × active cell × posterior draw) we compute the
# linear predictor η = X · β and apply the logistic link to get an
# extirpation probability. We aggregate per species and write
# `results/projection_headline.json` with the ranking + 95 % HDI per
# horizon.
#
# **Critical step — predictor scaling**
#
# The Tier-1 GLMM was fit on z-scored predictors. The standardisation
# constants (mean + ddof=1 SD on `TEI_bs`, `TEI_delta`, `PEI_bs`,
# `PEI_delta`, `sampling`) are recoverable from the unscaled raw
# columns of `dataGLMM_extinction.parquet` — Tier 1 wrote both the raw
# and `sc_` columns. Apply the same mean / SD to the future predictors
# **before** plugging into the design matrix; using the future-period
# mean / SD instead would silently invalidate the projection.

# %%
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
from scipy.special import expit

# %%
ROOT = Path("..").resolve()
PORT = ROOT / "soroye_port"
OUT_DIR = PORT / "outputs_iberia"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

HORIZONS = ["2046_2055", "2076_2085"]
N_DRAWS = 1000
RNG = np.random.default_rng(20260509)

# %% [markdown]
# ## Recover Tier-1 standardisation constants from `dataGLMM_extinction.parquet`
#
# Each `sc_<col>` column is `(<col> - mean) / sd` with sd computed as
# the unbiased sample SD (`ddof=1`, matching R's `scale()` and Tier-1
# `soroye_port/05_regression.py`). Recompute the mean + sd from the
# raw columns; cross-check against the saved `sc_` columns to confirm.

# %%
parquet_path = OUT_DIR / "dataGLMM_extinction.parquet"
df = pd.read_parquet(parquet_path)
print(f"Loaded {parquet_path.name}: {df.shape[0]:,} rows")

SCALED_COLS = ["sampling", "TEI_bs", "TEI_delta", "PEI_bs", "PEI_delta"]
scaling = {}
for col in SCALED_COLS:
    mean = float(df[col].mean())
    sd = float(df[col].std(ddof=1))
    # Cross-check: re-derive sc_ column from raw and compare.
    rederived = (df[col] - mean) / sd
    diff = float(np.max(np.abs(rederived - df[f"sc_{col}"])))
    scaling[col] = {"mean": mean, "sd": sd, "max_check_diff": diff}
    print(f"  {col}: mean={mean:.6f}, sd={sd:.6f}, max |sc_check - sc_orig| = {diff:.2e}")

assert all(v["max_check_diff"] < 1e-6 for v in scaling.values()), (
    "Standardisation constants do not match the saved sc_ columns. "
    "Tier-1 z-score recovery is invalid."
)


def _z(arr: np.ndarray, col: str) -> np.ndarray:
    """Apply the Tier-1 z-score for column ``col`` to a future-predictor array."""
    return (arr - scaling[col]["mean"]) / scaling[col]["sd"]


# %% [markdown]
# ## Load the Tier-1 bambi posterior (joint covariance)
#
# Stack chains × draws into a flat (n_samples, n_params) matrix.
# Subsample ``N_DRAWS`` rows for the projection.

# %%
idata = az.from_netcdf(RESULTS_DIR / "posterior_bambi.nc")
posterior = idata.posterior

PARAM_NAMES = [
    "Intercept",
    "sc_sampling",
    "sc_TEI_bs",
    "sc_TEI_delta",
    "sc_TEI_bs:sc_TEI_delta",
    "sc_PEI_bs",
    "sc_PEI_delta",
    "sc_PEI_bs:sc_PEI_delta",
    "sc_TEI_bs:sc_PEI_bs",
    "sc_TEI_delta:sc_PEI_delta",
]

# Confirm every coefficient is in the posterior — if Tier 1 dropped
# `continent` (single-continent Iberian fit) the formula is a strict
# subset of the canonical Soroye 2020 formula.
missing = [p for p in PARAM_NAMES if p not in posterior.data_vars]
if missing:
    raise KeyError(f"Posterior is missing fixed-effect terms: {missing}")

flat_chain_draw = posterior.stack(sample=("chain", "draw"))
n_samples_total = flat_chain_draw.sizes["sample"]
print(f"Total posterior samples (chain × draw): {n_samples_total}")

# Beta matrix: each fixed-effect var is (sample,) after the stack.
beta_full = np.column_stack([
    flat_chain_draw[name].values for name in PARAM_NAMES
])
print(f"Beta matrix shape: {beta_full.shape}  (samples × params)")

# Per-species random intercept (1|species). arviz orders the dims as
# (chain, draw, species_factor_dim), so after stack(sample=...) the
# resulting array is (species_factor, sample) — we transpose to
# (sample, species_factor) so row-indexing by draw is direct.
species_re_da = flat_chain_draw["1|species"]
print(f"Random intercept dims (post-stack): {species_re_da.dims}, shape {species_re_da.shape}")
# Locate the species-factor dim (the one that is not 'sample').
factor_dim = [d for d in species_re_da.dims if d != "sample"][0]
species_re = species_re_da.transpose("sample", factor_dim).values  # (sample, n_species_re)

# arviz/bambi typically uses 'species__factor_dim' or '1|species_dim'.
# Pick whichever coord matches the recovered factor_dim.
if factor_dim in posterior.coords:
    species_factor_levels = [str(x) for x in posterior.coords[factor_dim].values]
else:
    species_factor_levels = None
print(f"Random-intercept species levels: {species_factor_levels}")

# Subsample N_DRAWS draws for the projection.
draw_idx = RNG.choice(n_samples_total, size=N_DRAWS, replace=False)
beta = beta_full[draw_idx]                         # (N_DRAWS, 10)
species_re_draws = species_re[draw_idx]            # (N_DRAWS, n_species_re)
print(f"Subsampled to {N_DRAWS} draws. species_re_draws shape {species_re_draws.shape}")

# %% [markdown]
# ## Load the cell mask + sampling-effort term
#
# Per the Limitations field of the Outcome: sampling effort is held at
# the recent-period (2000–2014) mean per cell. This is a known
# documented assumption — extirpation in the GLMM is partly explained
# by sampling effort, so projections under different observer effort
# would shift the absolute probabilities (but not, in the linear
# predictor, the species ranking that the predictor difference
# induces).

# %%
sc = np.load(OUT_DIR / "sampling_continent.npz", allow_pickle=True)
sampling = sc["samp_total"]                  # (n_cells,)
active_mask = ~np.isnan(sampling)
n_active = int(active_mask.sum())
print(f"Active cells (sampled at least once): {n_active}")

# Tier-1 species ordering (from the parquet) drives the column order of
# the random-intercept matrix in the bambi posterior.
species_in_parquet = sorted(df["species"].unique().tolist())
n_spp_data = len(species_in_parquet)
print(f"Species in dataGLMM: {n_spp_data}")

# Map each species to its column in `species_re_draws`.
if species_factor_levels is None:
    print("  WARNING: could not recover species factor levels — assuming "
          "ordering matches sorted unique species in dataGLMM.")
    species_factor_levels = species_in_parquet
spp_to_re_col = {sp: species_factor_levels.index(sp) for sp in species_in_parquet
                 if sp in species_factor_levels}

# %% [markdown]
# ## Build the design matrix per (species × cell × horizon)
#
# Column order MUST match `PARAM_NAMES`:
#
#   [Intercept, sc_sampling, sc_TEI_bs, sc_TEI_delta,
#    sc_TEI_bs:sc_TEI_delta, sc_PEI_bs, sc_PEI_delta,
#    sc_PEI_bs:sc_PEI_delta, sc_TEI_bs:sc_PEI_bs,
#    sc_TEI_delta:sc_PEI_delta]
#
# The linear predictor is η = X·β + species_re. We compute X cell-by-cell
# but we vectorise across (cells × draws) within each species.

# %%

def build_design_row(
    sc_sampling: np.ndarray,
    sc_TEI_bs: np.ndarray,
    sc_TEI_delta: np.ndarray,
    sc_PEI_bs: np.ndarray,
    sc_PEI_delta: np.ndarray,
) -> np.ndarray:
    """Stack the 10 design columns for an array of cells. Shape (n, 10)."""
    n = len(sc_sampling)
    return np.column_stack([
        np.ones(n),                              # Intercept
        sc_sampling,                             # sc_sampling
        sc_TEI_bs,                               # sc_TEI_bs
        sc_TEI_delta,                            # sc_TEI_delta
        sc_TEI_bs * sc_TEI_delta,                # sc_TEI_bs:sc_TEI_delta
        sc_PEI_bs,                               # sc_PEI_bs
        sc_PEI_delta,                            # sc_PEI_delta
        sc_PEI_bs * sc_PEI_delta,                # sc_PEI_bs:sc_PEI_delta
        sc_TEI_bs * sc_PEI_bs,                   # sc_TEI_bs:sc_PEI_bs
        sc_TEI_delta * sc_PEI_delta,             # sc_TEI_delta:sc_PEI_delta
    ])


# Pre-scale the sampling effort once (same for all horizons / species).
sc_sampling_active = _z(sampling[active_mask], "sampling")

# %% [markdown]
# ## Run the projection per horizon

# %%
projection_summary = {"horizons": {}, "method": {}}

for horizon in HORIZONS:
    src = OUT_DIR / f"climate_tei_pei_future_{horizon}.npz"
    if not src.exists():
        print(f"[skip] {src.name} missing — run 06_destine_clean.py first.")
        continue
    fut = np.load(src, allow_pickle=True)

    species_fut = list(fut["species"])
    if species_fut != species_in_parquet:
        # Order may differ (sorted-unique vs presence/absence original).
        # Reindex to the dataGLMM species order.
        order = [species_fut.index(sp) for sp in species_in_parquet
                 if sp in species_fut]
    else:
        order = list(range(len(species_fut)))

    TEI_bs_fut = fut["TEI_bs"][order]            # (n_spp, n_cells)
    PEI_bs_fut = fut["PEI_bs"][order]
    TEI_delta_fut = fut["TEI_delta"][order]
    PEI_delta_fut = fut["PEI_delta"][order]

    # Per-cell community-mean probability accumulator (across species).
    p_per_cell_sum = np.zeros(n_active, dtype=np.float64)
    n_species_per_cell = np.zeros(n_active, dtype=np.int64)

    species_records = []

    for i, sp in enumerate(species_in_parquet):
        # Future TEI/PEI for this species at active cells.
        tei_bs_act = TEI_bs_fut[i, active_mask]
        tei_dl_act = TEI_delta_fut[i, active_mask]
        pei_bs_act = PEI_bs_fut[i, active_mask]
        pei_dl_act = PEI_delta_fut[i, active_mask]

        # Drop cells where any predictor is NaN.
        valid = (
            np.isfinite(tei_bs_act)
            & np.isfinite(tei_dl_act)
            & np.isfinite(pei_bs_act)
            & np.isfinite(pei_dl_act)
            & np.isfinite(sc_sampling_active)
        )
        n_valid = int(valid.sum())
        if n_valid == 0:
            print(f"  {sp}: no valid cells; skipped")
            continue

        # Apply Tier-1 z-score constants.
        sc_TEI_bs_v = _z(tei_bs_act[valid], "TEI_bs")
        sc_TEI_delta_v = _z(tei_dl_act[valid], "TEI_delta")
        sc_PEI_bs_v = _z(pei_bs_act[valid], "PEI_bs")
        sc_PEI_delta_v = _z(pei_dl_act[valid], "PEI_delta")
        sc_sampling_v = sc_sampling_active[valid]

        X = build_design_row(
            sc_sampling_v, sc_TEI_bs_v, sc_TEI_delta_v,
            sc_PEI_bs_v, sc_PEI_delta_v,
        )                                  # (n_valid, 10)

        # η = X · β.T → (n_valid, N_DRAWS)
        eta = X @ beta.T

        # Add species random intercept (per draw).
        if sp in spp_to_re_col:
            re_col = spp_to_re_col[sp]
            re_per_draw = species_re_draws[:, re_col]    # (N_DRAWS,)
            eta = eta + re_per_draw[np.newaxis, :]

        p = expit(eta)                       # (n_valid, N_DRAWS)

        # Per-species summary across cells × draws (unweighted across
        # cells in the species' historical range, posterior mean over
        # draws per cell first, then mean over cells).
        p_post_mean_per_cell = p.mean(axis=1)        # (n_valid,)
        p_post_per_draw = p.mean(axis=0)             # (N_DRAWS,) species-level mean p across cells per draw
        post_mean_p = float(p_post_per_draw.mean())
        hdi = az.hdi(p_post_per_draw[np.newaxis, :], hdi_prob=0.95)
        # az.hdi on a 1-d array → returns shape (2,)
        hdi_low, hdi_high = float(hdi[0]), float(hdi[1])

        species_records.append({
            "species": sp,
            "post_mean_p_extirpation": post_mean_p,
            "hdi95_low": hdi_low,
            "hdi95_high": hdi_high,
            "n_cells": n_valid,
        })

        # Accumulate community-mean probability per cell (posterior mean across draws).
        valid_active_idx = np.flatnonzero(valid)
        p_per_cell_sum[valid_active_idx] += p_post_mean_per_cell
        n_species_per_cell[valid_active_idx] += 1

    # Sort species by posterior-mean extirpation probability (descending).
    species_records.sort(
        key=lambda r: r["post_mean_p_extirpation"], reverse=True,
    )

    # Per-cell community-mean probability raster (gitignored .npy).
    with np.errstate(invalid="ignore"):
        p_per_cell_mean = np.where(
            n_species_per_cell > 0,
            p_per_cell_sum / np.maximum(n_species_per_cell, 1),
            np.nan,
        )
    # Embed back into the full (n_x*n_y) flat grid for figure plotting.
    n_x, n_y = int(fut["n_x"]), int(fut["n_y"])
    full_grid = np.full(n_x * n_y, np.nan, dtype=np.float64)
    full_grid[active_mask] = p_per_cell_mean
    np.save(
        RESULTS_DIR / f"projection_per_cell_{horizon}.npy",
        full_grid.reshape((n_y, n_x)),
    )
    print(f"  Saved per-cell raster → results/projection_per_cell_{horizon}.npy "
          f"(gitignored)")

    projection_summary["horizons"][horizon] = {
        "n_draws": N_DRAWS,
        "species_ranked": species_records,
    }

    # Tabular preview.
    print(f"\n--- Top 5 most-vulnerable species ({horizon}) ---")
    for rec in species_records[:5]:
        print(
            f"  {rec['species']:<14}  post-mean p = {rec['post_mean_p_extirpation']:.3f}  "
            f"95% HDI [{rec['hdi95_low']:.3f}, {rec['hdi95_high']:.3f}]  "
            f"n_cells = {rec['n_cells']}"
        )

# %% [markdown]
# ## Method block + write JSON

# %%
projection_summary["method"] = {
    "n_posterior_draws": N_DRAWS,
    "scaling_source": (
        "Tier-1 dataGLMM_extinction.parquet z-score constants "
        "(mean + ddof=1 SD per raw column; cross-checked against the "
        "saved sc_ columns to ≤ 1e-6)"
    ),
    "scaling_constants": scaling,
    "sampling_effort_assumption": (
        "held at recent-period (2000–2014) mean per cell"
    ),
    "data_source": (
        "DestinE Climate Digital Twin SSP3-7.0 "
        "(licence-restricted; per-cell rasters not redistributed)"
    ),
    "design_columns": PARAM_NAMES,
    "posterior_total_samples": int(n_samples_total),
}

out_json = RESULTS_DIR / "projection_headline.json"
with open(out_json, "w") as f:
    json.dump(projection_summary, f, indent=2)
print(f"\nWrote {out_json}")
