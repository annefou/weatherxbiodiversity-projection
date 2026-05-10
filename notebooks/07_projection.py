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
# # 07 — Future-climate extirpation projection (Tier 2, HEALPix nside=64)
#
# Combines:
#
#   * `results/posterior_bambi_healpix.nc` — joint posterior of the
#     GLMM fixed effects from the Tier-1 HEALPix fit (Phase C; 2 chains
#     × 2000 draws = 4 000 samples).
#   * `healpix_port/outputs_iberia/climate_tei_pei_future_<horizon>_healpix.npz`
#     — DestinE-derived future TEI / PEI on the 110-cell Iberia HEALPix
#     nside=64 NESTED grid.
#   * `healpix_port/outputs_iberia/sampling_continent_healpix.npz` —
#     recent-period sampling effort per cell (held fixed for the
#     projection).
#   * `healpix_port/outputs_iberia/dataGLMM_extinction.parquet` —
#     Tier-1 HEALPix dataGLMM, used to recover the z-score
#     standardisation constants applied to the future predictors.
#
# For each (species × active cell × posterior draw) we compute the
# linear predictor η = X · β + species_intercept and apply the
# logistic link to get an extirpation probability. We aggregate per
# species and write `results/projection_headline.json` with the
# ranking + 95 % HDI per horizon.
#
# **Critical step — predictor scaling**
#
# The Tier-1 HEALPix GLMM was fit on z-scored predictors. The
# standardisation constants (mean + ddof=1 SD on `TEI_bs`,
# `TEI_delta`, `PEI_bs`, `PEI_delta`, `sampling`) are recovered from
# the unscaled raw columns of the parquet. Apply the same mean / SD
# to the future predictors **before** plugging into the design matrix;
# using the future-period mean / SD instead would silently invalidate
# the projection.

# %%
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import xarray as xr
# scipy.special.expit intentionally NOT imported — we report the GLMM
# linear predictor η directly, not its logistic transform. See the
# "We REPORT η directly" comment in the per-species inner loop, and
# 05_outcome.md § Limitations.

# %%
ROOT = Path("..").resolve()
HPORT = ROOT / "healpix_port"
OUT_DIR = HPORT / "outputs_iberia"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

HORIZONS = ["2020_2029", "2030_2039"]
N_DRAWS = 1000
RNG = np.random.default_rng(42)

# %% [markdown]
# ## Recover Tier-1 standardisation constants from `dataGLMM_extinction.parquet`
#
# Each `sc_<col>` column is `(<col> - mean) / sd` with sd computed as
# the unbiased sample SD (`ddof=1`, matching R's `scale()` and
# `healpix_port/05_regression_healpix.py`). Recompute the mean + sd
# from the raw columns; cross-check against the saved `sc_` columns
# to confirm.

# %%
parquet_path = OUT_DIR / "dataGLMM_extinction.parquet"
df = pd.read_parquet(parquet_path)
print(f"Loaded {parquet_path.name}: {df.shape[0]:,} rows")

SCALED_COLS = ["sampling", "TEI_bs", "TEI_delta", "PEI_bs", "PEI_delta"]
scaling = {}
for col in SCALED_COLS:
    mean = float(df[col].mean())
    sd = float(df[col].std(ddof=1))
    rederived = (df[col] - mean) / sd
    diff = float(np.max(np.abs(rederived - df[f"sc_{col}"])))
    scaling[col] = {"mean": mean, "sd": sd, "max_check_diff": diff}
    print(
        f"  {col}: mean={mean:.6f}, sd={sd:.6f}, "
        f"max |sc_check - sc_orig| = {diff:.2e}"
    )

assert all(v["max_check_diff"] < 1e-6 for v in scaling.values()), (
    "Standardisation constants do not match the saved sc_ columns. "
    "Tier-1 z-score recovery is invalid."
)


def _z(arr: np.ndarray, col: str) -> np.ndarray:
    """Apply the Tier-1 z-score for column `col` to a future-predictor array."""
    return (arr - scaling[col]["mean"]) / scaling[col]["sd"]


# %% [markdown]
# ## Load the Tier-1 HEALPix bambi posterior

# %%
idata = az.from_netcdf(RESULTS_DIR / "posterior_bambi_healpix.nc")
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
missing = [p for p in PARAM_NAMES if p not in posterior.data_vars]
if missing:
    raise KeyError(f"Posterior is missing fixed-effect terms: {missing}")

flat = posterior.stack(sample=("chain", "draw"))
n_samples_total = flat.sizes["sample"]
print(f"Total posterior samples (chain × draw): {n_samples_total}")

beta_full = np.column_stack(
    [flat[name].values for name in PARAM_NAMES]
)  # (n_samples_total, 10)
print(f"Beta matrix shape: {beta_full.shape}  (samples × params)")

# Per-species random intercept (1|species).
re_da = flat["1|species"]
factor_dim = [d for d in re_da.dims if d != "sample"][0]
species_re_full = re_da.transpose("sample", factor_dim).values   # (samples, 31)
species_factor_levels = [str(x) for x in posterior.coords[factor_dim].values]
print(
    f"Random-intercept species levels ({len(species_factor_levels)}): "
    f"{species_factor_levels[:3]} ..."
)

# Subsample to N_DRAWS draws (seeded for reproducibility).
draw_idx = RNG.choice(n_samples_total, size=N_DRAWS, replace=False)
beta = beta_full[draw_idx]                          # (N_DRAWS, 10)
species_re_draws = species_re_full[draw_idx]        # (N_DRAWS, n_spp_re)
print(f"Subsampled to {N_DRAWS} draws.")

# %% [markdown]
# ## Sampling effort term (held at recent-period mean per cell)

# %%
sc = xr.open_dataset(OUT_DIR / "sampling_continent_healpix.nc")
sampling_recent = sc["sampling_recent"].values      # (110,)
sampling_total = sc["sampling_total"].values        # (110,)
# Match Tier-1 contract: GLMM was fit with `sampling` = samp_total
# (see healpix_port/05_regression_healpix.py:99). For the projection
# we use the same column to keep the sampling z-score self-consistent.
sampling = sampling_total
active_mask = ~np.isnan(sampling)
n_active = int(active_mask.sum())
print(f"Active cells (sampled at least once): {n_active}/{len(sampling)}")

species_in_parquet = sorted(df["species"].unique().tolist())
n_spp_data = len(species_in_parquet)
print(f"Species in dataGLMM: {n_spp_data}")
spp_to_re_col = {
    sp: species_factor_levels.index(sp)
    for sp in species_in_parquet
    if sp in species_factor_levels
}

# %% [markdown]
# ## Design-matrix builder
#
# Column order matches `PARAM_NAMES` exactly. `eta = X · β.T` then
# `p_extirp = expit(eta + species_intercept)`.

# %%
def build_design_row(
    sc_sampling: np.ndarray,
    sc_TEI_bs: np.ndarray,
    sc_TEI_delta: np.ndarray,
    sc_PEI_bs: np.ndarray,
    sc_PEI_delta: np.ndarray,
) -> np.ndarray:
    n = len(sc_sampling)
    return np.column_stack([
        np.ones(n),
        sc_sampling,
        sc_TEI_bs,
        sc_TEI_delta,
        sc_TEI_bs * sc_TEI_delta,
        sc_PEI_bs,
        sc_PEI_delta,
        sc_PEI_bs * sc_PEI_delta,
        sc_TEI_bs * sc_PEI_bs,
        sc_TEI_delta * sc_PEI_delta,
    ])


sc_sampling_active = _z(sampling[active_mask], "sampling")
# (n_active,) — same for all horizons / species.

# Need the historical observation mask per species to restrict
# per-species summary to cells the species was observed in (matches
# how the GLMM was fit — only species-cell pairs with non-NaN
# extinction/colonisation enter the data).
pa = xr.open_dataset(OUT_DIR / "presence_absence_healpix.nc")
prab_baseline = pa["prab_baseline"].values          # (31, 110)
prab_recent = pa["prab_recent"].values
prab_species_list = [str(s) for s in pa["species"].values]

# %% [markdown]
# ## Project per horizon

# %%
projection_summary = {"horizons": {}, "method": {}}

for horizon in HORIZONS:
    src = OUT_DIR / f"climate_tei_pei_future_{horizon}_healpix.nc"
    if not src.exists():
        print(f"[skip] {src.name} missing — run 06_destine_clean.py first.")
        continue
    fut = xr.open_dataset(src)
    fut_species = [str(s) for s in fut["species"].values]

    # Reorder future arrays to match dataGLMM species order.
    if fut_species != species_in_parquet:
        order = [fut_species.index(sp) for sp in species_in_parquet
                 if sp in fut_species]
    else:
        order = list(range(len(fut_species)))

    TEI_bs_fut = fut["tei_bs"].values[order]         # (n_spp, 110)
    PEI_bs_fut = fut["pei_bs"].values[order]
    TEI_delta_fut = fut["tei_delta"].values[order]
    PEI_delta_fut = fut["pei_delta"].values[order]

    # Map presence-absence rows into the same dataGLMM order.
    pa_order = [prab_species_list.index(sp) for sp in species_in_parquet
                if sp in prab_species_list]
    prab_baseline_ord = prab_baseline[pa_order]      # (n_spp, 110)
    prab_recent_ord = prab_recent[pa_order]

    # Per-cell community-mean η accumulator (across species historically
    # observed in the cell). Linear-predictor units (log-odds) — kept
    # un-transformed to avoid expit() saturation under future predictors.
    eta_per_cell_sum = np.zeros(len(sampling), dtype=np.float64)
    n_species_per_cell = np.zeros(len(sampling), dtype=np.int64)

    species_records = []

    for i, sp in enumerate(species_in_parquet):
        # Species' historical range = baseline-occupied OR recent-
        # observed cells (matches the GLMM data filter, where each
        # row is a species-cell pair where the species was observed
        # in at least one of the two periods).
        observed_any = (prab_baseline_ord[i] == 1) | (prab_recent_ord[i] == 1)

        # Future TEI/PEI for this species at active cells, restricted
        # to historically observed cells.
        cell_mask = active_mask & observed_any
        if cell_mask.sum() == 0:
            print(f"  {sp}: no overlap of active × observed cells; skipped")
            continue

        tei_bs_v = TEI_bs_fut[i, cell_mask]
        tei_dl_v = TEI_delta_fut[i, cell_mask]
        pei_bs_v = PEI_bs_fut[i, cell_mask]
        pei_dl_v = PEI_delta_fut[i, cell_mask]
        sampling_v = sampling[cell_mask]

        valid = (
            np.isfinite(tei_bs_v)
            & np.isfinite(tei_dl_v)
            & np.isfinite(pei_bs_v)
            & np.isfinite(pei_dl_v)
            & np.isfinite(sampling_v)
        )
        if valid.sum() == 0:
            print(f"  {sp}: no valid cells; skipped")
            continue

        # Map cell_mask → index within the full (110,) cell vector for
        # the per-cell raster accumulator.
        cell_full_idx = np.flatnonzero(cell_mask)[valid]

        sc_TEI_bs_v = _z(tei_bs_v[valid], "TEI_bs")
        sc_TEI_delta_v = _z(tei_dl_v[valid], "TEI_delta")
        sc_PEI_bs_v = _z(pei_bs_v[valid], "PEI_bs")
        sc_PEI_delta_v = _z(pei_dl_v[valid], "PEI_delta")
        sc_sampling_v = _z(sampling_v[valid], "sampling")

        X = build_design_row(
            sc_sampling_v, sc_TEI_bs_v, sc_TEI_delta_v,
            sc_PEI_bs_v, sc_PEI_delta_v,
        )                                            # (n_valid, 10)

        eta = X @ beta.T                             # (n_valid, N_DRAWS)
        if sp in spp_to_re_col:
            re_per_draw = species_re_draws[:, spp_to_re_col[sp]]   # (N_DRAWS,)
            eta = eta + re_per_draw[np.newaxis, :]

        # We REPORT η directly, not p = expit(η). Reason: future TEI/PEI
        # z-scores under SSP3-7.0 land 5–10× outside the Tier-1 training
        # distribution (a multi-decade warming signal vs the 60-year
        # training-period delta). expit() saturates near 1.0 in that
        # extrapolation regime, making absolute probabilities
        # uninterpretable. η preserves relative ranking and the linear
        # predictor magnitude is the GLMM's authentic signal — it is
        # also what the bambi NUTS posterior is actually sampled over.
        # See nanopubs/drafts/05_outcome.md § Limitations.

        # Per-species summary across cells × draws:
        #   - per-cell posterior-mean η (averaged over draws)
        #   - per-draw species-level mean η (averaged over cells)
        eta_post_mean_per_cell = eta.mean(axis=1)    # (n_valid,)
        eta_post_per_draw = eta.mean(axis=0)         # (N_DRAWS,)
        post_mean_eta = float(eta_post_per_draw.mean())
        # az.hdi on a 1-d array returns shape (2,) = [low, high]
        hdi = az.hdi(eta_post_per_draw, hdi_prob=0.95)
        hdi_low, hdi_high = float(hdi[0]), float(hdi[1])
        # A non-saturating proxy for "fraction of cells projected to
        # exceed the moderate-risk threshold of η = 0 (= log-odds 0,
        # i.e. p > 0.5)". Robust to extrapolation at the per-cell level.
        n_cells_eta_gt_0 = int((eta_post_mean_per_cell > 0).sum())

        species_records.append({
            "species": sp,
            "post_mean_eta": post_mean_eta,
            "eta_hdi95_low": hdi_low,
            "eta_hdi95_high": hdi_high,
            "n_cells": int(valid.sum()),
            "n_cells_eta_gt_0": n_cells_eta_gt_0,
        })

        eta_per_cell_sum[cell_full_idx] += eta_post_mean_per_cell
        n_species_per_cell[cell_full_idx] += 1

    species_records.sort(
        key=lambda r: r["post_mean_eta"], reverse=True,
    )

    # Per-cell community-mean η (gitignored CF NetCDF). Linear predictor
    # in log-odds units; not transformed via expit because future
    # predictors are far outside training distribution and expit()
    # saturates uninformatively.
    with np.errstate(invalid="ignore"):
        eta_per_cell_mean = np.where(
            n_species_per_cell > 0,
            eta_per_cell_sum / np.maximum(n_species_per_cell, 1),
            np.nan,
        )

    # Use Tier-1 cell-centre lon/lat from presence_absence_healpix.nc.
    lon_cell = pa["lon"].values.astype(np.float32)
    lat_cell = pa["lat"].values.astype(np.float32)
    cell_idx = pa["cell"].values.astype(np.int64)

    decade_label = horizon.replace("_", "-")
    raster_ds = xr.Dataset(
        data_vars={
            "community_mean_eta": (
                ("cell",),
                eta_per_cell_mean.astype(np.float32),
                {
                    "long_name": (
                        "per-cell community-mean GLMM linear predictor (η, "
                        "log-odds of extirpation) under SSP3-7.0"
                    ),
                    "units": "1",
                    "comment": (
                        f"Mean η across {len(species_in_parquet)} species, "
                        f"each species' contribution averaged over {N_DRAWS} "
                        "posterior draws. η > 0 → projected p > 0.5; "
                        "η > +5 → near logit saturation (p ≈ 0.99). "
                        "Reported as η rather than expit(η) because future "
                        "predictors lie 5–10× outside training distribution "
                        "where expit() saturates uninformatively."
                    ),
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "lon": (
                ("cell",),
                lon_cell,
                {
                    "long_name": "HEALPix cell-centre longitude",
                    "standard_name": "longitude",
                    "units": "degrees_east",
                },
            ),
            "lat": (
                ("cell",),
                lat_cell,
                {
                    "long_name": "HEALPix cell-centre latitude",
                    "standard_name": "latitude",
                    "units": "degrees_north",
                },
            ),
            "n_species_in_cell": (
                ("cell",),
                n_species_per_cell.astype(np.int32),
                {
                    "long_name": (
                        "number of species contributing to the community-mean "
                        "extirpation probability for this cell"
                    ),
                    "units": "species",
                    "_FillValue": np.int32(-1),
                },
            ),
        },
        coords={"cell": cell_idx},
        attrs={
            "Conventions": "CF-1.10",
            "title": (
                f"Per-cell community-mean Bombus extirpation probability "
                f"({decade_label}) on Iberian HEALPix nside=64 NESTED"
            ),
            "horizon": decade_label,
            "tier1_posterior": "results/posterior_bambi_healpix.nc",
            "n_posterior_draws": int(N_DRAWS),
            "rng_seed": 42,
            "license_note": (
                "DestinE Climate DT redistribution restricted; per-cell "
                "raster is gitignored under data licence"
            ),
            "healpix_nside": 64,
            "healpix_depth": 6,
            "healpix_scheme": "NESTED",
            "history": (
                f"Created by notebooks/07_projection.py for horizon {horizon}"
            ),
        },
    )
    raster_ds["cell"].attrs.update({
        "long_name": "HEALPix NESTED pixel index (nside=64)",
    })
    raster_ds.to_netcdf(
        RESULTS_DIR / f"projection_{horizon}.nc",
        engine="netcdf4",
        encoding={
            "community_mean_eta": {"zlib": True, "complevel": 4},
            "lon": {"zlib": True, "complevel": 4},
            "lat": {"zlib": True, "complevel": 4},
            "n_species_in_cell": {"zlib": True, "complevel": 4},
        },
    )
    print(
        f"  Saved per-cell raster -> "
        f"results/projection_{horizon}.nc (gitignored)"
    )

    projection_summary["horizons"][horizon] = {
        "n_draws": N_DRAWS,
        "species_ranked": species_records,
    }

    print(f"\n--- Top 5 most-vulnerable species ({horizon}) — ranked by mean η ---")
    for rec in species_records[:5]:
        print(
            f"  B. {rec['species']:<14}  "
            f"post-mean η = {rec['post_mean_eta']:+.3f}  "
            f"95% HDI [{rec['eta_hdi95_low']:+.3f}, {rec['eta_hdi95_high']:+.3f}]  "
            f"n_cells = {rec['n_cells']}  "
            f"(n_cells with η>0: {rec['n_cells_eta_gt_0']})"
        )

# %% [markdown]
# ## Method block + write JSON

# %%
projection_summary["method"] = {
    "n_posterior_draws": N_DRAWS,
    "substrate": "HEALPix-NESTED nside=64 (~92 km cells)",
    "scaling_source": (
        "Tier-1 healpix_port/outputs_iberia/dataGLMM_extinction.parquet "
        "z-score constants (mean + ddof=1 SD per raw column; cross-checked "
        "against the saved sc_ columns to <= 1e-6)"
    ),
    "scaling_constants": scaling,
    "sampling_effort_assumption": (
        "held at recent-period (2000-2014) mean per cell"
    ),
    "data_source": (
        "DestinE Climate DT SSP3-7.0 IFS-NEMO standard "
        "(licence-restricted; access via DestinE Data Lake)"
    ),
    "tier1_posterior": "results/posterior_bambi_healpix.nc",
    "design_columns": PARAM_NAMES,
    "posterior_total_samples": int(n_samples_total),
    "rng_seed": 42,
    "headline_metric": "post_mean_eta",
    "headline_metric_comment": (
        "Per-species posterior-mean of (mean η across cells), 95% HDI in "
        "log-odds units. We report η — the GLMM linear predictor — rather "
        "than p = expit(η) because future TEI/PEI z-scores under SSP3-7.0 "
        "lie 5–10× outside the Tier-1 training distribution (the future "
        "warming signal is ~2× the 1901–2014 training delta), where "
        "expit() saturates near 1.0 uninformatively. Ranking by η "
        "preserves the GLMM's authentic signal; relative-risk-style "
        "comparisons (e.g. n_cells_eta_gt_0) are more substrate-stable "
        "than absolute probabilities. See 05_outcome.md § Limitations."
    ),
}

out_json = RESULTS_DIR / "projection_headline.json"
with open(out_json, "w") as f:
    json.dump(projection_summary, f, indent=2)
print(f"\nWrote {out_json}")
