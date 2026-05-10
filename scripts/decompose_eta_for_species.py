"""Decompose the GLMM linear predictor η per species into the 10 fixed-effect
term contributions, plus the species random effect.

For each species' active cells under the projection horizon, computes:

    eta = beta_0 (Intercept)
        + beta_1 * sc_sampling
        + beta_2 * sc_TEI_bs
        + beta_3 * sc_TEI_delta
        + beta_4 * sc_TEI_bs * sc_TEI_delta
        + beta_5 * sc_PEI_bs
        + beta_6 * sc_PEI_delta
        + beta_7 * sc_PEI_bs * sc_PEI_delta
        + beta_8 * sc_TEI_bs * sc_PEI_bs
        + beta_9 * sc_TEI_delta * sc_PEI_delta
        + species_random_intercept

Reports the per-cell mean of each term (averaged over the species' cells)
so we can see which term is driving the per-species η. Useful for
explaining anomalies like terrestris's negative projected η.

Usage:
    python scripts/decompose_eta_for_species.py terrestris [--horizon 2030_2039]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
HPORT = ROOT / "healpix_port" / "outputs_iberia"
RESULTS = ROOT / "results"

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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("species", help="Lowercase Bombus epithet, e.g. 'terrestris'")
    ap.add_argument("--horizon", default="2030_2039",
                    choices=("2020_2029", "2030_2039"))
    args = ap.parse_args()
    sp = args.species
    horizon = args.horizon

    # Load future TEI/PEI for this horizon
    fut = xr.open_dataset(HPORT / f"climate_tei_pei_future_{horizon}_healpix.nc")
    species_arr = fut["species"].values.astype(str)
    if sp not in species_arr:
        raise SystemExit(f"Species '{sp}' not in {HPORT}/climate_tei_pei_future_{horizon}_healpix.nc")
    si = list(species_arr).index(sp)

    # Load standardisation constants from Tier-1 dataGLMM
    parquet = pd.read_parquet(HPORT / "dataGLMM_extinction.parquet")
    sampling_zarr = parquet["sampling"].values
    sampling_mu, sampling_sd = sampling_zarr.mean(), sampling_zarr.std(ddof=1)
    constants = {
        col: (parquet[col].mean(), parquet[col].std(ddof=1))
        for col in ("TEI_bs", "TEI_delta", "PEI_bs", "PEI_delta")
    }

    def z(x, col):
        mu, sd = constants[col]
        return (x - mu) / sd

    # VB posterior means (use these as point estimates for term decomposition;
    # for full posterior use the bambi NetCDF — kept simple here)
    post = pd.read_csv(HPORT / "posterior_vb_summary.csv", index_col=0)
    coef = {p: float(post.loc[p, "mean"]) for p in PARAM_NAMES}

    # Species random effect from bambi posterior
    idata = az.from_netcdf(RESULTS / "posterior_bambi_healpix.nc")
    re_ds = idata.posterior.get("1|species", None)
    if re_ds is None:
        # bambi sometimes calls it differently
        for k in idata.posterior.data_vars:
            if "species" in k.lower():
                re_ds = idata.posterior[k]
                break
    re_mean_per_species = re_ds.mean(dim=("chain", "draw")).values  # (n_species,)
    re_species_coords = (
        re_ds.coords.get("1|species_dim") or re_ds.coords.get("species_dim")
        or list(idata.posterior.coords.values())[-1]
    ).values.astype(str)
    if sp in re_species_coords:
        re_si = list(re_species_coords).index(sp)
        re_value = float(re_mean_per_species[re_si])
    else:
        re_value = 0.0
        print(f"  [warn] species '{sp}' not in random-effect coord; using RE = 0")

    # Future predictors at this species' active cells
    tei_bs = fut["tei_bs"].values[si, :]
    tei_delta = fut["tei_delta"].values[si, :]
    pei_bs = fut["pei_bs"].values[si, :]
    pei_delta = fut["pei_delta"].values[si, :]
    cells = fut["cell"].values

    valid = (np.isfinite(tei_bs) & np.isfinite(tei_delta)
             & np.isfinite(pei_bs) & np.isfinite(pei_delta))
    if valid.sum() == 0:
        raise SystemExit(f"No cells with all-finite predictors for {sp}")

    # Standardise (sampling held at recent-period mean = 0 in z-space)
    sc_TEI_bs = z(tei_bs[valid], "TEI_bs")
    sc_TEI_delta = z(tei_delta[valid], "TEI_delta")
    sc_PEI_bs = z(pei_bs[valid], "PEI_bs")
    sc_PEI_delta = z(pei_delta[valid], "PEI_delta")
    sc_sampling = np.zeros_like(sc_TEI_bs)

    # Per-cell term contributions
    contrib = {
        "Intercept": np.full_like(sc_TEI_bs, coef["Intercept"]),
        "sc_sampling": coef["sc_sampling"] * sc_sampling,
        "sc_TEI_bs": coef["sc_TEI_bs"] * sc_TEI_bs,
        "sc_TEI_delta": coef["sc_TEI_delta"] * sc_TEI_delta,
        "sc_TEI_bs:sc_TEI_delta":
            coef["sc_TEI_bs:sc_TEI_delta"] * sc_TEI_bs * sc_TEI_delta,
        "sc_PEI_bs": coef["sc_PEI_bs"] * sc_PEI_bs,
        "sc_PEI_delta": coef["sc_PEI_delta"] * sc_PEI_delta,
        "sc_PEI_bs:sc_PEI_delta":
            coef["sc_PEI_bs:sc_PEI_delta"] * sc_PEI_bs * sc_PEI_delta,
        "sc_TEI_bs:sc_PEI_bs":
            coef["sc_TEI_bs:sc_PEI_bs"] * sc_TEI_bs * sc_PEI_bs,
        "sc_TEI_delta:sc_PEI_delta":
            coef["sc_TEI_delta:sc_PEI_delta"] * sc_TEI_delta * sc_PEI_delta,
        "species_RE": np.full_like(sc_TEI_bs, re_value),
    }
    eta = sum(contrib.values())

    print(f"\n=== B. {sp} — {horizon} — η decomposition ===")
    print(f"Species cells with finite predictors: {valid.sum()}")
    print(f"Mean η across cells: {eta.mean():+.3f}\n")

    # Mean term contribution across cells (this is what drives the species' mean η)
    print(f"{'Term':<28}  {'β (VB mean)':>12}  {'mean cell value':>16}  {'mean contribution':>20}")
    print("-" * 84)
    rows = []
    for term in list(contrib.keys()):
        if term == "species_RE":
            beta = float("nan")
            mean_cell_value = float("nan")
        elif term == "Intercept":
            beta = coef["Intercept"]
            mean_cell_value = 1.0
        elif ":" in term:
            beta = coef[term]
            mean_cell_value = float("nan")  # not applicable to single value
        else:
            beta = coef[term]
            # Re-derive the mean of the standardised predictor for that term
            arr = {
                "sc_sampling": sc_sampling,
                "sc_TEI_bs": sc_TEI_bs,
                "sc_TEI_delta": sc_TEI_delta,
                "sc_PEI_bs": sc_PEI_bs,
                "sc_PEI_delta": sc_PEI_delta,
            }[term]
            mean_cell_value = float(arr.mean())
        mean_contrib = float(contrib[term].mean())
        rows.append((term, beta, mean_cell_value, mean_contrib))
        print(f"{term:<28}  {beta:>+12.4f}  "
              f"{(f'{mean_cell_value:+.3f}' if np.isfinite(mean_cell_value) else 'n/a'):>16}  "
              f"{mean_contrib:>+20.4f}")

    print("-" * 84)
    total = sum(r[3] for r in rows)
    print(f"{'TOTAL η (sum of contributions)':<28}  {'':>12}  {'':>16}  {total:>+20.4f}")

    # Highlight the dominant term(s)
    rows_sorted = sorted(rows, key=lambda r: -abs(r[3]))
    print(f"\nTop 3 contributing terms (by |contribution|):")
    for term, beta, _, contrib_val in rows_sorted[:3]:
        sign = "POSITIVE → ↑ η" if contrib_val > 0 else "NEGATIVE → ↓ η"
        print(f"  {term:<28} contributes {contrib_val:+.3f}  ({sign})")


if __name__ == "__main__":
    main()
