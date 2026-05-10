"""Cross-substrate per-species η decomposition: nside=64 vs nside=128.

For a given species + horizon, decomposes the GLMM linear predictor η
into its 10 fixed-effect term contributions + species random intercept,
at BOTH substrates side-by-side. The aim is to identify which term
explains the η shift between the two substrates' projections — a
species-level diagnostic of substrate-sensitivity.

Reads from BOTH the reference repo (this one, nside=64) and the sibling
repo (weatherxbiodiversity-projection-nside128).

Usage:
    python scripts/compare_substrates.py norvegicus
    python scripts/compare_substrates.py mucidus --horizon 2030_2039
"""
from __future__ import annotations

import argparse
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import xarray as xr

REF64 = Path("/Users/annef/Documents/ScienceLive/weatherxbiodiversity-projection")
REF128 = Path("/Users/annef/Documents/ScienceLive/weatherxbiodiversity-projection-nside128")

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


def load_substrate(root: Path, horizon: str, label: str):
    """Load all the pieces needed to decompose η at one substrate."""
    out = {}

    # GLMM coefficient posterior (use VB summary CSV — fast + same as projection)
    post = pd.read_csv(
        root / "healpix_port" / "outputs_iberia" / "posterior_vb_summary.csv",
        index_col=0,
    )
    out["coef"] = {p: float(post.loc[p, "mean"]) for p in PARAM_NAMES}

    # Random effects from the bambi posterior (mean per species)
    idata = az.from_netcdf(root / "results" / "posterior_bambi_healpix.nc")
    re_da = idata.posterior["1|species"]
    factor_dim = [d for d in re_da.dims if d not in ("chain", "draw")][0]
    out["re_levels"] = [str(x) for x in idata.posterior.coords[factor_dim].values]
    out["re_means"] = re_da.mean(dim=("chain", "draw")).values  # (n_species,)

    # Future climate (per-cell predictors)
    if root == REF128:
        # New repo: single climate file at nside=128
        fut_path = root / "healpix_port" / "outputs_iberia" / f"climate_tei_pei_future_{horizon}_healpix.nc"
    else:
        # Reference repo: nside=64 climate file
        fut_path = root / "healpix_port" / "outputs_iberia" / f"climate_tei_pei_future_{horizon}_healpix.nc"
    fut = xr.open_dataset(fut_path)
    out["species"] = [str(s) for s in fut["species"].values]
    out["tei_bs"] = fut["tei_bs"].values
    out["tei_delta"] = fut["tei_delta"].values
    out["pei_bs"] = fut["pei_bs"].values
    out["pei_delta"] = fut["pei_delta"].values

    # Standardisation constants (recover from dataGLMM parquet)
    parquet = pd.read_parquet(root / "healpix_port" / "outputs_iberia" / "dataGLMM_extinction.parquet")
    out["scaling"] = {
        col: (parquet[col].mean(), parquet[col].std(ddof=1))
        for col in ("TEI_bs", "TEI_delta", "PEI_bs", "PEI_delta", "sampling")
    }

    # Per-species niche limits (from the historical climate file)
    hist_path = root / "healpix_port" / "outputs_iberia" / "climate_tei_pei_healpix.nc"
    hist = xr.open_dataset(hist_path)
    out["T_min_spp"] = hist["T_min_spp"].values
    out["T_max_spp"] = hist["T_max_spp"].values
    out["P_min_spp"] = hist["P_min_spp"].values
    out["P_max_spp"] = hist["P_max_spp"].values

    # Sampling effort per cell
    sc = xr.open_dataset(root / "healpix_port" / "outputs_iberia" / "sampling_continent_healpix.nc")
    out["sampling"] = sc["sampling_total"].values

    # Per-species observation mask
    pa = xr.open_dataset(root / "healpix_port" / "outputs_iberia" / "presence_absence_healpix.nc")
    out["prab_baseline"] = pa["prab_baseline"].values
    out["prab_recent"] = pa["prab_recent"].values
    out["pa_species"] = [str(s) for s in pa["species"].values]

    out["label"] = label
    return out


def decompose(sub: dict, sp: str):
    """Compute per-term mean η contribution for one species at one substrate."""
    if sp not in sub["species"]:
        return None
    si = sub["species"].index(sp)

    pa_si = sub["pa_species"].index(sp) if sp in sub["pa_species"] else None
    if pa_si is None:
        return None

    observed = (sub["prab_baseline"][pa_si] > 0) | (sub["prab_recent"][pa_si] > 0)
    active = ~np.isnan(sub["sampling"])
    cell_mask = observed & active
    if cell_mask.sum() == 0:
        return None

    tei_bs = sub["tei_bs"][si, cell_mask]
    tei_dl = sub["tei_delta"][si, cell_mask]
    pei_bs = sub["pei_bs"][si, cell_mask]
    pei_dl = sub["pei_delta"][si, cell_mask]
    samp = sub["sampling"][cell_mask]
    valid = (np.isfinite(tei_bs) & np.isfinite(tei_dl)
             & np.isfinite(pei_bs) & np.isfinite(pei_dl)
             & np.isfinite(samp))
    if valid.sum() == 0:
        return None

    def z(arr, col):
        m, s = sub["scaling"][col]
        return (arr - m) / s

    sc_TEI_bs = z(tei_bs[valid], "TEI_bs")
    sc_TEI_delta = z(tei_dl[valid], "TEI_delta")
    sc_PEI_bs = z(pei_bs[valid], "PEI_bs")
    sc_PEI_delta = z(pei_dl[valid], "PEI_delta")
    sc_sampling = z(samp[valid], "sampling")

    coef = sub["coef"]
    contrib = {
        "Intercept": coef["Intercept"] * np.ones_like(sc_TEI_bs),
        "sc_sampling": coef["sc_sampling"] * sc_sampling,
        "sc_TEI_bs": coef["sc_TEI_bs"] * sc_TEI_bs,
        "sc_TEI_delta": coef["sc_TEI_delta"] * sc_TEI_delta,
        "sc_TEI_bs:sc_TEI_delta": coef["sc_TEI_bs:sc_TEI_delta"] * sc_TEI_bs * sc_TEI_delta,
        "sc_PEI_bs": coef["sc_PEI_bs"] * sc_PEI_bs,
        "sc_PEI_delta": coef["sc_PEI_delta"] * sc_PEI_delta,
        "sc_PEI_bs:sc_PEI_delta": coef["sc_PEI_bs:sc_PEI_delta"] * sc_PEI_bs * sc_PEI_delta,
        "sc_TEI_bs:sc_PEI_bs": coef["sc_TEI_bs:sc_PEI_bs"] * sc_TEI_bs * sc_PEI_bs,
        "sc_TEI_delta:sc_PEI_delta": coef["sc_TEI_delta:sc_PEI_delta"] * sc_TEI_delta * sc_PEI_delta,
    }
    re_value = 0.0
    if sp in sub["re_levels"]:
        re_idx = sub["re_levels"].index(sp)
        re_value = float(sub["re_means"][re_idx])
    contrib["species_RE"] = re_value * np.ones_like(sc_TEI_bs)

    return {
        "n_cells": int(valid.sum()),
        "raw": {
            "tei_bs_mean":    float(tei_bs[valid].mean()),
            "tei_delta_mean": float(tei_dl[valid].mean()),
            "pei_bs_mean":    float(pei_bs[valid].mean()),
            "pei_delta_mean": float(pei_dl[valid].mean()),
        },
        "z": {
            "sc_TEI_bs_mean":    float(sc_TEI_bs.mean()),
            "sc_TEI_delta_mean": float(sc_TEI_delta.mean()),
            "sc_PEI_bs_mean":    float(sc_PEI_bs.mean()),
            "sc_PEI_delta_mean": float(sc_PEI_delta.mean()),
        },
        "term_means": {k: float(v.mean()) for k, v in contrib.items()},
        "eta_mean": float(sum(v.mean() for v in contrib.values())),
        "T_min_spp": float(sub["T_min_spp"][si]),
        "T_max_spp": float(sub["T_max_spp"][si]),
        "P_min_spp": float(sub["P_min_spp"][si]),
        "P_max_spp": float(sub["P_max_spp"][si]),
    }


def compare(sp: str, horizon: str = "2030_2039"):
    sub64 = load_substrate(REF64, horizon, "nside=64")
    sub128 = load_substrate(REF128, horizon, "nside=128")
    d64 = decompose(sub64, sp)
    d128 = decompose(sub128, sp)
    if d64 is None or d128 is None:
        print(f"  [skip] {sp}: not present at one of the substrates")
        return

    print(f"\n{'='*84}")
    print(f"  B. {sp}  ({horizon})")
    print(f"{'='*84}")
    print(f"{'':<28} {'nside=64':>15} {'nside=128':>15} {'Δ':>15}")
    print("-" * 84)

    print(f"{'n_cells (in fit)':<28} {d64['n_cells']:>15} {d128['n_cells']:>15} "
          f"{d128['n_cells'] - d64['n_cells']:>+15}")

    print(f"{'T_min_spp (degC)':<28} {d64['T_min_spp']:>+15.2f} {d128['T_min_spp']:>+15.2f} "
          f"{d128['T_min_spp'] - d64['T_min_spp']:>+15.2f}")
    print(f"{'T_max_spp (degC)':<28} {d64['T_max_spp']:>+15.2f} {d128['T_max_spp']:>+15.2f} "
          f"{d128['T_max_spp'] - d64['T_max_spp']:>+15.2f}")
    print(f"{'T_range (degC)':<28} {d64['T_max_spp'] - d64['T_min_spp']:>+15.2f} "
          f"{d128['T_max_spp'] - d128['T_min_spp']:>+15.2f} "
          f"{(d128['T_max_spp'] - d128['T_min_spp']) - (d64['T_max_spp'] - d64['T_min_spp']):>+15.2f}")

    print()
    print(f"{'mean per-cell raw value':<28}")
    for col in ("tei_bs_mean", "tei_delta_mean", "pei_bs_mean", "pei_delta_mean"):
        print(f"  {col:<26} {d64['raw'][col]:>+15.4f} {d128['raw'][col]:>+15.4f} "
              f"{d128['raw'][col] - d64['raw'][col]:>+15.4f}")

    print()
    print(f"{'mean per-cell z-score':<28}")
    for col in ("sc_TEI_bs_mean", "sc_TEI_delta_mean", "sc_PEI_bs_mean", "sc_PEI_delta_mean"):
        print(f"  {col:<26} {d64['z'][col]:>+15.3f} {d128['z'][col]:>+15.3f} "
              f"{d128['z'][col] - d64['z'][col]:>+15.3f}")

    print()
    print(f"{'mean η term contribution':<28}")
    for term, v64 in d64["term_means"].items():
        v128 = d128["term_means"][term]
        d = v128 - v64
        flag = "  <--" if abs(d) > 1.0 else ""
        print(f"  {term:<26} {v64:>+15.3f} {v128:>+15.3f} {d:>+15.3f}{flag}")

    print("-" * 84)
    print(f"{'TOTAL η':<28} {d64['eta_mean']:>+15.3f} {d128['eta_mean']:>+15.3f} "
          f"{d128['eta_mean'] - d64['eta_mean']:>+15.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("species", nargs="+", help="Lowercase Bombus epithets")
    ap.add_argument("--horizon", default="2030_2039", choices=("2020_2029", "2030_2039"))
    args = ap.parse_args()
    for sp in args.species:
        compare(sp, args.horizon)
