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
# # 03h — Analysis (HEALPix-NESTED nside=64 substrate)
#
# **Phase C — substrate-robustness branch.** Runs the bambi MCMC GLMM
# and the statsmodels variational-Bayes GLMM on the HEALPix nside=64
# NESTED substrate, then writes a comparison JSON
# (`results/headline_statistic_healpix.json`) that lays the HEALPix
# fit side-by-side with weatherxbio v0.2.1's published CEA value
# (`+0.479`, the external reference).
#
# **Strict tolerance test** — the headline coefficient `sc_TEI_delta`
# from the HEALPix VB fit is checked against three conditions:
#
# 1. **Sign**: `mean > 0` (same direction as CEA).
# 2. **Significance**: `p_2sided < 0.05` (VB summary z-test).
# 3. **Magnitude**: `0.7 * 0.479 <= mean <= 1.3 * 0.479` (within +-30%
#    of the CEA reference of +0.479).
#
# All three must pass for a "Substrate-robust" verdict. If any fails,
# the verdict is "Substrate-divergent" and the HEALPix fit is reported
# transparently with the failed condition flagged. Honest negative
# results are publishable; overclaiming is not (DOMAIN.md).

# %%
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# %%
ROOT = Path("..").resolve()
PORT = ROOT / "healpix_port"
OUT_DIR = PORT / "outputs_iberia"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

env = {**os.environ}


def run(script: str, *, check: bool = True) -> int:
    print(f"\n=== {script} ===", flush=True)
    result = subprocess.run(
        [sys.executable, script],
        cwd=PORT,
        env=env,
        check=check,
    )
    return result.returncode


# %% [markdown]
# ## Bayesian GLMM (bambi / PyMC) on HEALPix nside=64
#
# Same NUTS sampler config as the CEA fit (2 chains x 2000 draws +
# 1000 tune). If bambi / pymc / pytensor is unavailable, the script
# logs the failure and exits 0; the statsmodels VB fit downstream still
# produces a usable posterior summary.

# %%
bambi_ok = run("05_regression_healpix.py", check=False) == 0

parquet_path = OUT_DIR / "dataGLMM_extinction.parquet"
if not parquet_path.exists():
    raise SystemExit(
        f"Expected {parquet_path} after running 05_regression_healpix.py, "
        "but it is missing. The parquet write must have failed."
    )
print(f"\n[ok] {parquet_path.name} present "
      f"({parquet_path.stat().st_size:,} bytes)")

upstream_idata = OUT_DIR / "posterior.nc"
results_idata = RESULTS_DIR / "posterior_bambi_healpix.nc"
if bambi_ok and upstream_idata.exists():
    shutil.copy2(upstream_idata, results_idata)
    print(f"[ok] copied bambi idata -> {results_idata} "
          f"({results_idata.stat().st_size:,} bytes)")
else:
    print("[warn] bambi MCMC did not produce posterior.nc -- "
          "falling back to VB-only headline.")

# %% [markdown]
# ## statsmodels variational-Bayes GLMM on HEALPix

# %%
run("05b_regression_statsmodels_healpix.py")

# %% [markdown]
# ## Load posterior summary + extract `sc_TEI_delta`

# %%
posterior_csv = OUT_DIR / "posterior_vb_summary.csv"
post = pd.read_csv(posterior_csv, index_col=0)
print("\nPosterior summary (HEALPix VB):")
print(post.round(4).to_string())

row = post.loc["sc_TEI_delta"]
hp_mean = float(row["mean"])
hp_sd = float(row["sd"])
hp_p = float(row["p_2sided"])
hp_ci_low = hp_mean - 1.96 * hp_sd
hp_ci_high = hp_mean + 1.96 * hp_sd

# weatherxbio v0.2.1 published CEA values (external reference)
CEA_MEAN = 0.479
CEA_SD = 0.109
CEA_CI_LOW = CEA_MEAN - 1.96 * CEA_SD
CEA_CI_HIGH = CEA_MEAN + 1.96 * CEA_SD

# %% [markdown]
# ## Substrate-robustness tolerance test

# %%
TOL_REL = 0.30                                  # +-30% relative magnitude band
mag_lo = (1 - TOL_REL) * CEA_MEAN               # 0.7 * 0.479 = 0.3353
mag_hi = (1 + TOL_REL) * CEA_MEAN               # 1.3 * 0.479 = 0.6227

sign_positive = bool(hp_mean > 0)
significant_p_lt_05 = bool(hp_p < 0.05)
magnitude_within = bool(mag_lo <= hp_mean <= mag_hi)
all_pass = bool(sign_positive and significant_p_lt_05 and magnitude_within)

verdict = "Substrate-robust" if all_pass else "Substrate-divergent"

if all_pass:
    comment = (
        f"HEALPix nside=64 NESTED fit recovers sc_TEI_delta = {hp_mean:+.3f} "
        f"(p = {hp_p:.2e}), within +-30% of weatherxbio v0.2.1's CEA value of "
        f"+0.479 -- mechanism is substrate-robust at matched ~92 km cell scale."
    )
else:
    failed = []
    if not sign_positive:
        failed.append("sign (HEALPix mean is non-positive)")
    if not significant_p_lt_05:
        failed.append(f"significance (p={hp_p:.3g} >= 0.05)")
    if not magnitude_within:
        failed.append(
            f"magnitude (HEALPix mean {hp_mean:+.3f} outside "
            f"[{mag_lo:+.3f}, {mag_hi:+.3f}])"
        )
    comment = (
        f"HEALPix nside=64 NESTED fit gives sc_TEI_delta = {hp_mean:+.3f} "
        f"(p = {hp_p:.2e}); failed condition(s): {'; '.join(failed)} -- "
        "Soroye's mechanism does not transfer cleanly to this substrate."
    )

# %% [markdown]
# ## MCMC posterior (if bambi succeeded)

# %%
mcmc_block = None
hp_hdi_low = hp_ci_low
hp_hdi_high = hp_ci_high
if results_idata.exists():
    import arviz as az
    idata = az.from_netcdf(results_idata)
    mcmc_summary = az.summary(
        idata, var_names=["sc_TEI_delta"], hdi_prob=0.95,
    )
    mcmc_mean = float(mcmc_summary.loc["sc_TEI_delta", "mean"])
    mcmc_sd = float(mcmc_summary.loc["sc_TEI_delta", "sd"])
    mcmc_hdi_low = float(mcmc_summary.loc["sc_TEI_delta", "hdi_2.5%"])
    mcmc_hdi_high = float(mcmc_summary.loc["sc_TEI_delta", "hdi_97.5%"])
    sd_inflation = mcmc_sd / hp_sd if hp_sd else float("nan")
    mcmc_block = {
        "sc_TEI_delta_mean": mcmc_mean,
        "sc_TEI_delta_sd": mcmc_sd,
        "sc_TEI_delta_hdi95_low": mcmc_hdi_low,
        "sc_TEI_delta_hdi95_high": mcmc_hdi_high,
        "vb_sd_underestimation_factor": sd_inflation,
        "comment": (
            "MCMC posterior on HEALPix nside=64 substrate -- VB sd is a "
            "lower bound on true posterior uncertainty."
        ),
    }
    # Use MCMC HDI for the JSON's hdi95 fields when available
    hp_hdi_low = mcmc_hdi_low
    hp_hdi_high = mcmc_hdi_high
    print(f"\nMCMC sc_TEI_delta: {mcmc_mean:+.3f} "
          f"[{mcmc_hdi_low:+.3f}, {mcmc_hdi_high:+.3f}]  "
          f"sd ratio MCMC/VB = {sd_inflation:.2f}")

# %% [markdown]
# ## Iberia HEALPix substrate sanity counts (n_cells, n_species)

# %%
pa = np.load(OUT_DIR / "presence_absence_healpix.npz", allow_pickle=True)
n_cells_iberia = int(pa["n_cells"][0])
n_species = int(len(pa["species"]))

data_ext = pd.read_parquet(parquet_path)
n_cells_in_fit = int(data_ext["site"].nunique())

# %% [markdown]
# ## Write `results/headline_statistic_healpix.json`

# %%
report = {
    "healpix_fit": {
        "nside": 64,
        "depth": 6,
        "scheme": "NESTED",
        "n_cells_iberia": n_cells_iberia,
        "n_cells_in_fit": n_cells_in_fit,
        "n_species": n_species,
        "sc_TEI_delta_mean": hp_mean,
        "sc_TEI_delta_sd": hp_sd,
        "p_2sided": hp_p,
        "hdi95_low": hp_hdi_low,
        "hdi95_high": hp_hdi_high,
        "marginal_R2": None,
    },
    "weatherxbio_v0_2_1_cea": {
        "sc_TEI_delta_mean": CEA_MEAN,
        "sc_TEI_delta_sd": CEA_SD,
        "ci95_low": CEA_CI_LOW,
        "ci95_high": CEA_CI_HIGH,
        "source": (
            "weatherxbio v0.2.1 / "
            "soroye_port/outputs_iberia/posterior_vb_summary.csv"
        ),
        "framing": (
            "external reference for substrate-robustness check "
            "(NOT this repo's local CEA rerun)"
        ),
    },
    "substrate_robustness": {
        "sign_positive": sign_positive,
        "significant_p_lt_05": significant_p_lt_05,
        "magnitude_within_30pct_of_cea": magnitude_within,
        "all_three_pass": all_pass,
        "tolerance_band_mean": [mag_lo, mag_hi],
        "verdict": verdict,
        "comment": comment,
    },
    "mcmc": mcmc_block,
}

out_path = RESULTS_DIR / "headline_statistic_healpix.json"
with open(out_path, "w") as f:
    json.dump(report, f, indent=2)
print(f"\nWrote {out_path}")
print(json.dumps(report, indent=2))

# Also persist the full posterior table inside results/
post.to_csv(RESULTS_DIR / "glmm_coefficients_healpix.csv")
print(f"Wrote {RESULTS_DIR / 'glmm_coefficients_healpix.csv'}")

# %% [markdown]
# ## One-line conclusion

# %%
print(
    f"\nHEALPix nside=64 sc_TEI_delta = {hp_mean:+.3f} "
    f"[{hp_ci_low:+.3f}, {hp_ci_high:+.3f}]  "
    f"vs CEA reference {CEA_MEAN:+.3f} -- verdict: {verdict}"
)
