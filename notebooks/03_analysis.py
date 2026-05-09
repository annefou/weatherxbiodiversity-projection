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
# # 03 — Analysis (Pass 1, Iberian baseline)
#
# Runs the vendored upstream **statsmodels variational-Bayes GLMM** on
# the Iberia presence/absence + climate stack produced by
# `02_data_clean.py`, then writes a comparison JSON
# (`results/headline_statistic.json`) that lays the re-run side-by-side
# with the upstream-v0.2.1 published numbers.
#
# The headline coefficient under test is `sc_TEI_delta` — the change in
# Thermal Exposure Index, which Soroye et al. (2020) identify as the
# strongest predictor of local extinction. The upstream-v0.2.1
# replication for Iberia reported
#
# > `sc_TEI_delta = +0.479 [0.265, 0.694]`  (mean ± 95 % CI from VB SD)
#
# A faithful re-run should land within 1.96 sigma of that value (the
# stochasticity of variational Bayes is small, but `statsmodels`
# updates may shift the posterior slightly between releases).

# %%
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# %%
ROOT = Path("..").resolve()
PORT = ROOT / "soroye_port"
OUT_DIR = PORT / "outputs_iberia"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

env = {**os.environ, "OUT_SUBDIR": "outputs_iberia"}


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
# ## Bayesian GLMM (bambi / PyMC) — full posterior covariance
#
# `05_regression.py` assembles the `dataGLMM_extinction.parquet` table
# from the three npz intermediates (presence/absence, sampling effort,
# climate TEI/PEI) at line ~158, then runs a `bambi` / `pymc` MCMC fit
# (NUTS, 2 chains × 2000 draws + 1000 tune) and writes the full
# posterior to `posterior.nc`. Tier 2 (DestinE projection) samples
# from this posterior to propagate coefficient uncertainty through the
# future-climate projection.
#
# The upstream noted macOS pytensor / CLT issues; on current
# pymc/pytensor (≥ 5.x) this is usually resolved. If the fit fails,
# the wrapper falls through and the statsmodels VB run downstream still
# produces a usable (but lower-fidelity) posterior summary.

# %%
import shutil

bambi_ok = run("05_regression.py", check=False) == 0

parquet_path = OUT_DIR / "dataGLMM_extinction.parquet"
if not parquet_path.exists():
    raise SystemExit(
        f"Expected {parquet_path} after running 05_regression.py, but it "
        "is missing. The parquet write at line ~158 must have failed."
    )
print(f"\n[ok] {parquet_path.name} present "
      f"({parquet_path.stat().st_size:,} bytes)")

# Promote the bambi posterior NetCDF into results/ so Tier 2 finds it
# at a stable, committable path. soroye_port/outputs_iberia/ is
# gitignored; results/posterior_bambi.nc is tracked.
upstream_idata = OUT_DIR / "posterior.nc"
results_idata = RESULTS_DIR / "posterior_bambi.nc"
if bambi_ok and upstream_idata.exists():
    shutil.copy2(upstream_idata, results_idata)
    print(f"[ok] copied bambi idata → {results_idata} "
          f"({results_idata.stat().st_size:,} bytes)")
else:
    print("[warn] bambi MCMC did not produce posterior.nc — "
          "Tier 2 projection will need a successful Tier 1 bambi fit.")

# %% [markdown]
# ## Run the statsmodels variational-Bayes regression
#
# `05b_regression_statsmodels.py` reads the parquet and writes
# `soroye_port/outputs_iberia/posterior_vb_summary.csv` with the
# variational-Bayes posterior means, standard deviations, z-scores and
# 2-sided p-values for every fixed effect.

# %%
run("05b_regression_statsmodels.py")

# %% [markdown]
# ## Load posterior summary + extract `sc_TEI_delta`

# %%
posterior_csv = OUT_DIR / "posterior_vb_summary.csv"
post = pd.read_csv(posterior_csv, index_col=0)
print("\nPosterior summary:")
print(post.round(4).to_string())

row = post.loc["sc_TEI_delta"]
rep_mean = float(row["mean"])
rep_sd = float(row["sd"])
rep_p = float(row["p_2sided"])
rep_ci_low = rep_mean - 1.96 * rep_sd
rep_ci_high = rep_mean + 1.96 * rep_sd

# Upstream v0.2.1 published numbers (from
# soroye_port/outputs_iberia/posterior_vb_summary.csv at tag v0.2.1).
upstream_mean = 0.479
upstream_sd = 0.109
upstream_ci_low = upstream_mean - 1.96 * upstream_sd
upstream_ci_high = upstream_mean + 1.96 * upstream_sd

# Two-sample z-difference between independent posterior means.
z_diff = (rep_mean - upstream_mean) / np.sqrt(rep_sd ** 2 + upstream_sd ** 2)
verdict = "Replicated" if abs(z_diff) < 1.96 else "Diverges"

# %% [markdown]
# ## Write the comparison JSON

# %%
# MCMC marginal for sc_TEI_delta (full posterior covariance available
# in results/posterior_bambi.nc for Tier 2). VB underestimates posterior
# variance — this comparison documents the magnitude.
mcmc_block = None
if results_idata.exists():
    import arviz as az
    idata = az.from_netcdf(results_idata)
    mcmc_summary = az.summary(idata, var_names=["sc_TEI_delta"], hdi_prob=0.95)
    mcmc_mean = float(mcmc_summary.loc["sc_TEI_delta", "mean"])
    mcmc_sd = float(mcmc_summary.loc["sc_TEI_delta", "sd"])
    mcmc_hdi_low = float(mcmc_summary.loc["sc_TEI_delta", "hdi_2.5%"])
    mcmc_hdi_high = float(mcmc_summary.loc["sc_TEI_delta", "hdi_97.5%"])
    sd_inflation_factor = mcmc_sd / rep_sd if rep_sd else float("nan")
    mcmc_block = {
        "sc_TEI_delta_mean": mcmc_mean,
        "sc_TEI_delta_sd": mcmc_sd,
        "sc_TEI_delta_hdi95_low": mcmc_hdi_low,
        "sc_TEI_delta_hdi95_high": mcmc_hdi_high,
        "vb_sd_underestimation_factor": sd_inflation_factor,
        "comment": (
            "MCMC posterior SD vs statsmodels VB posterior SD — ratio > 1 "
            "means VB underestimates uncertainty. Tier 2 samples 1000 draws "
            "from this posterior."
        ),
    }
    print(f"\nMCMC sc_TEI_delta: {mcmc_mean:+.3f} "
          f"[{mcmc_hdi_low:+.3f}, {mcmc_hdi_high:+.3f}]  "
          f"sd ratio MCMC/VB = {sd_inflation_factor:.2f}")

report = {
    "replication": {
        "sc_TEI_delta_mean": rep_mean,
        "sc_TEI_delta_sd": rep_sd,
        "sc_TEI_delta_ci95_low": rep_ci_low,
        "sc_TEI_delta_ci95_high": rep_ci_high,
        "p_2sided": rep_p,
    },
    "upstream_v0_2_1": {
        "sc_TEI_delta_mean": upstream_mean,
        "sc_TEI_delta_sd": upstream_sd,
        "sc_TEI_delta_ci95_low": upstream_ci_low,
        "sc_TEI_delta_ci95_high": upstream_ci_high,
    },
    "agreement": {
        "z_diff": float(z_diff),
        "verdict": verdict,
    },
    "mcmc": mcmc_block,
}

out_path = RESULTS_DIR / "headline_statistic.json"
with open(out_path, "w") as f:
    json.dump(report, f, indent=2)
print(f"\nWrote {out_path}")
print(json.dumps(report, indent=2))

# Also persist the full posterior table inside results/ for the
# Jupyter Book / nanopub draft authors.
post.to_csv(RESULTS_DIR / "glmm_coefficients.csv")
print(f"Wrote {RESULTS_DIR / 'glmm_coefficients.csv'}")

# %% [markdown]
# ## Conclusion

# %%
print(
    f"\nsc_TEI_delta replication: {rep_mean:+.3f} "
    f"[{rep_ci_low:+.3f}, {rep_ci_high:+.3f}] "
    f"vs upstream {upstream_mean:+.3f} "
    f"[{upstream_ci_low:+.3f}, {upstream_ci_high:+.3f}] — "
    f"verdict: {verdict}"
)
