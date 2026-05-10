# HEALPix-NESTED nside=64 port of soroye_port/05_regression.py — substrate-robustness branch.
"""
HEALPix-NESTED nside=64 port of script 05.

Same Bayesian binomial GLMM as the upstream / CEA pipeline, fitted on
the HEALPix substrate. Formula:

    extinction ~ continent + sc_sampling
      + sc_TEI_bs + sc_TEI_delta + sc_TEI_bs:sc_TEI_delta
      + sc_PEI_bs + sc_PEI_delta + sc_PEI_bs:sc_PEI_delta
      + sc_TEI_bs:sc_PEI_bs + sc_TEI_delta:sc_PEI_delta
      + (1|species)

If `continent` is constant in the subset (always the case for the
Iberia-only run), it is dropped from the formula to keep the design
matrix non-singular -- same as the upstream.

bambi/PyMC NUTS: 2 chains x 2000 draws + 1000 tune. If the bambi /
pymc / pytensor stack is unavailable on this host, the script logs the
failure and exits 0 -- the downstream statsmodels variational-Bayes
fit (script 05b) is enough to compute the headline `sc_TEI_delta`.
"""
from __future__ import annotations

from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / 'healpix_port' / 'outputs_iberia'

# ---------------------------------------------------------------------------
# 1. Load intermediates from scripts 02, 03, 04.

print('Loading HEALPix intermediates ...')
pa = xr.open_dataset(OUT_DIR / 'presence_absence_healpix.nc')
sc = xr.open_dataset(OUT_DIR / 'sampling_continent_healpix.nc')
cl = xr.open_dataset(OUT_DIR / 'climate_tei_pei_healpix.nc')

species_list = [str(s) for s in pa['species'].values]
prab_baseline = pa['prab_baseline'].values
prab_recent = pa['prab_recent'].values

sampling_baseline = sc['sampling_baseline'].values
sampling_recent = sc['sampling_recent'].values
sampling_total = sc['sampling_total'].values
# Continent is int8 with -1 as missing; convert to float w/ NaN to
# preserve the upstream isnan-based gating.
continent_int = sc['continent'].values.astype(np.int16)
continent = np.where(continent_int < 0, np.nan, continent_int).astype(np.float32)
n_cells = int(pa.sizes['cells'])

avgtemp_bs = cl['avgtemp_bs'].values
avgtemp_delta = cl['avgtemp_delta'].values
avgprecip_bs = cl['avgprecip_bs'].values
avgprecip_delta = cl['avgprecip_delta'].values
TEI_bs = cl['tei_bs'].values
TEI_delta = cl['tei_delta'].values
PEI_bs = cl['pei_bs'].values
PEI_delta = cl['pei_delta'].values

n_spp = prab_baseline.shape[0]
print(f'  {n_spp} species x {n_cells} cells')

# ---------------------------------------------------------------------------
# 2. Build dataGLMM -- one row per (species, cell).

print('\nAssembling dataGLMM ...')
rows = []
for s, sp in enumerate(species_list):
    bs = prab_baseline[s]
    rc = prab_recent[s]
    with np.errstate(invalid='ignore'):
        pr_change = 2 * rc - bs
        extinction = np.where(
            pr_change == -1, 1.0,
            np.where((pr_change == 2) | (pr_change == 1), 0.0, np.nan),
        )
        colonization = np.where(
            pr_change == 2, 1.0,
            np.where((pr_change == -1) | (pr_change == 1), 0.0, np.nan),
        )
    for c in range(n_cells):
        if np.isnan(extinction[c]) and np.isnan(colonization[c]):
            continue
        if np.isnan(continent[c]) or np.isnan(sampling_total[c]):
            continue
        if np.isnan(TEI_bs[s, c]) or np.isnan(TEI_delta[s, c]):
            continue
        if np.isnan(PEI_bs[s, c]) or np.isnan(PEI_delta[s, c]):
            continue
        rows.append({
            'species': sp,
            'site': c,
            'extinction': extinction[c],
            'colonization': colonization[c],
            'continent': int(continent[c]),
            'sampling': sampling_total[c],
            'TEI_bs': TEI_bs[s, c],
            'TEI_delta': TEI_delta[s, c],
            'PEI_bs': PEI_bs[s, c],
            'PEI_delta': PEI_delta[s, c],
            'avgtemp_bs': avgtemp_bs[c],
            'avgtemp_delta': avgtemp_delta[c],
            'avgprecip_bs': avgprecip_bs[c],
            'avgprecip_delta': avgprecip_delta[c],
        })

dataGLMM = pd.DataFrame(rows)
print(f'  dataGLMM: {len(dataGLMM):,} rows')

# ---------------------------------------------------------------------------
# 3. Z-score scaling (ddof=1, matches R's scale()).


def z(col: str) -> pd.Series:
    m = dataGLMM[col].mean()
    s = dataGLMM[col].std(ddof=1)
    return (dataGLMM[col] - m) / s


dataGLMM['sc_sampling'] = z('sampling')
dataGLMM['sc_TEI_bs'] = z('TEI_bs')
dataGLMM['sc_TEI_delta'] = z('TEI_delta')
dataGLMM['sc_PEI_bs'] = z('PEI_bs')
dataGLMM['sc_PEI_delta'] = z('PEI_delta')
dataGLMM['sc_avgtemp_bs'] = z('avgtemp_bs')
dataGLMM['sc_avgtemp_delta'] = z('avgtemp_delta')
dataGLMM['sc_avgprecip_bs'] = z('avgprecip_bs')
dataGLMM['sc_avgprecip_delta'] = z('avgprecip_delta')

dataGLMM['continent'] = dataGLMM['continent'].map({1: 'North America', 2: 'Europe'})
dataGLMM['continent'] = pd.Categorical(
    dataGLMM['continent'], categories=['North America', 'Europe'],
)

data_ext = dataGLMM[dataGLMM['extinction'].notna()].copy()
print(f'  extinction subset: {len(data_ext):,} rows')
print(f'    extinctions: {int((data_ext.extinction == 1).sum()):,}')
print(f'    persistences: {int((data_ext.extinction == 0).sum()):,}')
print(f'    species: {data_ext.species.nunique()}, '
      f'cells: {data_ext.site.nunique()}')

data_ext.to_parquet(OUT_DIR / 'dataGLMM_extinction.parquet')

# ---------------------------------------------------------------------------
# 4. Bayesian GLMM via bambi (optional -- statsmodels VB fallback in 05b).

print('\nFitting Bayesian GLMM via bambi ...')
try:
    import bambi as bmb

    formula = (
        'extinction ~ continent + sc_sampling'
        ' + sc_TEI_bs + sc_TEI_delta + sc_TEI_bs:sc_TEI_delta'
        ' + sc_PEI_bs + sc_PEI_delta + sc_PEI_bs:sc_PEI_delta'
        ' + sc_TEI_bs:sc_PEI_bs + sc_TEI_delta:sc_PEI_delta'
        ' + (1|species)'
    )
    if data_ext['continent'].nunique() < 2:
        print('  continent is constant -- dropping from formula')
        formula = formula.replace('continent + ', '')

    model = bmb.Model(formula, data_ext, family='bernoulli')
    print(model)

    idata = model.fit(
        draws=2000, tune=1000, chains=2, target_accept=0.95,
        progressbar=True, idata_kwargs={'log_likelihood': False},
    )

    import arviz as az
    summary = az.summary(idata, hdi_prob=0.95)
    print('\n=== Posterior summary (bambi MCMC) ===')
    print(summary.to_string())

    summary.to_csv(OUT_DIR / 'posterior_summary.csv')

    if 'sc_TEI_delta' in summary.index:
        key_row = summary.loc['sc_TEI_delta']
        print('\n=== KEY TARGET COEFFICIENT ===')
        print(f'sc_TEI_delta posterior mean: {key_row["mean"]:.4f}')
        print(f'  95% HDI: [{key_row["hdi_2.5%"]:.4f}, {key_row["hdi_97.5%"]:.4f}]')
        print(f'  sign: {"POSITIVE" if key_row["mean"] > 0 else "NEGATIVE"}')

    idata.to_netcdf(OUT_DIR / 'posterior.nc')
    print(f'\nSaved bambi posterior summary + netcdf to {OUT_DIR}')
except Exception as e:
    import traceback
    traceback.print_exc()
    print(f'\n[warn] bambi MCMC failed: {e}')
    print('Falling back to statsmodels variational-Bayes fit '
          '(see 05b_regression_statsmodels_healpix.py).')
