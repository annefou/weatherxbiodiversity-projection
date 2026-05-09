# HEALPix-NESTED nside=64 port of soroye_port/05b_regression_statsmodels.py — substrate-robustness branch.
"""
HEALPix-NESTED nside=64 port of script 05b.

Variational-Bayes binomial mixed-effects GLMM via statsmodels'
`BinomialBayesMixedGLM.fit_vb()`. Reads `dataGLMM_extinction.parquet`
(written by 05_regression_healpix.py) and writes
`posterior_vb_summary.csv` -- a compact mean / sd / z / p_2sided table
that the analysis notebook reads to compute the headline
`sc_TEI_delta`. Same approach the upstream pipeline uses to fall
through bambi when pytensor / CLT issues block MCMC on macOS.
"""
from __future__ import annotations

from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / 'healpix_port' / 'outputs_iberia'

data_ext = pd.read_parquet(OUT_DIR / 'dataGLMM_extinction.parquet')
print(f'Data: {len(data_ext):,} rows, '
      f'{data_ext["species"].nunique()} species, '
      f'{data_ext["site"].nunique()} cells')
print(f'Extinction rate: {data_ext["extinction"].mean():.3f}')

FORMULA = (
    'extinction ~ continent + sc_sampling'
    ' + sc_TEI_bs + sc_TEI_delta + sc_TEI_bs:sc_TEI_delta'
    ' + sc_PEI_bs + sc_PEI_delta + sc_PEI_bs:sc_PEI_delta'
    ' + sc_TEI_bs:sc_PEI_bs + sc_TEI_delta:sc_PEI_delta'
)
if data_ext['continent'].nunique() < 2:
    print('  continent is constant -- dropping from formula')
    FORMULA = FORMULA.replace('continent + ', '')

# ---------------------------------------------------------------------------
# 1. Plain logistic regression (no random effect) -- reference / sanity check.

print('\n=== 1) Plain logistic regression (no species random effect) ===')
try:
    glm = smf.logit(FORMULA, data=data_ext).fit(disp=False, maxiter=200)
    print(glm.summary().tables[1])
except Exception as e:
    print(f'  Plain logistic failed: {e}')
    glm = None

# ---------------------------------------------------------------------------
# 2. BinomialBayesMixedGLM (variational Bayes) with species random effect.

print('\n=== 2) statsmodels BinomialBayesMixedGLM (species random effect) ===')
try:
    from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

    md0 = smf.glm(FORMULA, data=data_ext, family=sm.families.Binomial())
    exog_names = md0.exog_names
    print(f'  Fixed effects ({len(exog_names)}): {exog_names}')

    md = BinomialBayesMixedGLM.from_formula(
        FORMULA,
        vc_formulas={'species': '0 + C(species)'},
        data=data_ext,
    )
    print('  Fitting variational Bayes (fast approximation) ...')
    vb_result = md.fit_vb()
    coef_df = pd.DataFrame({
        'mean': vb_result.fe_mean,
        'sd':   vb_result.fe_sd,
    }, index=exog_names)
    coef_df['z'] = coef_df['mean'] / coef_df['sd']
    from scipy.stats import norm
    coef_df['p_2sided'] = 2 * (1 - norm.cdf(coef_df['z'].abs()))
    print(coef_df.round(4).to_string())
    coef_df.to_csv(OUT_DIR / 'posterior_vb_summary.csv')
except Exception as e:
    import traceback
    traceback.print_exc()
    print(f'  Mixed Bayesian failed: {e}')
    vb_result = None

# ---------------------------------------------------------------------------
# 3. Key coefficient report.

if glm is not None:
    coef = glm.params.get('sc_TEI_delta', np.nan)
    pval = glm.pvalues.get('sc_TEI_delta', np.nan)
    print('\n=== PLAIN LOGIT -- sc_TEI_delta ===')
    print(f'  coef: {coef:+.4f}  p = {pval:.4g}')

if vb_result is not None:
    fe_names = list(md.exog_names)
    idx = fe_names.index('sc_TEI_delta')
    m = vb_result.fe_mean[idx]
    s = vb_result.fe_sd[idx]
    print('\n=== VB MIXED MODEL -- sc_TEI_delta ===')
    print(f'  posterior mean: {m:+.4f}   sd: {s:.4f}')
    print(f'  95% approx CI: [{m - 1.96 * s:.4f}, {m + 1.96 * s:.4f}]')

print('\nDone.')
