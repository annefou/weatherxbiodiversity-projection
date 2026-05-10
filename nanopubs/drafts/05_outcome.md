# 05 — FORRT Replication Outcome

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Verify the actual numerical results first** by reading `results/` and `notebooks/03_analysis.py`. Don't quote numbers from memory. See `docs/verify-before-drafting.md`.

## Field-by-field draft

### Short URI suffix for outcome ID (text input, required)

```
soroye2020-iberian-bombus-substrate-robust-validated
```

### Plain-text label for the outcome (text input, required)

```
Soroye et al. 2020 TEI-based extirpation mechanism replicates on Iberian Bombus (CEA + HEALPix nside=64; SSP3-7.0 projection)
```

### Search for a FORRT replication study (search/select, required)

URI of the Replication Study published in step 04. Pull from `nanopubs/PUBLISHED.md` after step 04 is live.

```
<replace-with-published-Study-URI-from-step-04>
```

### Repository URL (text input, required)

```
https://github.com/annefou/weatherxbiodiversity-projection
```

### Completion date (date picker, required)

```
2026-05-10
```

### Validation status (dropdown, required)

- [x] **Validated**
- [ ] PartiallySupported
- [ ] Contradicted

Maps to CiTO `confirms` in step 06. Soroye's central biological claim — that thermal-niche exceedance increases extirpation probability — replicates on Iberian *Bombus* using Soroye's own CRU TS 3.24.01 climate inputs and own-issued GBIF download DOIs, at two independent spatial substrates within ±30%.

### Confidence level (dropdown, required)

```
High
```

(Posterior 95% HDI on the headline coefficient excludes zero at both substrates; same sign and order of magnitude as the original; confirmed independently by two GLMM fitting strategies — variational-Bayes via statsmodels `BinomialBayesMixedGLM` and full-posterior NUTS via bambi/PyMC.)

### Describe the overall conclusion about the original claim (textarea, required)

```
Soroye et al.'s (2020) headline claim — that thermal-niche exceedance (TEI) increases extirpation probability for bumble bees — REPLICATES on Iberian Bombus.

The GLMM coefficient on standardised TEI_delta is positive, large, and credible at both spatial substrates tested in this replication: +0.479 at the original CEA grid (~100 km cells) and +0.454 at HEALPix nside=64 (~92 km cells, 95% HDI [+0.130, +0.751]). Both estimates are within ±30 percent of each other, same sign, same order of magnitude as Soroye 2020's continental-scale fit. A higher-resolution sibling study at HEALPix nside=128 (weatherxbiodiversity-projection-nside128) finds +0.347, HDI [+0.139, +0.533] — confirming substrate-robustness across three independent pixelisations.

The Tier-2 SSP3-7.0 projection (DestinE Climate DT, 2030–2039 horizon) shows a systematic shift in mean future TEI from ~0.43 historically to ~0.45–0.50, indicating that Iberian Bombus species are drifting toward their upper thermal niche edge. None of the studied species hit TEI > 1 (Soroye's actual extirpation threshold) at any currently-occupied cell over the next 15 years — the extirpation event Soroye projects is a longer-timescale phenomenon, with 2030–2039 being the early-warning period rather than the extirpation period itself.

When the per-species ranking is filtered to species with at least 10 occupied cells per substrate AND uses the GLMM main-effects-only η at projection time (per the substrate-sensitivity diagnostic in the sibling repo weatherxbiodiversity-substrate-sensitivity, ρ_Spearman = +0.97 between substrates at this filter), the highest-risk Iberian species are B. humilis, B. muscorum, and B. ruderarius — short-tongued grassland specialists already documented in decline elsewhere in Europe. The lowest-risk species are B. terrestris, B. pascuorum, and B. pratorum — the dominant Western Palearctic generalists. The substrate-stable Iberian projection independently identifies the same European high-risk species the long-term monitoring data flag.
```

### Describe the evidence that supports your conclusion (textarea, required)

```
TIER 1 — historical fit (1901–2014 CRU TS 3.24.01 from Soroye Figshare; GBIF Iberian Bombus occurrences via own download DOI).

GLMM headline coefficient sc_TEI_delta (standardised TEI change between baseline 1901–1974 and recent 2000–2014):
  - CEA (~100 km cells):       +0.479    (this run, replicates published v0.2.1)
  - HEALPix nside=64 (~92 km): +0.454    (95% HDI [+0.130, +0.751])
  - HEALPix nside=128 (~46 km): +0.347   (95% HDI [+0.139, +0.533]; sibling repo)
  All three within ±30 percent. Sign, magnitude, and HDI agree with Soroye 2020.

Cross-fit sanity check: variational-Bayes via statsmodels BinomialBayesMixedGLM agrees with full-posterior NUTS via bambi/PyMC on β_sc_TEI_delta to within VB-underestimation noise (VB posterior SD ~1.4× smaller than NUTS posterior SD on the same data; report MCMC HDIs as authoritative). See results/headline_statistic_healpix.json.

TIER 2 — SSP3-7.0 projection (DestinE Climate DT, IFS-NEMO standard, native HEALPix nside=128, 2020–2029 + 2030–2039 horizons).

Mean future TEI per species, 2030–2039 horizon (substrate-invariant physical metric, Spearman ρ between substrates = +0.66 at n≥1, +0.90 at n≥10):
  - Most-shifted species: B. humilis (mean fut TEI ~0.50), B. muscorum (~0.49), B. ruderarius (~0.48)
  - Least-shifted: B. terrestris (~0.45), B. pascuorum (~0.45), B. pratorum (~0.45)

Per-species community-mean η under SSP3-7.0 (filtered to species with n_cells≥10 at both substrates, main-effects-only η to avoid extrapolation amplification — ρ_Spearman between substrates = +0.97 mid-term, +0.96 near-term):
  - Highest projected η (most at-risk): B. humilis (η = +1.58), B. muscorum (+1.50), B. ruderarius (+1.54)
  - Lowest projected η (most robust): B. terrestris (-2.26), B. pascuorum (-0.69), B. pratorum (-0.56)

Fraction of currently-occupied cells where future TEI exceeds 1.0 (Soroye's actual extirpation threshold) at 2030–2039: 0.00 for every studied species — the extirpation threshold is not yet crossed under SSP3-7.0 over the next 15 years.

Files: results/headline_statistic.json (CEA), results/headline_statistic_healpix.json (HEALPix nside=64), results/projection_headline.json (Tier-2 per-species rank), results/posterior_bambi.nc + posterior_bambi_healpix.nc (full posteriors); figures figures/main_result_healpix.png + figures/projection_summary.png.
```

### Describe what limits the conclusions of the study (textarea, optional)

```
1. Substrate-coupling at projection time. Per-species ranking under SSP3-7.0 is grid-coupled at the projection step (NOT at fit time) for species observed in fewer than ~10 historical Iberian cells. This affects 18 of the 31 species in the study, including narrowly-distributed Pyrenean specialists (B. pyrenaeus, B. mucidus, B. mendax, B. wurflenii, B. monticola). The full mechanism (interaction-term amplification of substrate-specific predictor standardisation at extrapolation), validated and refuted hypotheses, and recommended reporting protocol are documented in the methodological sibling repo weatherxbiodiversity-substrate-sensitivity. The figures and rankings reported in this Outcome are filtered per that diagnostic (n_cells ≥ 10, main-effects-only η).

2. Reporting on η, not p_extirp. The Tier-2 projection reports the GLMM linear predictor η (log-odds of extirpation) rather than its logistic transform p = expit(η). Future TEI/PEI z-scores under SSP3-7.0 lie 2–4σ outside the Tier-1 training distribution; expit(η) saturates uninformatively at ~1.0 in many cells where the actual signal is the GLMM's authentic relative ordering. Absolute extirpation probabilities derived from this projection are NOT interpretable as-is; relative ranking by η and substrate-invariant physical metrics (mean future TEI, fraction TEI_future>0.5) are the substrate-stable summaries.

3. Negative-η species reflect random intercepts, not climate-driven benefit. The GLMM's species random effect dominates the projected η for B. terrestris (RE ~ -2.6) and B. pascuorum (RE ~ -1.0), encoding lower-than-average historical extirpation susceptibility from the 1901–2014 GBIF training data. All climate-driven term contributions are POSITIVE for both species. The correct reading is "historically robust species, with SSP3-7.0 adding moderate climate forcing on top" — NOT "projected to benefit from SSP3-7.0".

4. Sampling effort held at recent-period mean for the projection — assumes monitoring intensity in 2020–2039 mirrors the 2000–2014 baseline. Biases the projection if survey effort actually changes.

5. DestinE Climate DT archive coverage. SSP3-7.0 IFS-NEMO archive populated through 2039 only at time of analysis; mid- and end-of-century horizons (2046–2055, 2076–2085) deferred to a follow-up Outcome when the archive extends past 2050.

6. Daily Tmax/Tmin from 4-times-daily 2t snapshots. Notebook 06_destine_clean approximates daily extremes from the 00/06/12/18 UTC samples rather than from true daily max/min (which DestinE Climate DT does not archive — see notebooks/05_destine_download.py § Why no native max/min). Defensible at decadal mean; biased high on Tmax extremes by an unknown but small amount.

7. tp 1-time-per-day approximation in the DestinE retrieve — flagged in 06_destine_clean as TP_HOURLY_TO_DAILY_FACTOR; defensible at decadal mean, invalid for any sub-monthly statistic.

8. VB underestimation of posterior variance. statsmodels variational-Bayes posterior SD is ~1.4× smaller than the bambi/PyMC NUTS posterior SD on the same data; this Outcome's HDIs are NUTS-based.

9. Geographic scope. Iberian peninsula only. Extension to other latitudes (especially boreal/alpine systems where Bombus diversity peaks) is not addressed by this replication.
```

> **Drafter notes — additional context for finalising this Outcome.**
>
> - Validation Status = `Validated` because Soroye's TEI coefficient REPLICATES at three substrates within ±30%. The substrate-coupling at projection time is a downstream methodological observation (covered by the substrate-sensitivity sibling chain), not a refutation of the original mechanism.
> - The CiTO intention in step 06 is therefore `cito:confirms` Soroye 2020, plus an additional `cito:extends` to the substrate-sensitivity sibling chain that documents the projection-time caveat.

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 05.
