# 05 — FORRT Replication Outcome

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Verify the actual numerical results first** by reading `results/` and `notebooks/03_analysis.py`. Don't quote numbers from memory. See `docs/verify-before-drafting.md`.

## Field-by-field draft

### Short URI suffix for outcome ID (text input, required)

Slug. Use kebab-case.

```

```

### Plain-text label for the outcome (text input, required)

Descriptive title.

```

```

### Search for a FORRT replication study (search/select, required)

URI of the Replication Study published in step 04. Pull from `nanopubs/PUBLISHED.md`.

```

```

### Repository URL (text input, required)

```
https://github.com/annefou/weatherxbiodiversity-projection
```

### Completion date (date picker, required)

```
2026-05-09
```

### Validation status (dropdown, required)

- [ ] Validated
- [ ] PartiallySupported
- [ ] Contradicted

This dropdown maps to the CiTO intention in step 06: Validated → `confirms`, PartiallySupported → `qualifies`, Contradicted → `disputes`.

### Confidence level (dropdown, required)

_Vocabulary not yet captured._

```

```

### Describe the overall conclusion about the original claim (textarea, required)

Substantive interpretation. Headline comparison: replication's number vs the paper's number, sign + significance.

```

```

### Describe the evidence that supports your conclusion (textarea, required)

Numerical results, test statistics, model coefficients. Read directly from `results/`.

```

```

### Describe what limits the conclusions of the study (textarea, optional)

Honest caveats. If the result is partial or contradicted, say so plainly. Don't overclaim.

```
While Soroye 2020's TEI→extirpation mechanism is substrate-robust between sphere-HEALPix and WGS84-ellipsoidal HEALPix at matched cell scale (Tier 1 Outcome [URI]), per-species ranking near the logit saturation boundary (P > 0.7) is sensitive to grid choice. A 9% change in cell coverage between substrates shifted the top-3 species list. Practitioners should treat absolute-probability rankings near saturation with caution and prefer substrate-stable summary statistics (e.g. relative risk, count of species above a moderate threshold) for prioritisation decisions.

The Tier-2 projection therefore reports the GLMM linear predictor η (log-odds of extirpation) rather than its logistic transform p = expit(η). Future TEI/PEI z-scores under SSP3-7.0 lie 5–10× outside the Tier-1 training distribution (the 2020–2039 vs 1901–1974 warming signal is roughly twice the 2000–2014 vs 1901–1974 signal the GLMM was calibrated against), driving most species' linear predictors past η = +5 in many cells where expit(η) saturates uninformatively at ≈ 1.0. The η-based ranking preserves the GLMM's authentic signal: most-vulnerable species under SSP3-7.0 are cold-adapted, range-edge taxa (B. norvegicus, B. mendax, B. mucidus, B. pyrenaeus); B. terrestris and B. pascuorum show negative projected η (slight decrease in extirpation tendency), driven by precipitation predictors. Absolute extirpation probabilities derived from this projection are not interpretable as-is; relative ranking by η, fraction of cells with η > 0, and direction-of-effect are the substrate-stable summaries.
```

> **Drafter notes — additional limitations to weave into the textarea above when finalising this Outcome.** These are NOT to be pasted as-is; condense as needed.
>
> - **VB underestimation of posterior variance** — the statsmodels variational-Bayes posterior SD is ~1.4× smaller than the bambi/PyMC NUTS posterior SD on the same data; report MCMC HDIs (not VB CIs) wherever uncertainty is consequential.
> - **DestinE Climate DT archive coverage** — populated through 2039 only at time of analysis; mid-/end-of-century horizons (2046–2055, 2076–2085) deferred to a follow-up Outcome when the archive extends past 2050.
> - **Sampling effort held at recent-period mean** for the projection — assumes monitoring intensity in 2020-2039 mirrors the 2000-2014 baseline; biases the projection if survey effort actually changes.
> - **Daily Tmax/Tmin from 4-times-daily 2t snapshots** — 06_destine_clean approximates daily extremes from the 00/06/12/18 UTC samples rather than from true daily max/min (which DestinE Climate DT does not archive — see notebooks/05_destine_download.py § Why no native max/min).
> - **`tp` 1-time-per-day approximation in the DestinE retrieve** — flagged in 06_destine_clean as TP_HOURLY_TO_DAILY_FACTOR; defensible at decadal mean, invalid for any sub-monthly statistic.
> - **Small-N in the most-vulnerable taxa** — top-ranked species (norvegicus, mendax, mucidus, pyrenaeus) are observed in 1–7 Iberian cells in the 1901–1974 baseline period; high η reflects strong climate forcing at those few cells, not breadth of impact across the species' range. Treat the per-species η ranking jointly with `n_cells` and `n_cells_eta_gt_0` from `results/projection_headline.json`.
> - **Negative-η species reflect random intercepts, not climate-driven benefit** — the GLMM's species random effect dominates the projected η for *B. terrestris* (RE ≈ −2.6) and *B. pascuorum* (RE ≈ −0.97), encoding lower-than-average historical extirpation susceptibility from the 1901–2014 GBIF training data. **All climate-driven term contributions are POSITIVE for both species** (terrestris: +0.92 net across climate terms; pascuorum: +0.42). The correct reading is "historically robust species, with SSP3-7.0 adding moderate climate forcing on top" — NOT "projected to benefit from SSP3-7.0". Verified by per-term η decomposition in `scripts/decompose_eta_for_species.py`.

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 05.
