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
```

> **Drafter notes — additional limitations to weave into the textarea above when finalising this Outcome.** These are NOT to be pasted as-is; condense as needed.
>
> - **VB underestimation of posterior variance** — the statsmodels variational-Bayes posterior SD is ~1.4× smaller than the bambi/PyMC NUTS posterior SD on the same data; report MCMC HDIs (not VB CIs) wherever uncertainty is consequential.
> - **DestinE Climate DT archive coverage** — populated through 2039 only at time of analysis; mid-/end-of-century horizons (2046–2055, 2076–2085) deferred to a follow-up Outcome when the archive extends past 2050.
> - **Sampling effort held at recent-period mean** for the projection — assumes monitoring intensity in 2020-2039 mirrors the 2000-2014 baseline; biases the projection if survey effort actually changes.
> - **Daily Tmax/Tmin from 4-times-daily 2t snapshots** — 06_destine_clean approximates daily extremes from the 00/06/12/18 UTC samples rather than from true daily max/min (which DestinE Climate DT does not archive — see notebooks/05_destine_download.py § Why no native max/min).
> - **`tp` 1-time-per-day approximation in the DestinE retrieve** — flagged in 06_destine_clean as TP_HOURLY_TO_DAILY_FACTOR; defensible at decadal mean, invalid for any sub-monthly statistic.

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 05.
