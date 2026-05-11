# 08 — Research Synthesis (optional)

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> Use this template only when this chain is **one of several** testing facets of a shared underlying property. The Synthesis names the cross-cutting conclusion and lists the multiple Outcomes as supporting sources.

**Form heading:** *"Science Live Research Synthesis — Synthesise findings across multiple replication outcomes with conclusions, recommendations, conditions, and limitations."*

## Field-by-field draft

### Short URI suffix for synthesis ID (text input, required)

Slug. Use kebab-case.

```
soroye2020-iberian-bombus-tei-mechanism-substrate-robust-three-substrates
```

### Label of the synthesis (text input, required)

A one-line summary.

```
Soroye et al. 2020 TEI-based extirpation mechanism is substrate-robust on Iberian Bombus across CEA, HEALPix nside=64, and HEALPix nside=128
```

### Conclusion of the synthesis (textarea, required)

The aggregate finding across the underlying outcomes.

```
Across three independent spatial pixelisations on Iberian Bombus — CEA (~100 km, this chain's Tier-1 Pass-1), HEALPix nside=64 (~92 km, this chain's Tier-1 Pass-2), and HEALPix nside=128 (~46 km, sibling chain weatherxbiodiversity-projection-nside128) — the GLMM coefficient on standardised TEI_delta is positive, large, and credibly above zero: +0.479 at CEA, +0.454 at nside=64 (95% HDI [+0.130, +0.751]), +0.347 at nside=128 (95% HDI [+0.139, +0.533]). All three estimates are within ±30% of each other, with the same sign and order of magnitude as Soroye et al.'s (2020) continental fit. Soroye's central biological claim — that thermal-niche exceedance increases extirpation probability — replicates on Iberian Bombus and is substrate-robust at the GLMM-fit step. The per-species PROJECTION ranking under SSP3-7.0 is grid-coupled for species observed in fewer than ~10 historical Iberian cells (see sibling chain weatherxbiodiversity-substrate-sensitivity); this synthesis is about the substrate-robustness of the FIT, not the projection.
```

### Recommendations (textarea, required)

Actionable guidance for practitioners.

```
1. Replicators of Soroye-style TEI-based extirpation models should expect substrate-robust headline coefficients within ±30% across factor-of-2 changes in cell resolution (HEALPix nside=64 ↔ nside=128) and across HEALPix-vs-CEA substrate choices. Substrate sensitivity at fit time is small.

2. Use the substrate that matches your projection forcing data's native pixelisation when possible (e.g. HEALPix nside=128 for DestinE Climate DT IFS-NEMO standard). This eliminates one source of cross-substrate aggregation noise between fit and projection.

3. When projecting to future climate, follow the protocol established in weatherxbiodiversity-substrate-sensitivity: filter species to those with ≥ 10 occupied + active cells per substrate, and drop the GLMM interaction terms at projection time. The substrate-robustness of the fit does NOT carry through to substrate-robustness of the per-species projection ranking without this filter.

4. Cite Soroye et al. (2020) with cito:confirms for any substrate-robust Iberian or Western Palearctic Bombus replication based on the GLMM fit alone; use cito:qualifies if your contribution is about projection-time substrate-coupling.
```

### Conditions under which the synthesis applies (textarea, required)

Scope: data types, methods, domains, regions, time periods.

```
- Region: Iberian peninsula only.
- Species set: Bombus species observed in GBIF on the Iberian peninsula in the 1901–2014 period (31 species in the combined fit).
- Climate forcing for fit: CRU TS 3.24.01 monthly temperature and precipitation, identical to Soroye et al. 2020 (sourced from their Figshare deposit).
- Spatial substrates tested: equal-area cylindrical (CEA, ~100 km cells), HEALPix-NESTED nside=64 (~92 km), HEALPix-NESTED nside=128 (~46 km) — all on the WGS84 ellipsoid for the HEALPix variants.
- GLMM specification: Soroye et al. 2020's full formula including all main effects, four predictor interaction terms, and per-species random intercept.
- Inference: bambi/PyMC NUTS, 4 chains × 2000 samples, authoritative posterior HDIs reported.
- Synthesis statement applies to the GLMM FIT step, i.e., the substrate-robustness of the headline coefficient on TEI_delta at the historical training period.
```

### Limitations of the synthesis (textarea, required)

What was not tested? What might not generalise?

```
1. Three substrates only (CEA, HEALPix nside=64, HEALPix nside=128). Whether the substrate-robustness extends to coarser HEALPix levels (nside=32, nside=16) or to non-HEALPix grids (EASE-Grid 2.0, S2, etc.) was not directly tested. The within-HEALPix factor-of-2 resolution change is the strongest case where substrate stability is expected; finding stability there is necessary but not sufficient for substrate generality.

2. One region only (Iberian peninsula). The same substrate-robustness pattern may or may not hold for high-latitude Bombus systems, boreal/alpine systems where species' niche margins are different, or other Bombus-rich regions (Andes, Caucasus).

3. One climate dataset (CRU TS 3.24.01 from Soroye's Figshare). Substituting another reanalysis (ERA5, CHELSA, WorldClim) could shift the headline coefficient by an unknown amount; this synthesis does NOT address climate-input robustness.

4. Fit step only. This synthesis concerns the GLMM FIT — the substrate-robustness of the trained model parameters. It does NOT claim substrate-robustness of the future projection step, which is documented separately and is qualified by the substrate-sensitivity sibling chain.

5. The CEA Pass-1 replication uses statsmodels variational-Bayes; the HEALPix passes use both VB and full NUTS. Cross-substrate comparison of HDIs assumes both fitting strategies target the same Bayesian posterior, which holds for these data within the VB-underestimation noise.
```

### Completion date (date picker, required)

```
2026-05-11
```

### Supporting sources (repeatable group, required ≥1)

Each entry is a URL — typically the FORRT Outcome URIs being synthesised. Pull from `nanopubs/PUBLISHED.md` (and/or registries from sibling repos).

- This chain's Outcome (step 05): `<replace-with-published-Outcome-URI-from-step-05>`
- Sibling nside=128 chain's Outcome: `<replace-with-nside128-sibling-Outcome-URI>` (or, before publication: https://doi.org/10.5281/zenodo.20113780)
- Sibling substrate-sensitivity chain's Outcome: `<replace-with-substrate-sensitivity-Outcome-URI>` (or, before publication: https://doi.org/10.5281/zenodo.20113786)
- Soroye, Newbold & Kerr (2020): https://doi.org/10.1126/science.aax8591

### Search topics (Wikidata) (multi-select, optional)

Provide labels (not QIDs).

- bumblebee
- climate change
- species distribution model
- generalized linear mixed model
- HEALPix
- replication
- Iberian Peninsula

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 08.
