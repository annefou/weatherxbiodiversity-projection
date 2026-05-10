# Substrate-sensitivity findings — TEI projection at HEALPix nside=64 vs nside=128

**Inputs:** two single-substrate replications of Soroye et al. (2020) for Iberian *Bombus*:

- `weatherxbiodiversity-projection` — HEALPix nside=64 (~92 km cells), full GLMM refit
- `weatherxbiodiversity-projection-nside128` — HEALPix nside=128 (~46 km cells), full GLMM refit

**Climate inputs.** Tier-1 (historical fit + niche limits) uses CRU TS 3.24.01 from Soroye's own Figshare deposit — the exact climate data Soroye fit the original GLMM on, kept unchanged so the replication varies only region and pixelisation. Tier-2 (future projection) uses DestinE Climate DT SSP3-7.0 (IFS-NEMO standard, native HEALPix nside=128) via polytope on LUMI — the only place new climate data enters, and only for the future horizon.

**Occurrences.** GBIF download DOIs (own-issued, Iberia subset of *Bombus*).

**Question:** when projecting per-species extirpation risk to a future climate, how stable is the per-species ranking under a change of pixelisation grid?

## Headline result

Per-species ranking under SSP3-7.0 is **substrate-sensitive** for the full GLMM η reported in either single-substrate study, but the sensitivity is **predictable**, **mechanistic**, and **fixable** by either (i) restricting the report to species observed in ≥10 cells per substrate AND dropping the GLMM interaction terms at projection time, or (ii) reporting a substrate-invariant physical metric (mean future TEI, or fraction of currently-occupied cells with future TEI > 0.5).

The substrate-coupling is **not** caused by per-species random-intercept refit (which is substrate-stable, ΔRE ≤ 0.6 logits across all species), nor by per-species niche-limit shifts (which are modest, ΔT_range = 0–3°C). It is caused by:

1. **Per-species sample size at projection time.** Species observed in fewer than ~5 historical cells per substrate produce per-cell-mean η that is dominated by extrapolation noise from very few cells. Above ~10 cells, the ranking signal is robust.
2. **The GLMM interaction terms (sc_TEI_delta:sc_PEI_delta in particular) compound substrate-specific standardisation quadratically.** Future-period climate deltas extrapolate 2–4σ outside the training distribution; the same physical signal z-scores to opposite tails of the distribution under different substrates' (μ, σ); the interaction term scales as the product, amplifying the divergence. This is the largest single contributor to projection η at SSP3-7.0 horizons (typically +1.4 to +7.7 logits per species), so its substrate-sensitivity dominates the total η ranking.

## Five projection variants tested

For each substrate, four candidate per-species rankings were computed, plus a fifth substrate-invariant metric:

| Variant | Predictor at projection | Standardisation | Notes |
|---|---|---|---|
| (a) Full GLMM η | full | within-substrate (μ_64 / μ_128) | the headline statistic in each repo's Outcome |
| (b) Main-effects only η | drop interactions | within-substrate | tests whether extrapolation amplification via interactions is the driver |
| (c) Full η, shared CEA reference | full | shared CEA pool (μ_CEA, σ_CEA) | tests whether the standardisation reference alone resolves divergence |
| (d) Mean future TEI | physical mean | none | substrate-invariant in physical units; the natural Soroye mechanism |
| (d2) Frac TEI_future > 0.5 | threshold count | none | substrate-invariant; mid-range thermal exceedance |

## Cross-substrate concordance (Spearman ρ)

Horizon **2030–2039**, n = 31 species in both substrates:

| Variant | ρ at n≥1 | ρ at n≥5 | ρ at n≥10 | ρ at n≥20 |
|---|---:|---:|---:|---:|
| (a) Full GLMM η | +0.27 | +0.51 | +0.59 | +0.77 |
| (b) Main-effects only η | +0.40 | +0.52 | **+0.97** | **+0.98** |
| (c) Shared CEA ref η | +0.27 | +0.49 | +0.52 | +0.55 |
| (d) Mean future TEI | +0.66 | +0.69 | +0.90 | +0.82 |
| (d2) Frac TEI_future > 0.5 | +0.66 | +0.71 | +0.88 | +0.83 |

Horizon **2020–2029** shows the same qualitative pattern: (b) at n≥10 hits ρ = +0.96; (d) at n≥10 hits +0.90; (a) and (c) lag.

## Mechanism — confirmed and refuted hypotheses

### Confirmed

- **Sample size dominates.** All variants improve monotonically as the n_cells filter tightens; even (a) reaches ρ = +0.77 at n≥20.
- **GLMM interaction terms amplify substrate divergence at extrapolation.** Dropping them (variant b) lifts ρ from +0.59 to +0.97 at n≥10. The mechanism is quadratic compounding of substrate-specific z-scores: when future TEI_delta lies 3σ outside training distribution, an σ_64 ≠ σ_128 mismatch becomes a ~3× larger gap in the product term.
- **Physical-unit metrics decouple substrates well.** Mean future TEI agrees at ρ ≥ 0.66 at all sample sizes — better than any GLMM-based variant at n < 10.

### Refuted

- **Per-species random-intercept shrinkage is NOT the driver.** ΔRE across substrates is small (≤ 0.6 logits) for every species in the comparison set. The GLMM hyperparameters are substrate-stable; the substrate-coupling is not in the model parameters.
- **Per-species niche limits (T_min_spp, T_max_spp) are NOT the driver.** Range shifts are 0–3°C; even species with identical T_range across substrates (norvegicus: ΔT_range = +0.01°C) show 8.8-logit η swings. Niche refit is not enough to explain divergence.
- **Re-standardising future predictors against a shared reference (variant c) does NOT fix it.** ρ stays near (a) at all n filters, because the GLMM β coefficients were fit against substrate-local σ — just rescaling the predictors at projection time without refitting β does not decouple the model. To properly test this fix, the GLMM would need to be refit with shared standardisation (out of scope here).

## Recommended reporting protocol

For any future TEI-based extirpation projection that is to be compared across spatial substrates:

1. **Report only on species with ≥ 10 currently-occupied + active cells per substrate.** Below this, per-species mean η is dominated by per-cell extrapolation noise, regardless of variant choice.
2. **Drop the GLMM interaction terms at projection time** (keep them in the fit; they are part of the original Soroye specification, but they should not be used to extrapolate). Report η using main effects only: Intercept + sc_sampling + sc_TEI_bs + sc_TEI_delta + sc_PEI_bs + sc_PEI_delta + species_RE.
3. **Cross-check against a substrate-invariant physical metric** — mean future TEI per species, or fraction of cells where future TEI > 0.5. If the GLMM main-effects ranking and the physical ranking disagree for a given species, the GLMM is unreliable for that species at that substrate.

## Scientific findings (beyond the substrate-sensitivity diagnostic)

### F1 — Soroye's TEI-based extirpation mechanism replicates on Iberian *Bombus*

Using Soroye's own CRU TS 3.24.01 climate inputs (Figshare deposit) and own-issued GBIF download DOIs for Iberia, the GLMM coefficient on `sc_TEI_delta` is positive, large, and credible at all three pixelisations:

| Substrate | sc_TEI_delta β | 95% HDI |
|---|---:|---|
| CEA (~100 km) | +0.479 | excludes 0 |
| HEALPix nside=64 (~92 km) | +0.454 | [+0.130, +0.751] |
| HEALPix nside=128 (~46 km) | +0.347 | [+0.139, +0.533] |

All three within ±30% of CEA. **Soroye's central biological claim — that thermal-niche exceedance increases extirpation probability — replicates in a new region (Iberia) at three independent pixelisations, holding climate inputs and model specification identical to the original.** This is the headline positive replication.

### F2 — At 2030–2039 the SSP3-7.0 signal is *drift toward* the niche edge, not exceedance

For every Iberian *Bombus* species in this study, the fraction of currently-occupied cells where future TEI exceeds 1.0 (Soroye's actual extirpation threshold) is **0.00 at 2030–2039**. None of these species hits its upper thermal niche at any currently-occupied cell over the next 15 years under SSP3-7.0. What's happening instead is a systematic shift in mean future TEI (~0.43 historically → ~0.45–0.50 by 2039). **The "extirpation event" Soroye projects is a longer-timescale phenomenon — 2030–2039 is the early-warning period, not the extirpation period itself.** The honest framing is "X species are drifting fastest toward their upper niche edge," not "X species will be extirpated by 2039."

### F3 — Substrate-stable high-risk and low-risk Iberian *Bombus*

Filtering to species with ≥10 occupied cells per substrate AND using main-effects-only η (the trustworthy projection variant per the substrate-sensitivity analysis), the top-3 highest-risk and bottom-3 lowest-risk species **agree perfectly between nside=64 and nside=128**:

| Rank | Species | n@64 | n@128 |
|---|---|---:|---:|
| Highest risk | *B. humilis* | 23 | 43 |
| 2 | *B. muscorum* | 25 | 43 |
| 3 | *B. ruderarius* | 16 | 29 |
| ... | | | |
| Safest | *B. terrestris* | 57 | 179 |
| 2nd safest | *B. pascuorum* | 46 | 111 |
| 3rd safest | *B. pratorum* | 25 | 40 |

This is biologically coherent: *B. terrestris* is the dominant Western Palearctic generalist with the broadest realised niche; *B. humilis*, *muscorum*, *ruderarius* are short-tongued grassland specialists already in decline elsewhere in Europe (UK BBS, Belgium, Netherlands). **The substrate-stable Iberian projection independently identifies the same European high-risk species the long-term monitoring data flag.**

### F4 — The species we cannot rank reliably are the ones we most worry about

The 18 species with n_cells < 10 per substrate (which include Pyrenean specialists *B. pyrenaeus*, *B. mucidus*, *B. mendax*, *B. wurflenii*, *B. monticola*, etc.) are precisely the narrowly-distributed alpine endemics. Their per-species ranking under SSP3-7.0 swings 5–9 logits between nside=64 and nside=128. **The TEI projection cannot give a substrate-stable rank for the species the public most worries about.** This is the most uncomfortable finding.

## Implication for the two single-substrate Outcomes

The **Replicated** verdict at CEA (`weatherxbiodiversity-projection` Tier 1) and the **Substrate-robust** verdict at HEALPix nside=64 / nside=128 (Tier 1 of both repos) are unaffected — the GLMM coefficient ratios (sc_TEI_delta/sc_PEI_delta β) agree across substrates within ±30%. The substrate-coupling shows up only at **projection time**, not at fit time.

The Tier-2 per-species ranking under SSP3-7.0 is the substrate-coupled artefact. Both repos' Outcomes should therefore:

- Lead with the substrate-robust GLMM headline (Tier 1 sc_TEI_delta) as the validated result
- Report the per-species projection ranking with the **n_cells ≥ 10 filter and main-effects-only η** as the trustworthy version
- Add a Limitation pointing to this substrate-sensitivity analysis as the basis for the recommendations

## Artefacts generated

- `scripts/compare_substrates.py` — per-species per-term η decomposition (used to identify the interaction-term mechanism).
- `scripts/compare_variants.py` — five-variant cross-substrate concordance, with n_cells filtering.
- `scripts/plot_variant_concordance.py` — heatmap + scatter pair plots.
- `results/substrate_comparison_decomposition_2030_2039.txt` — 12-species decomposition output.
- `results/variant_comparison_2020_2029.{csv,json}` — near-term horizon.
- `results/variant_comparison_2030_2039.{csv,json}` — mid-term horizon.
- `figures/variant_concordance_{2020_2029,2030_2039}.png` — Spearman ρ heatmaps.
- `figures/variant_pairs_{2020_2029,2030_2039}.png` — per-variant scatter, coloured by min(n_cells_64, n_cells_128).
