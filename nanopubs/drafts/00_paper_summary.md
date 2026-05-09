# Paper summary

> This is a working scratchpad for the paper-analysis phase. The output of this file feeds the Quote / AIDA / Claim drafts. It is not itself a nanopub.

**Reference paper:** Climate change contributes to widespread declines among bumble bees across continents

**DOI:** 10.1126/science.aax8591

**Authors:** Peter Soroye, Tim Newbold, Jeremy Kerr

**Year:** 2020

**Venue:** *Science* 367(6478): 685–688, 7 February 2020.

## Headline claim

The single sentence in the paper that this replication tests. Should be one of the paper's core empirical assertions, not a definition or framing statement.

> Verbatim from the paper PDF, p. 685, Abstract:

> "Increasing frequency of hotter temperatures predicts species' local extinction risk, chances of colonizing a new area, and changing species richness."

(149 characters — fits whole-text Quote mode.)

This sentence is the paper's headline empirical assertion: it states the mechanism (frequency of hotter temperatures exceeding historical tolerances) and the three outcomes the analysis predicts (local extinction, colonization, richness change). It is the sentence the replication will test, because the Iberian future-projection study applies the same mechanism (climatic position relative to species' historical thermal niche → occupancy / extinction probability) to a different region, dataset, and time window.

A nearby alternative — *"Increasing frequencies of temperatures that exceed historically observed tolerances help explain widespread bumble bee decline."* (also Abstract, p. 685) — is more conclusion-like but mixes mechanism with the historical decline observation. The chosen sentence is preferred because it cleanly states the mechanism the replication propagates forward in time.

## Methodology summary

- **Data sources:** ~550,000 georeferenced bumble bee occurrence records for 66 *Bombus* species across North America and Europe, harmonised into 100 km × 100 km quadrats. Climate covariates were monthly maximum / minimum temperatures and total precipitation; "historically observed tolerances" were derived from each species' baseline-period occupancy. (See main text p. 685–686 and refs 1, 9.)
- **Time windows:** baseline period 1901–1974 vs recent period 2000–2014 (Fig. 1, Fig. 3). Note: the abstract / Fig. 1 caption use "2000–2015" while the main text and Fig. 3 use "2000–2014" — minor inconsistency in the paper.
- **Statistical model:** detection-corrected occupancy models per species per period; then phylogenetic generalized linear mixed models (GLMMs) of change in probability of site occupancy on a community-averaged "climatic position index" (thermal and precipitation), where position = where local climate sits between each species' historically observed cold/dry and hot/wet limits (0 = cold/dry limit, 1 = hot/wet limit). Predictors: baseline thermal & precipitation position, change-since-baseline in each, and their interaction; controlled for continent. Species richness modelled separately via analysis of covariance on detection-corrected richness.
- **Sample sizes:** 66 species total (35 North American + 36 European, 5 shared species); occupancy quadrats split at first/third quantile into near-cold-limit (n = 969 sites), middle (n = 11,793 sites), near-hot-limit (n = 2,244 sites). Marginal R² = 0.11 for occupancy probability; marginal R² = 0.53 to 0.87 for extirpation/colonization separately.
- **Headline numerical results the replication must compare against:**
  - Mean probability of site occupancy declined by **46% (±3.3% SE) in North America** and **17% (±4.9% SE) in Europe** relative to baseline (p. 686, Results; Fig. 2).
  - Climatic-position model fit: **marginal R² = 0.11** for occupancy; **marginal R² = 0.53 to 0.87** for extirpation/colonization modelled separately (p. 686, Results).
  - Climatic-position predictors outperformed mean / max / min / precipitation predictors by **2.6% lower to 23% higher marginal R²**, ΔDIC = **98.7 to 241.9** (p. 687, Discussion).
  - Pagel's λ = 0.12 (weak phylogenetic signal in occupancy response) (p. 686).

## Replication design choice

- [ ] **Reproduction Study** — direct reproduction: same methodology, same tools.
- [x] **Replication Study** — replication with different methodology or conditions.
- [ ] **Reproduction/Replication Study** — both.

**Justification.** This is a Replication Study, not a Reproduction. We propagate the *mechanism* validated in Soroye et al. 2020 — climatic position relative to a species' historical thermal niche → local extinction / occupancy probability — forward in time to a future-projection setting: Iberian Peninsula bumble bees under SSP3-7.0, driven by the Destination Earth Climate Digital Twin (~5 km, daily) rather than the ~50 km CRU TS reanalysis the original paper used. What is **preserved**: the climatic-position index construction, the linkage from exceedance of historical thermal tolerances to occupancy decline, and the use of long-term occurrence records to derive species-specific limits. What is **different**: (i) data source — DestinE Climate DT instead of CRU TS; (ii) region — Iberian Peninsula instead of continent-scale North America + Europe; (iii) temporal scope — forward projection under a future scenario rather than retrospective baseline-vs-recent contrast; (iv) species pool — Iberian *Bombus* assemblage only; (v) spatial resolution — ~5 km daily climate vs ~50 km monthly. The replication therefore tests transportability of the mechanism, not numerical reproducibility of the original 46% / 17% headline.

## Notes for downstream drafts

- The headline sentence chosen for `01_quote.md` is the **mechanism** sentence from the abstract, not the **conclusion** sentence. The Personal Comment field should make clear the replication propagates this mechanism forward into a future projection on a different continent — i.e., we do not expect to recover the 46% / 17% numbers; we expect the *direction* and *climatic-position dependence* to hold (or not).
- Author list confirmed from p. 685 byline: Peter Soroye¹, Tim Newbold², Jeremy Kerr¹. Affiliations: ¹University of Ottawa, ²University College London. Corresponding author: peter.soroye@gmail.com.
- Data repository for the original paper: figshare (ref 21, doi: 10.6084/m9.figshare.9956471). Not required for this Replication Study, but a Reproduction-Study extension could pull from there.
- Note the time-window inconsistency in the original paper (abstract / Fig. 1: "2000–2015"; Fig. 3 caption and main text: "2000–2014") — flag in `05_outcome.md` if relevant; do not silently fix.
- The original paper's marginal R² of 0.11 for the occupancy model is modest. The replication should not over-interpret a similarly modest fit on Iberian data as a "weak" outcome — it would be consistent with the original paper's effect size.
