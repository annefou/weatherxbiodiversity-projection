# 04 — FORRT Replication Study

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Verify code first:** read the actual reproduction script in `notebooks/03_analysis.py` before writing the methodology field. See `docs/verify-before-drafting.md`.

## Field-by-field draft

### Short URI suffix for study ID (text input, required)

Slug. Use kebab-case.

```
soroye2020-iberian-bombus-cea-and-healpix-nside64
```

### Label/name of replication study (text input, required)

Human-readable title.

```
Soroye et al. 2020 TEI mechanism replicated on Iberian Bombus — CEA grid + HEALPix nside=64 substrate extension
```

### Study type (dropdown, required)

- [ ] Reproduction Study — direct reproduction: same methodology, same tools.
- [x] **Replication Study** — replication with different methodology or conditions.
- [ ] Reproduction/Replication Study — both.

The CEA fit (Tier 1, pass 1) is a direct reproduction of Soroye on a new region (Iberia). The HEALPix nside=64 fit (Tier 1, pass 2) extends the methodology to a different spatial substrate (HEALPix-NESTED on the WGS84 ellipsoid, ~92 km cells) — a methodological replication. The dominant character of this Study is the substrate extension, hence Replication Study.

### Search for a FORRT claim (search/select, required)

URI of the Claim published in step 03. Pull from `nanopubs/PUBLISHED.md`.

```
<replace-with-published-Claim-URI-from-step-03>
```

### Describe what part of the claim is reproduced/replicated (textarea, required)

The **scope** of the claim being tested. Which aspect, what's in/out of scope. NOT methodology. NOT results. See `docs/pico-study-outcome-levels.md`.

```
SCOPE: the GLMM coefficient on standardised TEI_delta (Soroye et al. 2020's "climatic position index change" between baseline 1901–1974 and recent 2000–2014), restricted to Bombus species observed on the Iberian peninsula in GBIF.

IN SCOPE
  - Soroye 2020's GLMM specification: extinction ~ continent + sc_sampling + sc_TEI_bs + sc_TEI_delta + sc_PEI_bs + sc_PEI_delta + (TEI:PEI interactions) + (1|species).
  - Soroye 2020's CRU TS 3.24.01 climate inputs (Figshare bundle, identical NetCDFs).
  - Two spatial substrates: original CEA grid (~100 km cells, Pass 1) and HEALPix-NESTED nside=64 on the WGS84 ellipsoid (~92 km cells, Pass 2).
  - Tier 1 historical fit on the 1901–1974 baseline period and 2000–2014 recent period.

OUT OF SCOPE for this Replication Study (handled by separate chains)
  - Higher-resolution substrates (HEALPix nside=128) — see weatherxbiodiversity-projection-nside128.
  - Cross-substrate methodological diagnostic (substrate-coupling at projection time) — see weatherxbiodiversity-substrate-sensitivity.
  - Tier 2 SSP3-7.0 future projection — reported as a separate Outcome contribution within this same chain (see step 05), not as a separate Study.
  - Bombus species outside the Iberian peninsula (Soroye's original ranges in North America + whole-Europe).
```

### Describe how the claim is reproduced/replicated (textarea, required)

The **method** in plain prose. Read `notebooks/03_analysis.py` and any config files first. NOT exact numerical results.

```
METHOD

Tier 1 — historical fit.

(1) Climate inputs. Soroye's bundled CRU TS 3.24.01 NetCDFs (mean monthly temperature and total monthly precipitation, 1901–2014, 0.5° lat-lon grid) are downloaded from Soroye's Figshare deposit and used unchanged. No substitution with cdsapi or other reanalysis sources at fit time — by design, to keep the climate forcing identical to the original.

(2) Occurrence inputs. GBIF Iberian Bombus occurrences are issued via own-DOI download (cited in the Outcome step 05 "evidence" field). Filtering follows Soroye's species list (genus Bombus, valid taxonomy, occurrences geo-referenced to the Iberian peninsula bounding box).

(3) GLMM specification. Replicates Soroye's exact formula:
    extinction ~ continent + sc_sampling + sc_TEI_bs + sc_TEI_delta + sc_TEI_bs:sc_TEI_delta + sc_PEI_bs + sc_PEI_delta + sc_PEI_bs:sc_PEI_delta + sc_TEI_bs:sc_PEI_bs + sc_TEI_delta:sc_PEI_delta + (1|species)
where TEI = (T_obs − T_min_spp) / (T_max_spp − T_min_spp), PEI is the analogous precipitation index, and sc_ is within-substrate z-scoring.

(4) Two spatial substrates.
    Pass 1 (CEA): cell coverage and per-species niche limits computed on Soroye's original CEA grid (~100 km cells); GLMM fit on this grid.
    Pass 2 (HEALPix nside=64): cell coverage re-computed on HEALPix-NESTED nside=64 cells of the WGS84 ellipsoid (using the healpix-geo Python library); per-species niche limits and the GLMM are refit at this substrate.

(5) Inference. Two independent fitting strategies on the same Pass-2 design matrix:
    (a) variational-Bayes mean-field via statsmodels.BinomialBayesMixedGLM (fast first pass);
    (b) full-posterior NUTS via bambi/PyMC, 4 chains × 2000 samples (authoritative HDIs).

Tier 2 — SSP3-7.0 future projection — runs on the DestinE Jupyter platform and is part of this same chain's Outcome (step 05), not a separate Study. Briefly: DestinE Climate DT IFS-NEMO standard SSP3-7.0 GRIB files are retrieved via polytope on LUMI for the 2020–2029 and 2030–2039 horizons, decoded with eccodes (Python API, NESTED-aware), subset to pre-computed Iberian HEALPix cells, aggregated to monthly means matching the CRU TS units, and used to compute future TEI_delta and PEI_delta per species per cell. Per-species ranking is reported following the protocol established in the methodological sibling chain.

Code: notebooks/01_data_download.py (CRU TS + GBIF), notebooks/02_data_clean.py (CEA Pass 1 cleaning), notebooks/02h_data_clean_healpix.py (HEALPix Pass 2 cleaning), notebooks/03_analysis.py (CEA fit), notebooks/03h_analysis_healpix.py (HEALPix fit), notebooks/04_figures.py + 04h_figures_healpix.py (figures).
```

### Describe any deviations from original methodology (textarea, optional)

What's different from the original method. Verify against the actual code, don't guess.

```
1. Pixelisation. Soroye 2020 fit at the CEA grid (~100 km, equal-area cylindrical projection at the equator). This Study additionally fits at HEALPix-NESTED nside=64 (~92 km cells on the WGS84 ellipsoid). Cell area is comparable (within ~10%); cell shape and topology differ.

2. Region. Soroye 2020 fit on North America + whole Europe (all Western Palearctic Bombus). This Study fits on Iberian peninsula only. Iberia is a subset of Soroye's European fitting domain; the same GLMM specification applies.

3. Inference engine. Soroye 2020 used MCMCglmm (R). This Study uses two independent Python implementations (statsmodels variational-Bayes for the fast first pass; bambi/PyMC NUTS for the authoritative posterior HDIs reported in the Outcome). Both target the same Bayesian posterior.

4. Prior choices. Soroye 2020's MCMCglmm priors (default informative inverse-Wishart on variance components) are not exactly reproducible across to bambi (Half-StudentT default on group SDs). This is the largest non-trivial methodological deviation; it should not affect the headline coefficient sign or order of magnitude (verified: it does not). The full prior table is in notebooks/03h_analysis_healpix.py.

5. Climate inputs. Soroye 2020's CRU TS 3.24.01 from his Figshare deposit is used unchanged — NOT a deviation, this is the SAME climate dataset.

6. Sampling-effort proxy. Same as Soroye (cells_per_decade × records_per_cell from the cleaned occurrence table) — NOT a deviation.
```

### Search keywords (Wikidata) (multi-select, optional)

Provide labels (not QIDs) — the Wikidata search picks up labels.

- bumblebee
- climate change
- generalized linear mixed model
- HEALPix
- species distribution model
- Iberian Peninsula
- replication study

### Search discipline (Wikidata) (search, optional)

Provide labels.

- macroecology
- biogeography
- conservation biology

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 04.
