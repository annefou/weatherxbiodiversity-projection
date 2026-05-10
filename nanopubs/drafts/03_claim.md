# 03 — FORRT Claim

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Form heading:** *"FORRT Claim — Declare an original claim according to FORRT, linking it to an AIDA sentence with a specific FORRT type."*

## Field-by-field draft

### Short URI suffix as claim ID (text input, required)

Slug becomes part of the nanopub URI. Use kebab-case.

```
soroye2020-tei-delta-positive-iberia-nside64
```

### Label of the claim (text input, required)

A descriptive title (not a sentence). Used for searches/discovery.

```
TEI_delta GLMM coefficient is positive and credibly above zero on Iberian Bombus at HEALPix nside=64
```

### Search for an AIDA sentence (search/select, required)

URI of the AIDA published in step 02. Pull from `nanopubs/PUBLISHED.md`.

> _If the AIDA was published via Nanodash (`w3id.org/np/...` namespace), the platform's search may not find it — paste the URI manually._

```
<replace-with-published-AIDA-URI-from-step-02>
```

### Type of FORRT claim (dropdown, required)

Pick one. See `docs/claim-type-vocabulary.md` for the seven options and how to choose.

- [ ] computational performance
- [ ] scalability
- [ ] data quality
- [ ] data governance
- [ ] descriptive pattern
- [ ] model performance
- [x] **statistical significance**

The claim is about a posterior credible interval on a GLMM coefficient: the 95% HDI on β_sc_TEI_delta excludes zero with the same sign and order of magnitude as Soroye 2020. This is canonically a "statistical significance" FORRT claim (it asserts a parameter is reliably non-zero in a Bayesian sense). It is not "model performance" because we are not asserting predictive skill, and it is not "descriptive pattern" because the claim is about a fitted parameter, not raw data.

### Source URI (text input, optional)

Full URL form: `https://doi.org/...` (NOT bare DOI).

```
https://doi.org/10.1126/science.aax8591
```

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 03.
