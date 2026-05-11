# 02 — AIDA Sentence

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Form heading:** *"AIDA Sentence — Make structured scientific claims following the AIDA model"*

## Field-by-field draft

### AIDA sentence (textarea, required)

Atomic, Independent, Declarative, Absolute. One empirical finding. Must end with a full stop.

> _If your draft AIDA contains "and" linking two distinct findings, split into two AIDA nanopubs._

```
On Iberian Bombus, the GLMM coefficient on standardised TEI_delta — the change in the climatic position index (Soroye et al. 2020) between baseline 1901–1974 and recent 2000–2014 climate — is positive and credibly greater than zero at HEALPix-NESTED nside=64 (cell area approximately 92 km), with posterior mean +0.454 and 95% highest-density interval [+0.130, +0.751].
```

### Select related topics/tags (dropdown, optional)

Predefined topic vocabulary — list the labels you intend to pick from the dropdown.

```
biodiversity
climate change
species distribution modelling
replication
```

### Relates to this nanopublication (text input, required)

URI of the nanopub the AIDA derives from.

- For paper-rooted chains: the Quote-with-comment URI (from step 01).
- For question-rooted chains: the PICO or PCC URI (from step 01).

Pull the URI from `nanopubs/PUBLISHED.md`.

```
<replace-with-published-Quote-URI-from-step-01>
```

### Supported by datasets (repeatable group, optional)

DOIs/URLs of datasets that ground the AIDA claim.

- CRU TS 3.24.01 climate (Soroye Figshare bundle, kept identical to original)
- GBIF Iberian Bombus occurrence download DOI (own-issued; record the DOI here once minted)
- This repo's Zenodo concept DOI: `10.5281/zenodo.20113777`

### Supported by other publications (repeatable group, optional)

DOIs/URLs of publications that support the AIDA claim — e.g. peer-reviewed methods papers, or the original paper if not already cited via the Quote.

- 10.1126/science.aax8591 (Soroye, Newbold & Kerr 2020 — original paper; cited via the Quote in step 01)

> **Known platform bug (2026-04-26):** if both *Supported by datasets* AND *Supported by other publications* are populated and publishing fails, fall back to publishing this AIDA via Nanodash. The URI namespace becomes `https://w3id.org/np/...` (still valid and citable).

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 02.
