# 03 — FORRT Claim

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Form heading:** *"FORRT Claim — Declare an original claim according to FORRT, linking it to an AIDA sentence with a specific FORRT type."*

## Field-by-field draft

### Short URI suffix as claim ID (text input, required)

Slug becomes part of the nanopub URI. Use kebab-case.

```

```

### Label of the claim (text input, required)

A descriptive title (not a sentence). Used for searches/discovery.

```

```

### Search for an AIDA sentence (search/select, required)

URI of the AIDA published in step 02. Pull from `nanopubs/PUBLISHED.md`.

> _If the AIDA was published via Nanodash (`w3id.org/np/...` namespace), the platform's search may not find it — paste the URI manually._

```

```

### Type of FORRT claim (dropdown, required)

Pick one. See `docs/claim-type-vocabulary.md` for the seven options and how to choose.

- [ ] computational performance
- [ ] scalability
- [ ] data quality
- [ ] data governance
- [ ] descriptive pattern
- [ ] model performance
- [ ] statistical significance

### Source URI (text input, optional)

Full URL form: `https://doi.org/...` (NOT bare DOI).

```
https://doi.org/{{PAPER_DOI}}
```

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 03.
