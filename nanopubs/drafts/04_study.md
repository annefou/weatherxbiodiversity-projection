# 04 — FORRT Replication Study

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Verify code first:** read the actual reproduction script in `notebooks/03_analysis.py` before writing the methodology field. See `docs/verify-before-drafting.md`.

## Field-by-field draft

### Short URI suffix for study ID (text input, required)

Slug. Use kebab-case.

```

```

### Label/name of replication study (text input, required)

Human-readable title.

```

```

### Study type (dropdown, required)

- [ ] Reproduction Study — direct reproduction: same methodology, same tools.
- [ ] Replication Study — replication with different methodology or conditions.
- [ ] Reproduction/Replication Study — both.

### Search for a FORRT claim (search/select, required)

URI of the Claim published in step 03. Pull from `nanopubs/PUBLISHED.md`.

```

```

### Describe what part of the claim is reproduced/replicated (textarea, required)

The **scope** of the claim being tested. Which aspect, what's in/out of scope. NOT methodology. NOT results. See `docs/pico-study-outcome-levels.md`.

```

```

### Describe how the claim is reproduced/replicated (textarea, required)

The **method** in plain prose. Read `notebooks/03_analysis.py` and any config files first. NOT exact numerical results.

```

```

### Describe any deviations from original methodology (textarea, optional)

What's different from the original method. Verify against the actual code, don't guess.

```

```

### Search keywords (Wikidata) (multi-select, optional)

Provide labels (not QIDs) — the Wikidata search picks up labels.

- _Label 1: ___
- _Label 2: ___

### Search discipline (Wikidata) (search, optional)

Provide labels.

- _Discipline label: ___

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 04.
