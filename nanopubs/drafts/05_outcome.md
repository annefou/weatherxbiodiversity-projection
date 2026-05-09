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
https://github.com/{{REPO_ORG}}/{{REPO_NAME}}
```

### Completion date (date picker, required)

```
{{RELEASE_DATE}}
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

```

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 05.
