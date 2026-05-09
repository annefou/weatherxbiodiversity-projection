# 07 — Research Software (optional)

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Scope check:** Research Software nanopubs describe **reusable software artefacts** — tools people would `pip install` or `git clone` to use in their own work. They do NOT describe one-off demo / reproduction repos. If your repo is a reproduction of someone else's paper, the reusable artefact is the *upstream library* it uses (e.g. `foscat`, `planktonclas`), not your reproduction repo. Author the Research Software nanopub for the upstream tool, not the demo. See `CLAUDE.md` § Layered architecture: FORRT vs Research Software.

**Form heading:** *"Research Software — Describe research software with metadata including repository, supporting publications, and related resources."*

## Field-by-field draft

### URI of published software (text input, required)

Zenodo concept DOI URL when available, or a GitHub URL. Full URL form.

```
{{ZENODO_DOI}}
```

### Software Title (text input, required)

The full name or title of the software.

```

```

### Repository URL (text input, required)

```
https://github.com/{{REPO_ORG}}/{{REPO_NAME}}
```

### Research Project (text input, optional)

URI of the FORRT Claim or PCC question this software is associated with — pull from `nanopubs/PUBLISHED.md`. This is the back-link to the FORRT chain.

```

```

### License (text input, optional)

```
https://spdx.org/licenses/MIT.html
```

### Related Datasets (repeatable group, optional)

Input data DOIs (Zenodo data records, dataset DOIs, ESA product DOIs).

- _Dataset URL 1: ___
- _Dataset URL 2: ___

### Related Publications (repeatable group, optional)

One-way back-links to the FORRT Outcome URI(s) the software implements, plus any cited methods papers.

- _Publication URL 1 (FORRT Outcome from step 05): ___
- _Publication URL 2 (methods paper, optional): ___

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 07.
