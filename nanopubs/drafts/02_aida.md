# 02 — AIDA Sentence

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Form heading:** *"AIDA Sentence — Make structured scientific claims following the AIDA model"*

## Field-by-field draft

### AIDA sentence (textarea, required)

Atomic, Independent, Declarative, Absolute. One empirical finding. Must end with a full stop.

> _If your draft AIDA contains "and" linking two distinct findings, split into two AIDA nanopubs._

```

```

### Select related topics/tags (dropdown, optional)

Predefined topic vocabulary — list the labels you intend to pick from the dropdown.

```

```

### Relates to this nanopublication (text input, required)

URI of the nanopub the AIDA derives from.

- For paper-rooted chains: the Quote-with-comment URI (from step 01).
- For question-rooted chains: the PICO or PCC URI (from step 01).

Pull the URI from `nanopubs/PUBLISHED.md`.

```

```

### Supported by datasets (repeatable group, optional)

DOIs/URLs of datasets that ground the AIDA claim.

- _DOI 1: ___
- _DOI 2: ___

### Supported by other publications (repeatable group, optional)

DOIs/URLs of publications that support the AIDA claim — e.g. peer-reviewed methods papers, or the original paper if not already cited via the Quote.

- _DOI 1: ___
- _DOI 2: ___

> **Known platform bug (2026-04-26):** if both *Supported by datasets* AND *Supported by other publications* are populated and publishing fails, fall back to publishing this AIDA via Nanodash. The URI namespace becomes `https://w3id.org/np/...` (still valid and citable).

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 02.
