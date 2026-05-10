# 06 — CiTO Citation

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Description:** *"Declare citations between papers or other works, using Citation Typing Ontology"*

## Field-by-field draft

### Identifier for the citing creative work (text input, required)

URI of the Outcome published in step 05. Pull from `nanopubs/PUBLISHED.md`.

```
<replace-with-published-Outcome-URI-from-step-05>
```

### List citations (repeatable group, required ≥1)

#### Citation 1 — back to the original paper (Soroye 2020)

##### Citation Type (dropdown)

- [x] **`confirms`**

(Outcome verdict is Validated, which maps to CiTO `confirms` per `docs/forrt-form-fields.md`.)

##### DOI or other URL of the cited work (text input)

```
https://doi.org/10.1126/science.aax8591
```

#### Citation 2 — extends the substrate-extension sibling chain

##### Citation Type (dropdown)

- [x] **`extends`**

(The nside=128 sibling chain is a methodological substrate extension that confirms the same TEI mechanism at finer resolution; this Outcome's substrate-robustness conclusion is supported by both substrates jointly. Cite the nside=128 sibling's Outcome URI once published.)

##### DOI or other URL of the cited work (text input)

```
<replace-with-nside128-sibling-Outcome-URI-from-its-PUBLISHED.md>
```

(Or: cite the sibling repo's Zenodo concept DOI: `<replace-with-nside128-Zenodo-concept-DOI>`.)

#### Citation 3 — extends the methodological substrate-sensitivity sibling chain

##### Citation Type (dropdown)

- [x] **`extends`**

(The substrate-sensitivity sibling chain documents the projection-time grid-coupling diagnostic and the recommended reporting protocol that this Outcome's per-species ranking is filtered against.)

##### DOI or other URL of the cited work (text input)

```
<replace-with-substrate-sensitivity-Outcome-URI-from-its-PUBLISHED.md>
```

(Or: cite that repo's Zenodo concept DOI: `<replace-with-substrate-sensitivity-Zenodo-concept-DOI>`.)

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 06.

This completes the six-step FORRT chain. Optional next layers:

- **Research Software** (`drafts/07_research_software.md`) — if the repo *produces* a reusable software artefact.
- **Research Synthesis** (`drafts/08_synthesis.md`) — if this chain is one of several testing facets of a shared property.
