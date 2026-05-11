# 07 — Research Software (optional)

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Scope check:** Research Software nanopubs describe **reusable software artefacts** — tools people would `pip install` or `git clone` to use in their own work. They do NOT describe one-off demo / reproduction repos. If your repo is a reproduction of someone else's paper, the reusable artefact is the *upstream library* it uses (e.g. `foscat`, `planktonclas`), not your reproduction repo. Author the Research Software nanopub for the upstream tool, not the demo. See `CLAUDE.md` § Layered architecture: FORRT vs Research Software.

**Form heading:** *"Research Software — Describe research software with metadata including repository, supporting publications, and related resources."*

## Field-by-field draft

### URI of published software (text input, required)

Zenodo concept DOI URL when available, or a GitHub URL. Full URL form.

```
https://doi.org/10.5281/zenodo.20113777
```

### Software Title (text input, required)

The full name or title of the software.

```
weatherxbiodiversity-projection — Iberian Bombus extirpation under DestinE Climate DT SSP3-7.0 (canonical replication of Soroye et al. 2020)
```

### Repository URL (text input, required)

```
https://github.com/annefou/weatherxbiodiversity-projection
```

### Research Project (text input, optional)

URI of the FORRT Claim or PCC question this software is associated with — pull from `nanopubs/PUBLISHED.md`. This is the back-link to the FORRT chain.

```
<replace-with-published-Claim-URI-from-step-03>
```

### License (text input, optional)

```
https://spdx.org/licenses/MIT.html
```

### Related Datasets (repeatable group, optional)

Input data DOIs (Zenodo data records, dataset DOIs, ESA product DOIs).

- Soroye Figshare deposit (CRU TS 3.24.01 + Kerr species list): https://doi.org/10.6084/m9.figshare.10058340
- GBIF Iberian Bombus occurrence download DOI: <replace-with-GBIF-download-DOI-once-issued>
- DestinE Climate DT SSP3-7.0 (IFS-NEMO standard, native HEALPix nside=128) — accessed via polytope on LUMI, not a fixed DOI; cite the DestinE Data Lake catalogue entry

### Related Publications (repeatable group, optional)

One-way back-links to the FORRT Outcome URI(s) the software implements, plus any cited methods papers.

- FORRT Outcome from step 05 of this chain: `<replace-with-published-Outcome-URI-from-step-05>`
- Soroye, Newbold & Kerr (2020) — original paper whose mechanism this software replicates: https://doi.org/10.1126/science.aax8591

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 07.
