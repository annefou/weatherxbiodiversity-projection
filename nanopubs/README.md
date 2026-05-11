# `nanopubs/` — FORRT nanopublication chain workspace

This directory holds the field-by-field drafts of the FORRT chain, plus the registry of published URIs. The chain is published manually on `https://platform.sciencelive4all.org` — Claude does not publish; Claude drafts each step in `drafts/`, you review and copy-paste into the platform UI, and the resulting URI goes into `PUBLISHED.md`.

## The chain

A complete paper-rooted FORRT chain is six steps published in order:

1. `drafts/01_quote.md` → Quote-with-comment
2. `drafts/02_aida.md` → AIDA Sentence
3. `drafts/03_claim.md` → FORRT Claim
4. `drafts/04_study.md` → FORRT Replication Study
5. `drafts/05_outcome.md` → FORRT Replication Outcome
6. `drafts/06_citation.md` → CiTO Citation

For question-rooted chains, replace step 1 with `drafts/01_pico.md` or `drafts/01_pcc.md`. See `docs/chain-decision-tree.md` for which to choose.

## Drafting workflow

For each step:

1. Run the pre-flight checklist in `docs/forrt-form-fields.md` — before drafting any content.
2. Open the matching draft file in `drafts/`.
3. Write each form field verbatim, in form order. Required fields: provide a value. Optional fields: provide a value, or write `*(skip — optional)*`.
4. Have the user review.
5. The user copy-pastes into the platform UI on `https://platform.sciencelive4all.org` and publishes.
6. The user records the resulting URI in `PUBLISHED.md`.
7. The next step's draft references the just-published URI (e.g. the AIDA's *Relates to this nanopublication* field is the Quote URI).

## Order matters

Each step in the chain references the URI of the previous step:

- AIDA *Relates to* → Quote URI (or PICO/PCC URI)
- FORRT Claim *Search for an AIDA* → AIDA URI
- Replication Study *Search for a FORRT claim* → Claim URI
- Replication Outcome *Search for a FORRT replication study* → Study URI
- CiTO Citation *citing creative work* → Outcome URI

Don't try to publish out of order. Don't forget to update `PUBLISHED.md` as you go — downstream drafts pull from there.

## Optional layer

Once the six-step chain is published, one optional further nanopub applies for this chain:

- **Research Software** — describes this repo as a reusable research-software artefact (Zenodo-archived, citable, MIT-licensed). Cites back to the FORRT Claim URI as `Research Project`. Drafted in `drafts/07_research_software.md`. See `docs/forrt-form-fields.md` § Research Software and `CLAUDE.md` § Layered architecture: FORRT vs Research Software.

A **Research Synthesis** nanopub is not published from this chain — the natural place for the cross-chain synthesis is the methodological sibling repository `weatherxbiodiversity-substrate-sensitivity`, which combines this Outcome with the nside=128 sibling Outcome and its own diagnostic Outcome into a single Research Synthesis. See that repo's `drafts/08_synthesis.md`.

## What lives in `drafts/` after publishing

After a draft has been published, it stays on disk as the historical record of what was submitted. The `PUBLISHED.md` registry is the canonical source for URIs; the drafts are the canonical source for content (including any optional fields skipped).

If a republish is needed (e.g. the platform rejected the original), edit the draft and republish. Don't delete drafts.
