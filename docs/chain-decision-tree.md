# `docs/chain-decision-tree.md` — paper-rooted vs question-rooted chains

A FORRT nanopublication chain on Science Live always has six steps. The first step depends on whether your replication has an **upstream paper** that makes the claim being tested.

## Decision

Ask: *"Is there a sentence in someone else's paper that we are testing or extending?"*

- **Yes** → **paper-rooted chain.** Start with **Quote-with-comment**.
- **No** → **question-rooted chain.** Start with **PICO** (comparative) or **PCC** (descriptive).

For a literal *replication* (tests an existing paper's claim with new data and/or new methods), the answer is almost always yes → paper-rooted.

For first-of-its-kind methodology work where no prior paper makes the equivalent claim → question-rooted.

## Chain shapes

### Paper-rooted (Quote-anchored)

```
Quote-with-comment  →  AIDA  →  FORRT Claim  →  Replication Study  →  Replication Outcome  →  CiTO Citation
```

- **Quote-with-comment** — verbatim sentence from the paper + your comment on why it matters. The DOI of the paper goes in `Cited DOI`.
- **AIDA** — atomic, independent, declarative, absolute version of the quoted claim, stripped of method language and pronouns. The Quote URI goes in `Relates to this nanopublication`.
- **FORRT Claim** — names the claim genre (`docs/claim-type-vocabulary.md`) and links the AIDA URI.
- **Replication Study** — design: scope (what part of the claim) + method (how) + deviations from the original.
- **Replication Outcome** — result: conclusion + evidence + limitations + repo URL.
- **CiTO Citation** — citing work = the Outcome URI; cited work = the original paper DOI; intention = `confirms` / `qualifies` / `disputes` based on the Outcome's validation status.

### Question-rooted (PICO/PCC-anchored)

```
PICO or PCC  →  AIDA  →  FORRT Claim  →  Replication Study  →  Replication Outcome  →  CiTO Citation
```

- **PICO** if the question has a clear comparator (X versus Y). **PCC** if the question is descriptive/scoping.
- **AIDA** — even for question-rooted chains, an AIDA sentence is required because the FORRT Claim form's "Search for an AIDA sentence" field is required. Write the AIDA as the predicted/established declarative answer to the PICO/PCC question.
- The rest of the chain is identical to the paper-rooted shape.

The CiTO Citation in a question-rooted chain still points to *something* — typically the foundational reference paper(s) for the methodology, with intention `usesMethodIn` or `citesAsAuthority`.

## Why two chain shapes

The Science Live FORRT chain is built around a **claim** as the central unit. A claim must come from somewhere:

- A paper makes the claim → cite it as a quote.
- The researcher poses the claim themselves → frame it as an answered research question.

Without one of these anchors, the FORRT Claim has no upstream provenance and the chain is unmoored.

## Common mistake — using PCC for a paper-rooted chain

Don't treat PCC as a generic "research question" template that you can use whenever convenient. PCC encodes a **scoping review** — Population, Concept, Context — for synthesising literature, not for testing a single paper's claim. If you have a paper, the right anchor is the Quote, not a PCC.

The mirror mistake — using Quote-with-comment to anchor a question-rooted chain — happens when there's a *related* paper but no specific sentence we're testing. In that case, cite the related paper at the AIDA step's *Supported by other publications* group, but anchor the chain in PICO/PCC.

## What happens after Phase 5

Once a single chain is published, you have three optional layers:

- **Research Software nanopub** — for any reusable software artefact the work produced (typically an upstream library, not the demo repo). Cites back to the FORRT Claim URI as `Research Project`. See `docs/forrt-form-fields.md` § Research Software.
- **Research Synthesis nanopub** — when this chain is one of several testing facets of a shared underlying property. The Synthesis lists multiple Outcomes as `Supporting sources` and names the cross-cutting conclusion.
- **Wikidata integration** — Outcomes with `cito:confirms` / `cito:disputes` citations are pickup-eligible by the Wikidata Scholia pipeline.

## When in doubt

If you're not sure whether your work is paper-rooted or question-rooted:

- *Can you point to a specific sentence in a paper?* → paper-rooted.
- *Are you the first to publish a claim of this shape?* → question-rooted.
- *Are you between the two?* (e.g. extending one paper using methods from another, with no single sentence being tested) → paper-rooted, with the *primary* paper providing the Quote and the secondary paper(s) appearing in the AIDA's *Supported by other publications* group.

Stop and ask the user before drafting if the answer isn't clear from the paper PDF in `paper/`.
