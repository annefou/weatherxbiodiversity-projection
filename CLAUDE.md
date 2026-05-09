# CLAUDE.md — operating manual

@DOMAIN.md
@USER_PREFERENCES.md

This file is the universal operating manual for replication studies built from `forrt-replication-template`. It is **never edited per project**. Domain-specific conventions live in `DOMAIN.md`; per-user style lives in `USER_PREFERENCES.md`. Both are imported above.

## Purpose of this repository

This repository is a self-contained replication study. By the time it is finished, it produces:

1. A reproducible computational pipeline (notebooks + Snakefile + Docker).
2. A Zenodo-archived release with a citable concept DOI.
3. A FORRT-tagged nanopublication chain on the Science Live platform that records the claim being tested, the replication design, the outcome, and the citation relationship to the original paper.

Your job, as Claude, is to help the user move from "freshly cloned template + a paper PDF" through all four phases without needing context outside this directory.

## Standards alignment — FAIR4RS

This template is designed to satisfy the **FAIR Principles for Research Software (FAIR4RS)** (Chue Hong et al. 2022, [RDA00068](https://doi.org/10.15497/RDA00068)). FAIR4RS adapts FAIR (Findable, Accessible, Interoperable, Reusable) to software's specific lifecycle: identifiers, metadata, qualified references, license, provenance, and community standards.

Every artefact in this repo has a FAIR4RS purpose. The mapping is in `docs/fair4rs-checklist.md` — read it before declaring a phase "done", because FAIR4RS conformance is part of what makes a replication credible, not an optional polish step. When in doubt about a release, check the checklist.

Key FAIR4RS handles in this template:

- **F (Findable)** — Zenodo concept DOI minted at first release; `CITATION.cff` + `codemeta.json` populated; GitHub topics suggested in `README.md`.
- **A (Accessible)** — public GitHub repo, MIT license, GHCR image, Zenodo archive of the source tarball and (optionally) the Docker image, so the software remains accessible even if GitHub or GHCR disappears.
- **I (Interoperable)** — `environment.yml` is the single declarative source of truth; standard formats (jupytext `.py` notebooks, parquet/CSV/NetCDF data, OCI Docker image); cross-references to upstream paper DOI, dataset DOIs, and FORRT nanopub URIs use qualified relations.
- **R (Reusable)** — explicit MIT license, detailed provenance via the FORRT chain (every claim traceable to data + method + author ORCID), `CITATION.cff` lists upstream paper as a `references` entry, `codemeta.json` lists `softwareRequirements` and `referencePublication`.

The FORRT nanopublication chain itself is what makes **R1.2 (provenance)** machine-actionable: it does what `CITATION.cff` cannot, by separating *claims about the world* from *facts about the software*, with cryptographically signed nanopubs at each step.

## First-run guard — DO NOT SKIP

Before doing any other work in this repository, run this check:

```bash
grep -r --include='*.md' --include='*.yml' --include='*.json' --include='*.yaml' --include='*.cff' --include='*.toml' '{{[A-Z_]\+}}' . 2>/dev/null | grep -v '^./.claude/' | grep -v '^./CLAUDE.md' | head
```

If the output contains any unsubstituted `{{...}}` token, the template has not been initialised. **Stop**, tell the user:

> "This repository still has unsubstituted placeholder tokens from the template. Run `/init-template` to bootstrap it (you'll be asked for author identity, paper DOI, etc.), then we can proceed."

…and offer to invoke the `init-template` skill (`.claude/skills/init-template/SKILL.md`).

## Workflow phases

A complete replication moves through six phases. Each phase has artefacts and exit criteria; do not promote a phase as "done" without the artefacts.

### Phase 0 — Bootstrap

- `init-template` skill has been run; no `{{...}}` tokens remain.
- `paper/` contains the source paper PDF.
- `USER_PREFERENCES.md` reflects this user's identity and style.
- A `git remote` exists pointing at the GitHub repo created from the template.

Exit: `git status` is clean and the first commit is on `main`.

### Phase 1 — Paper analysis

- The paper has been read end-to-end (use `Read` on the PDF, page by page if large).
- The **headline claim sentence** has been quoted verbatim from the paper into `nanopubs/drafts/01_quote.md`. This is the sentence the replication will test or extend. Verbatim means character-for-character; paraphrase is forbidden (`docs/verify-before-drafting.md`).
- A short **methodology summary** has been written into `nanopubs/drafts/00_paper_summary.md`: data sources, statistical model, sample sizes, headline numerical result.
- The replication's **design choice** has been recorded: are we doing a Reproduction Study (same data, same methods), a Replication Study (different data and/or different methods), or both?

Hand the user the `paper-analyst` agent (`.claude/agents/paper-analyst.md`) when starting this phase.

Exit: a complete `nanopubs/drafts/01_quote.md` with verified verbatim quote and DOI.

### Phase 2 — Code & data port

- `environment.yml` lists every dependency the notebooks import (`docs/cicd-conventions.md` § environment.yml is the single source of truth).
- `notebooks/01_data_download.py` fetches all input data — no manual steps.
- `notebooks/02_data_clean.py` produces a tidy intermediate dataset.
- `notebooks/03_analysis.py` reproduces the paper's headline statistic.
- `notebooks/04_figures.py` produces the figures used in the Jupyter Book.

Anti-patterns to avoid in this phase:

- **Don't claim notebook code "works" without running it** (`docs/verify-before-drafting.md` § Test before claiming ready). If a notebook has not been executed, say so explicitly.
- **Don't extrapolate framework choices** for the Replication Study draft — read the actual code that ran the analysis (`docs/verify-before-drafting.md` § Verify code before drafting).
- **Don't embed `matplotlib.use('Agg')`** in jupytext notebooks. It blocks inline plot display, which breaks the Jupyter Book build. Use `plt.show()` after `fig.savefig()`.

Exit: `snakemake --cores 1` runs end-to-end on a fresh checkout.

### Phase 3 — Local results & figures

- The headline number from the original paper has been compared against the replication's number.
- One main figure (`figures/main_result.png`) shows the comparison.
- The result has been honestly characterised: validated, partially supported, or contradicted (`docs/claim-type-vocabulary.md`). Negative results are publishable; overclaiming is not (`DOMAIN.md` § Scientific validity standard).

Exit: `nanopubs/drafts/05_outcome.md` is written with the conclusion sentence, the supporting numerical evidence, and the limitations.

### Phase 4 — Release & archive

- `CITATION.cff` and `codemeta.json` have correct version + date.
- A GitHub release is cut with a Zenodo-facing description (no internal ops detail, no bot signatures — `docs/cicd-conventions.md` § Release notes are Zenodo descriptions).
- Zenodo mints a concept DOI; the value is written into `CITATION.cff` and `codemeta.json`.
- The Docker image is pushed to GHCR via `.github/workflows/docker.yml` and (optionally) archived on Zenodo.

Exit: the release page is live, the Zenodo record exists, and `nanopubs/PUBLISHED.md` lists the source + image DOIs.

### Phase 5 — FORRT nanopublication chain

A FORRT chain has six steps in order. Each step is published as a separate nanopublication on `https://platform.sciencelive4all.org`. The chain shape depends on whether there is an upstream paper (`docs/chain-decision-tree.md`):

- **Paper-rooted** (almost always for a replication): Quote-with-comment → AIDA → FORRT Claim → Replication Study → Replication Outcome → CiTO Citation.
- **Question-rooted** (no upstream paper): PCC or PICO question → AIDA → FORRT Claim → Replication Study → Replication Outcome → CiTO Citation.

For each step, draft field-by-field into `nanopubs/drafts/0X_<step>.md` using the exact form structure documented in `docs/forrt-form-fields.md`. **Never invent field names**. **Never ship a draft that contains only the headline content without the full field enumeration.** The pre-flight checklist in `docs/forrt-form-fields.md` is non-negotiable; run it before every draft.

Once each draft is approved, the user copies the fields into the Science Live UI and publishes. The published URI goes into `nanopubs/PUBLISHED.md` and is referenced by downstream steps in the chain.

Exit: all six URIs are listed in `nanopubs/PUBLISHED.md`, the Jupyter Book embeds at least one of them, and the chain is browsable from `index.md`.

## Universal anti-patterns (apply across all phases)

These rules apply regardless of domain. Domain-specific anti-patterns live in `DOMAIN.md`.

### Verify quotes verbatim before drafting

The Quote-with-comment template requires **verbatim** text from the source. Generating an "approximate" quote from training-data memory is forbidden, even when the paraphrase is faithful. Always read the PDF in `paper/` first and copy the sentence character-for-character. See `docs/verify-before-drafting.md`.

### Verify code before drafting Study/Outcome

Before drafting the Methodology or Deviations fields of a FORRT Replication Study or Outcome, read the actual reproduction script in `notebooks/`. Don't infer the framework, library, or hyperparameters from prior memory or from what feels "natural" given the paper. See `docs/verify-before-drafting.md`.

### Test before claiming ready

Don't describe untested code as "works" or "ready". Run it first, or label it explicitly: "written, untested" or "written, runs in the Docker container per the API docs, not yet executed." See `docs/verify-before-drafting.md`.

### Atomic AIDA — one empirical finding per AIDA sentence

AIDA = Atomic, Independent, Declarative, Absolute. One AIDA sentence states one empirical finding. If your draft AIDA contains "and" linking two distinct findings, split it into two AIDA nanopubs anchored on two separate Claims. See `docs/forrt-form-fields.md`.

### Don't conflate PICO with method

PICO is the **question**. The Replication Study form's "what is reproduced" is the **scope**. The Replication Study form's "how is reproduced" is the **method**. Results live only in the Outcome. Stuffing methodology into PICO, or numerical results into the Study, makes the chain non-composable. See `docs/pico-study-outcome-levels.md`.

### Don't poll long-running experiments

If an analysis takes more than ~5 minutes, launch it as a background process (`nohup … &`) and tell the user the estimated completion time. Don't burn conversation context polling a results file every few seconds.

### Layered architecture: FORRT vs Research Software

The FORRT chain (Quote → AIDA → Claim → Study → Outcome → CiTO) makes claims about the world. A Research Software nanopub describes a reusable software artefact (a `pip install`-able tool, not a one-off demo repo) and cites the FORRT Claim one-way. Don't pull tooling URIs into the FORRT chain; publish a separate Research Software nanopub if the repo *produces* a reusable tool. See `docs/forrt-form-fields.md` § Research Software.

### Atomic publishing order

Publish nanopubs in the order they cite each other: Quote → AIDA → Claim → Study → Outcome → CiTO. Each step needs the URI of the previous step as input. Don't try to publish out of order.

## Document map

The documents under `docs/` are the load-bearing reference material; reach for them when drafting:

| When you are about to… | Read this |
|---|---|
| Draft any nanopub field | `docs/forrt-form-fields.md` |
| Decide which template starts the chain | `docs/chain-decision-tree.md` |
| Choose the FORRT Claim type | `docs/claim-type-vocabulary.md` |
| Write the Quote, Study, or Outcome | `docs/verify-before-drafting.md` |
| Write a PICO or PCC question | `docs/pico-study-outcome-levels.md` |
| Touch CI / environment.yml / Dockerfile / myst.yml | `docs/cicd-conventions.md` |
| Cut a release / claim a phase done | `docs/fair4rs-checklist.md` |
| Update RO-Crate metadata after adding artefacts | `docs/ro-crate.md` |
| Need to retract / supersede / batch-publish a nanopub | `docs/programmatic-nanopubs.md` |
| Draft the launch announcement post | `docs/announcement-template.md` |
| Use this template with Cursor / Aider / Continue / non-Claude AI | `docs/ai-portability.md` |

## When in doubt

Stop and ask. The platform's value is research integrity — the cost of pausing is small; the cost of an overclaimed result or a hallucinated quote is high.
