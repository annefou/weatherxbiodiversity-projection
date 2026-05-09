---
name: replication-study
description: End-to-end orchestration skill for a FORRT replication study — guides the user through Phases 0-5 of CLAUDE.md, dispatches to specialist agents (paper-analyst, replication-coder, nanopub-drafter), and tracks progress. Use this when starting a new replication or resuming one mid-phase.
---

# /replication-study

You're orchestrating a complete FORRT replication study from the freshly-initialised template through publication. This skill walks the user through the six phases in `CLAUDE.md`, dispatching to specialist agents at each stage.

## Procedure

### Step 1 — Detect current phase

Inspect repo state to determine which phase the user is currently in:

| Indicator | Phase |
|---|---|
| `{{...}}` tokens present | Phase 0 not done — run `/init-template` first |
| `paper/` empty | Phase 0 partial — ask user to drop the PDF |
| `nanopubs/drafts/00_paper_summary.md` placeholder | Phase 1 not done |
| `notebooks/03_analysis.py` is the unmodified scaffold | Phase 2 not done |
| `results/` empty or stale | Phase 3 not done |
| No GitHub release tag exists (`gh release list`) | Phase 4 not done |
| `nanopubs/PUBLISHED.md` shows "_not yet published_" | Phase 5 not done |

Report the detected phase and ask the user to confirm before proceeding.

### Step 2 — Run the relevant phase

For each phase, dispatch to the right specialist or guide the user manually:

| Phase | What happens | Handoff |
|---|---|---|
| **0 Bootstrap** | Run `/init-template`. Drop paper PDF in `paper/`. | — |
| **1 Paper analysis** | Use the `paper-analyst` agent. Output: `nanopubs/drafts/00_paper_summary.md` + `nanopubs/drafts/01_quote.md` Quoted Text. | `Agent({subagent_type: "paper-analyst"})` |
| **2 Code & data port** | Use the `replication-coder` agent. Update `environment.yml`, notebooks, Snakefile. | `Agent({subagent_type: "replication-coder"})` |
| **3 Local results** | Run `snakemake --cores 1`. Compare headline number to paper. Write `nanopubs/drafts/05_outcome.md` placeholders. | manual |
| **4 Release** | Run `docs/fair4rs-checklist.md` pre-release checklist. Cut a `gh release` with a Zenodo-formatted body. | manual |
| **5 FORRT chain** | Use the `nanopub-drafter` agent for each step. User publishes each draft on platform.sciencelive4all.org and pastes the URI into `nanopubs/PUBLISHED.md`. | `Agent({subagent_type: "nanopub-drafter"})` ×6 |

### Step 3 — Phase exit checks

Before promoting a phase as "done", verify the exit criteria from `CLAUDE.md`:

- Phase 1: `01_quote.md` has a verified verbatim quote.
- Phase 2: `snakemake --cores 1` runs end-to-end on a fresh checkout.
- Phase 3: `nanopubs/drafts/05_outcome.md` has a conclusion sentence.
- Phase 4: GitHub release page is live, Zenodo record exists, `CITATION.cff`/`codemeta.json` show the concept DOI.
- Phase 5: All six chain-step URIs are in `nanopubs/PUBLISHED.md`, embedded in the Jupyter Book, browsable from `index.md`.

If any exit criterion fails, do not promote the phase. Tell the user what's missing.

### Step 4 — Final FAIR4RS check

After Phase 5, run the checklist in `docs/fair4rs-checklist.md`. Tell the user which boxes are ticked and which need attention. The replication is "complete" when the checklist is fully ticked.

## Anti-patterns

- **Don't skip phases.** Each phase has artefacts the next one depends on.
- **Don't promote a phase without the exit criteria.** "Mostly done" leaves landmines for the user later.
- **Don't draft nanopubs before the code is verified to work.** The Outcome's evidence comes from `results/`, which doesn't exist until Phase 3 is done.
- **Don't publish nanopubs.** That's the user's job, on platform.sciencelive4all.org. You draft; they publish.

## Output

Phase-by-phase progress reports. Each phase ends with: "Phase N complete. Next: Phase N+1, which does X. Shall I continue?"
