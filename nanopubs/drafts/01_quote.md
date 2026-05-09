# 01 — Quote-with-comment (paper-rooted chains)

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> If this is a question-rooted chain, use `01_pico.md` or `01_pcc.md` instead — see `docs/chain-decision-tree.md`.

**Form heading:** *"Annotate a paper quotation — Annotating a paper quotation with personal interpretation"*

## Field-by-field draft

### Cited DOI (text input)

Format: starts with `10.` — bare DOI, **NOT** `https://doi.org/...` form.

```
10.1126/science.aax8591
```

### Quote mode (radio button)

- [x] **Quote whole text (less than 500 characters)**
- [ ] Quote start/end *(use this if the quote exceeds 500 chars)*

### Quoted Text (textarea, required)

Verbatim from the paper PDF in `paper/`. Character-for-character. ≤ 500 chars in whole-text mode.

Source: Soroye, Newbold, Kerr (2020), *Science* 367(6478): 685, **Abstract** (p. 685, second-to-last sentence of the abstract paragraph).

```
Increasing frequency of hotter temperatures predicts species' local extinction risk, chances of colonizing a new area, and changing species richness.
```

Character count: 149 / 500.

### Comment (textarea, required)

Subtitle: *"Our interpretation or explanation of why this quotation is relevant."*

Why this quote matters and what the replication tests. Connect the paper's claim to the work this repo does. Don't repeat the quote.

```
This sentence states the mechanism — frequency of temperatures exceeding species' historical thermal tolerances drives local extinction, colonization, and richness change — that our replication propagates forward in time. The original paper validated this mechanism retrospectively for 66 Bombus species across North America and Europe at ~50 km / monthly resolution. We test its transportability to the Iberian Peninsula under SSP3-7.0, using Destination Earth Climate Digital Twin output at ~5 km daily resolution. We anchor on this sentence (rather than on the 46% / 17% North America / Europe occupancy-decline headline) because the replication tests the mechanism, not the historical numbers.
```

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 01.
