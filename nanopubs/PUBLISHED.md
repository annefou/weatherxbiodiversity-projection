# Published FORRT nanopub chain — this repository

This is the URI registry for the **canonical Iberian *Bombus* replication chain** (CEA grid + HEALPix nside=64). The chain anchors on Soroye et al. 2020 ([10.1126/science.aax8591](https://doi.org/10.1126/science.aax8591)) and confirms its TEI-based extirpation mechanism on Iberian *Bombus*: **Outcome = Validated**, CiTO citation `confirms` Soroye 2020.

For the **full three-chain constellation view** — including the sibling HEALPix nside=128 substrate extension and the methodological substrate-sensitivity diagnostic, with a graph of how all 18 nanopubs interlink — see [the constellation chapter in `weatherxbiodiversity-substrate-sensitivity`](https://annefou.github.io/weatherxbiodiversity-substrate-sensitivity/nanopubs/published).

## Chain graph

```{mermaid}
graph TB
    Soroye(["Soroye et al. 2020<br/>10.1126/science.aax8591"]):::paper
    Quote(["01 — Quote<br/>RAErLL…"]):::shared
    AIDA(["02 — AIDA: TEI_delta positive<br/>RAgb6p…"]):::shared
    Claim(["03 — Claim: statistical significance<br/>RAh7NY…"]):::shared
    Study(["04 — Study: CEA + nside=64<br/>RAybO8c8…"]):::nside64
    Outcome(["05 — Outcome: Validated<br/>RAPZMgc…"]):::nside64
    CiTO(["06 — CiTO: confirms<br/>RALbHA-…"]):::nside64
    RS(["07 — Research Software<br/>RAKH9X…"]):::nside64

    Soroye --> Quote --> AIDA --> Claim
    Claim --> Study --> Outcome --> CiTO
    Claim -.-> RS

    classDef paper fill:#fff3b0,stroke:#7a5901,stroke-width:2px,color:#000
    classDef shared fill:#e8e8e8,stroke:#444,color:#000
    classDef nside64 fill:#cfe8fc,stroke:#0072B2,color:#000
```

## URI registry

### Chain (six required steps)

| Step | Template | URI |
|---|---|---|
| 01 | Quote-with-comment | <https://w3id.org/sciencelive/np/RAErLL_QSe3e0pKBxHkUHH5v49F66fFVuS2OmYMJz02OY> |
| 02 | AIDA Sentence — TEI_delta positive on Iberian *Bombus* | <https://w3id.org/sciencelive/np/RAgb6pxwyANh-jpPdiY3H5k-fGWGgCmN72UrV_zAJcSMI> |
| 03 | FORRT Claim — statistical significance | <https://w3id.org/sciencelive/np/RAh7NYjme8dajwxnoBfbOjsd1L76LQfN-pMEajIwiRDJE> |
| 04 | Replication Study — CEA + HEALPix nside=64 | <https://w3id.org/sciencelive/np/RAybO8c8qx0p5bz9lMhMxzNsXhp0aXyd8GHnGC3i53vQY> |
| 05 | Replication Outcome — **Validated** | <https://w3id.org/sciencelive/np/RAPZMgcYbScSAXnrnSySQwZzgSA_rn-xodlMxNlwwQYY8> |
| 06 | CiTO Citation — `confirms` Soroye 2020, `extends` siblings | <https://w3id.org/sciencelive/np/RALbHA-r6wIFOFPFlfIpwYqJEpzCFqeJ082iChgdfvhNM> |

### Optional layer

| Step | Template | URI |
|---|---|---|
| 07 | Research Software | <https://w3id.org/sciencelive/np/RAKH9XeZn3CUr9WaFKMC3O2pT_HJJ96c3jTa6v6dWEE3c> |

A Research Synthesis nanopub is not published from this chain — it is published from the sibling `weatherxbiodiversity-substrate-sensitivity` repository, which combines this Outcome with the nside=128 sibling Outcome and the substrate-sensitivity diagnostic Outcome into a single cross-chain synthesis.

## Sibling chains

- **`weatherxbiodiversity-projection-nside128`** — substrate extension at HEALPix nside=128. [Zenodo](https://doi.org/10.5281/zenodo.20113780) · [Jupyter Book](https://annefou.github.io/weatherxbiodiversity-projection-nside128/)
- **`weatherxbiodiversity-substrate-sensitivity`** — methodological diagnostic + cross-chain Research Synthesis. [Zenodo](https://doi.org/10.5281/zenodo.20113786) · [Jupyter Book](https://annefou.github.io/weatherxbiodiversity-substrate-sensitivity/)

## How to view a nanopub

Open any URI directly in your browser. The Science Live viewer renders the four named graphs (Head, Assertion, Provenance, PublicationInfo). If a direct link doesn't resolve, wrap the URI:

```
https://platform.sciencelive4all.org/np/?uri=<full-URI>
```

Nanopubs are immutable once published. To correct a published nanopub, publish a retraction or supersession (see `docs/programmatic-nanopubs.md`).
