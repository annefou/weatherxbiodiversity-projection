# Upstream provenance for `soroye_port/`

The Python files in this directory are **vendored verbatim** from a
prior repository — they are not authored or maintained inside this
projection repo.

## Source

- Repository: [annefou/weatherxbiodiversity](https://github.com/annefou/weatherxbiodiversity)
- Version: `v0.2.1`
- Zenodo concept DOI: [10.5281/zenodo.19756173](https://doi.org/10.5281/zenodo.19756173)
- Zenodo version record: 19762723

## Licence

The upstream code is released under the MIT licence. The unmodified
upstream `LICENSE` text is preserved verbatim at
[`soroye_port/LICENSE.upstream`](LICENSE.upstream) as required by the MIT
licence's notice clause.

## Citation

If you use these scripts, cite the upstream release alongside this
projection repo:

> Fouilloux, A. (2026). WeatherXBiodiversity: Soroye et al. (2020) Replication for Iberian Bombus (v0.2.1). Zenodo. https://doi.org/10.5281/zenodo.19756173

## Why vendor instead of `pip install`?

The upstream `weatherxbiodiversity` repo is a one-shot replication
study, not a published Python package — it has no PyPI release and no
public API surface. Pinning the exact code (six scripts and one plotting
helper) inside this projection repo is the only way to guarantee that
running `snakemake --cores 1` here produces the headline statistic
`sc_TEI_delta = +0.479 [0.265, 0.694]` reported in v0.2.1.

The vendored files are otherwise unmodified. The single allowed change
is a one-line `# Vendored from …` comment at the top of each file
pointing back to this README.

## Files

| File | Role |
|---|---|
| `01_clean_data_iberia.py` | Phase 3 — clean GBIF Iberia download, apply Kerr-2015 species filter and IUCN exclusion list |
| `02_presence_absence.py` | Build 100km CEA grid, infer presence/absence per (species × period × season) |
| `03_sampling_continent.py` | Per-cell sampling-effort raster (distinct LYIDs) and continent code |
| `04_climate_tei_pei.py` | Bilinearly interpolate CRU TS 3.24.01 onto the CEA grid, compute Thermal & Precipitation Exposure Indices |
| `05b_regression_statsmodels.py` | statsmodels variational-Bayes GLMM producing `posterior_vb_summary.csv` |
| `plot_forest.py` | Forest plots of the GLMM coefficients per phase |
| `LICENSE.upstream` | Verbatim copy of the upstream MIT licence text |
