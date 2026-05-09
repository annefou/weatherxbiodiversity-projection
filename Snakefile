# Snakefile — orchestrates the Pass 1 (Iberian baseline) replication
# pipeline end-to-end + the Tier 2 (DestinE Climate DT future-climate
# projection) extension.
#
# Tier-1 rules (default) cover the four headline notebooks; each rule
# converts the jupytext .py to .ipynb and executes it in place
# (per docs/cicd-conventions.md § jupyter execute --inplace) so the
# downstream MyST Jupyter Book build picks up cell outputs.
#
# Tier-2 rules wrap the four DestinE-projection notebooks. Each Tier-2
# notebook starts with `import _tier2_guard; _tier2_guard.ensure_destine_or_skip()`
# so on a host without DestinE credentials the notebooks exit 0 with a
# clear message — Snakemake will then see the empty execute and the
# expected outputs missing. Tier-2 outputs are therefore NOT in
# `rule all`; the user invokes `snakemake --cores 1 tier2` only when
# running on the DestinE Jupyter platform.
#
# Usage:
#   snakemake --cores 1                       # Tier 1 (default)
#   snakemake --cores 1 -n                    # dry run, Tier 1
#   snakemake --cores 1 figures/main_result.png
#   snakemake --cores 1 tier2                 # Tier 2 (DestinE only)
#   snakemake --cores 1 tier2 -n              # dry run, Tier 2

NOTEBOOKS = "notebooks"
DATA = "data"
RESULTS = "results"
FIGURES = "figures"
PORT = "soroye_port"


rule all:
    input:
        f"{FIGURES}/main_result.png",
        f"{RESULTS}/headline_statistic.json",


# ---------- 01: Data download ----------
# Self-contained: the GBIF Iberia download (DOI 10.15468/dl.3frmsq) and
# the Soroye Figshare deposit (DOI 10.6084/m9.figshare.9956471). No
# credentials needed — both endpoints are publicly accessible.
rule data_download:
    output:
        f"{DATA}/gbif_dl/0006204-260423192947929.csv",
        f"{DATA}/gbif_bombus_iberia_metadata.json",
        directory("reference/Bumblebee_repo_wbombusdat/0_ClimateData"),
    log:
        f"{RESULTS}/logs/01_data_download.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 01_data_download.py && "
        "jupyter execute --inplace 01_data_download.ipynb 2>&1 | tee ../{log}"


# ---------- 02: Data clean ----------
# Wraps soroye_port/01_clean_data_iberia.py + 02_presence_absence.py +
# 03_sampling_continent.py + 04_climate_tei_pei.py.
rule data_clean:
    input:
        f"{DATA}/gbif_dl/0006204-260423192947929.csv",
        directory("reference/Bumblebee_repo_wbombusdat/0_ClimateData"),
    output:
        f"{PORT}/outputs_iberia/bombus_clean.csv",
        f"{PORT}/outputs_iberia/presence_absence.npz",
        f"{PORT}/outputs_iberia/sampling_continent.npz",
        f"{PORT}/outputs_iberia/climate_tei_pei.npz",
    log:
        f"{RESULTS}/logs/02_data_clean.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 02_data_clean.py && "
        "jupyter execute --inplace 02_data_clean.ipynb 2>&1 | tee ../{log}"


# ---------- 03: Analysis ----------
# Wraps soroye_port/05b_regression_statsmodels.py + writes the
# upstream-vs-replication comparison JSON.
rule analysis:
    input:
        f"{PORT}/outputs_iberia/presence_absence.npz",
        f"{PORT}/outputs_iberia/sampling_continent.npz",
        f"{PORT}/outputs_iberia/climate_tei_pei.npz",
    output:
        f"{RESULTS}/headline_statistic.json",
        f"{RESULTS}/glmm_coefficients.csv",
        f"{RESULTS}/posterior_bambi.nc",
        f"{PORT}/outputs_iberia/posterior_vb_summary.csv",
        f"{PORT}/outputs_iberia/dataGLMM_extinction.parquet",
    log:
        f"{RESULTS}/logs/03_analysis.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 03_analysis.py && "
        "jupyter execute --inplace 03_analysis.ipynb 2>&1 | tee ../{log}"


# ---------- 04: Figures ----------
# Writes the side-by-side upstream-v0.2.1-vs-re-run comparison plot at
# figures/main_result.png. Does not invoke soroye_port/plot_forest.py
# (which expects both global Phase-2 + Phase-3 posteriors, and we only
# run the Iberia path).
rule figures:
    input:
        f"{PORT}/outputs_iberia/posterior_vb_summary.csv",
    output:
        f"{FIGURES}/main_result.png",
    log:
        f"{RESULTS}/logs/04_figures.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 04_figures.py && "
        "jupyter execute --inplace 04_figures.ipynb 2>&1 | tee ../{log}"


# =============================================================================
# Tier 2 — DestinE Climate DT projection (opt-in, DestinE platform only)
# =============================================================================
# Aggregate target: `snakemake --cores 1 tier2`. All Tier-2 notebooks
# self-skip when DestinE credentials are absent — see
# `notebooks/_tier2_guard.py`.

DATA_DESTINE = f"{DATA}/destine"


rule tier2:
    input:
        f"{RESULTS}/projection_headline.json",
        f"{FIGURES}/projection_species_rank.png",
        f"{FIGURES}/projection_risk_map_2020_2029.png",
        f"{FIGURES}/projection_risk_map_2030_2039.png",
        f"{FIGURES}/projection_summary.png",


# ---------- 05: DestinE download ----------
# Fetch SSP3-7.0 daily fields (tmax / tmin / total precip) for two
# decade slices over the Iberian bounding box. NetCDFs land in
# `data/destine/` and are gitignored — DestinE Climate DT is
# licence-locked.
rule destine_download:
    output:
        # 2m temperature (instantaneous, 4×/day) + total precipitation
        # (accumulated, 1×/day) per horizon — DestinE Climate DT has no
        # native daily max/min, derived from t2m hourly samples in 06.
        f"{DATA_DESTINE}/destine_iberia_2020_2029_t2m.nc",
        f"{DATA_DESTINE}/destine_iberia_2020_2029_tp.nc",
        f"{DATA_DESTINE}/destine_iberia_2030_2039_t2m.nc",
        f"{DATA_DESTINE}/destine_iberia_2030_2039_tp.nc",
    log:
        f"{RESULTS}/logs/05_destine_download.log",
    shell:
        "mkdir -p $(dirname {log}) && mkdir -p " + DATA_DESTINE + " && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 05_destine_download.py && "
        "jupyter execute --inplace 05_destine_download.ipynb 2>&1 | tee ../{log}"


# ---------- 06: DestinE clean ----------
# Aggregate to the Tier-1 CEA grid + recompute TEI / PEI per species
# under future-decade climate, holding species' historical niche
# limits fixed.
rule destine_clean:
    input:
        f"{DATA_DESTINE}/destine_iberia_2020_2029_t2m.nc",
        f"{DATA_DESTINE}/destine_iberia_2020_2029_tp.nc",
        f"{DATA_DESTINE}/destine_iberia_2030_2039_t2m.nc",
        f"{DATA_DESTINE}/destine_iberia_2030_2039_tp.nc",
        f"{PORT}/outputs_iberia/climate_tei_pei.npz",
    output:
        f"{PORT}/outputs_iberia/climate_tei_pei_future_2020_2029.npz",
        f"{PORT}/outputs_iberia/climate_tei_pei_future_2030_2039.npz",
    log:
        f"{RESULTS}/logs/06_destine_clean.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 06_destine_clean.py && "
        "jupyter execute --inplace 06_destine_clean.ipynb 2>&1 | tee ../{log}"


# ---------- 07: Projection ----------
# Sample 1000 draws from the Tier-1 bambi posterior, apply to future
# TEI/PEI, write per-species posterior-mean extirpation probability +
# 95% HDI to projection_headline.json and the per-cell community-mean
# raster to results/projection_per_cell_<horizon>.npy (gitignored).
rule projection:
    input:
        f"{RESULTS}/posterior_bambi.nc",
        f"{PORT}/outputs_iberia/dataGLMM_extinction.parquet",
        f"{PORT}/outputs_iberia/sampling_continent.npz",
        f"{PORT}/outputs_iberia/climate_tei_pei_future_2020_2029.npz",
        f"{PORT}/outputs_iberia/climate_tei_pei_future_2030_2039.npz",
    output:
        f"{RESULTS}/projection_headline.json",
    log:
        f"{RESULTS}/logs/07_projection.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 07_projection.py && "
        "jupyter execute --inplace 07_projection.ipynb 2>&1 | tee ../{log}"


# ---------- 08: Projection figures ----------
# Risk-rank chart (per horizon) + per-cell community-mean risk maps.
rule projection_figures:
    input:
        f"{RESULTS}/projection_headline.json",
    output:
        f"{FIGURES}/projection_species_rank.png",
        f"{FIGURES}/projection_risk_map_2020_2029.png",
        f"{FIGURES}/projection_risk_map_2030_2039.png",
        f"{FIGURES}/projection_summary.png",
    log:
        f"{RESULTS}/logs/08_projection_figures.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "cd " + NOTEBOOKS + " && "
        "jupytext --to notebook 08_projection_figures.py && "
        "jupyter execute --inplace 08_projection_figures.ipynb 2>&1 | tee ../{log}"
