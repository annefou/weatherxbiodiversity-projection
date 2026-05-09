# Snakefile — orchestrates the Pass 1 (Iberian baseline) replication
# pipeline end-to-end. There are four rules, one per notebook; each
# rule converts the jupytext .py to .ipynb and executes it in place
# (per docs/cicd-conventions.md § jupyter execute --inplace) so the
# downstream MyST Jupyter Book build picks up cell outputs.
#
# The notebooks themselves are thin orchestration wrappers around the
# vendored upstream pipeline in `soroye_port/` — see soroye_port/UPSTREAM.md.
#
# Usage:
#   snakemake --cores 1            # run everything (cached steps skip)
#   snakemake --cores 1 -n         # dry run
#   snakemake --cores 1 figures/main_result.png

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
