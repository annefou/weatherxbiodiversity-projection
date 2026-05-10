"""Four-variant projection comparison: nside=64 vs nside=128.

Each substrate produces four candidate per-species rankings under SSP3-7.0:

  (a) Full GLMM η                       — current; within-substrate (μ, σ)
  (b) Main-effects-only η               — drop interaction terms at projection
  (c) Shared CEA-reference η            — within-substrate β, but predictors
                                          standardised against CEA pool (μ, σ)
  (d) Raw TEI > 1 cell fraction         — substrate-invariant physical metric

For each variant, we compute:
  - per-species mean η (or fraction for d) across that species' active cells
  - per-species rank (1 = highest risk)
  - Spearman rank correlation between nside=64 and nside=128

Hypothesis: (b), (c), (d) decouple substrate-coupling; (a) does not.

Usage:
    python scripts/compare_variants.py
    python scripts/compare_variants.py --horizon 2020_2029
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr

REF64 = Path("/Users/annef/Documents/ScienceLive/weatherxbiodiversity-projection")
REF128 = Path("/Users/annef/Documents/ScienceLive/weatherxbiodiversity-projection-nside128")
CEA_PARQUET = REF64 / "soroye_port" / "outputs_iberia" / "dataGLMM_extinction.parquet"

PARAM_NAMES = [
    "Intercept",
    "sc_sampling",
    "sc_TEI_bs",
    "sc_TEI_delta",
    "sc_TEI_bs:sc_TEI_delta",
    "sc_PEI_bs",
    "sc_PEI_delta",
    "sc_PEI_bs:sc_PEI_delta",
    "sc_TEI_bs:sc_PEI_bs",
    "sc_TEI_delta:sc_PEI_delta",
]
INTERACTION_TERMS = {
    "sc_TEI_bs:sc_TEI_delta",
    "sc_PEI_bs:sc_PEI_delta",
    "sc_TEI_bs:sc_PEI_bs",
    "sc_TEI_delta:sc_PEI_delta",
}


def load_substrate(root: Path, horizon: str):
    coef = pd.read_csv(
        root / "healpix_port" / "outputs_iberia" / "posterior_vb_summary.csv",
        index_col=0,
    )
    out = {"coef": {p: float(coef.loc[p, "mean"]) for p in PARAM_NAMES}}

    idata = az.from_netcdf(root / "results" / "posterior_bambi_healpix.nc")
    re_da = idata.posterior["1|species"]
    factor_dim = [d for d in re_da.dims if d not in ("chain", "draw")][0]
    out["re_levels"] = [str(x) for x in idata.posterior.coords[factor_dim].values]
    out["re_means"] = re_da.mean(dim=("chain", "draw")).values

    fut = xr.open_dataset(
        root / "healpix_port" / "outputs_iberia" / f"climate_tei_pei_future_{horizon}_healpix.nc"
    )
    out["species"] = [str(s) for s in fut["species"].values]
    out["tei_bs"] = fut["tei_bs"].values
    out["tei_delta"] = fut["tei_delta"].values
    out["pei_bs"] = fut["pei_bs"].values
    out["pei_delta"] = fut["pei_delta"].values

    parquet = pd.read_parquet(root / "healpix_port" / "outputs_iberia" / "dataGLMM_extinction.parquet")
    out["scaling_local"] = {
        col: (float(parquet[col].mean()), float(parquet[col].std(ddof=1)))
        for col in ("TEI_bs", "TEI_delta", "PEI_bs", "PEI_delta", "sampling")
    }

    sc = xr.open_dataset(root / "healpix_port" / "outputs_iberia" / "sampling_continent_healpix.nc")
    out["sampling"] = sc["sampling_total"].values

    pa = xr.open_dataset(root / "healpix_port" / "outputs_iberia" / "presence_absence_healpix.nc")
    out["prab_baseline"] = pa["prab_baseline"].values
    out["prab_recent"] = pa["prab_recent"].values
    out["pa_species"] = [str(s) for s in pa["species"].values]
    return out


def cea_reference_scaling():
    cea = pd.read_parquet(CEA_PARQUET)
    return {
        col: (float(cea[col].mean()), float(cea[col].std(ddof=1)))
        for col in ("TEI_bs", "TEI_delta", "PEI_bs", "PEI_delta", "sampling")
    }


def per_cell_predictors(sub, sp):
    if sp not in sub["species"]:
        return None
    si = sub["species"].index(sp)
    if sp not in sub["pa_species"]:
        return None
    pa_si = sub["pa_species"].index(sp)
    observed = (sub["prab_baseline"][pa_si] > 0) | (sub["prab_recent"][pa_si] > 0)
    active = ~np.isnan(sub["sampling"])
    cell_mask = observed & active
    if cell_mask.sum() == 0:
        return None
    tei_bs = sub["tei_bs"][si, cell_mask]
    tei_dl = sub["tei_delta"][si, cell_mask]
    pei_bs = sub["pei_bs"][si, cell_mask]
    pei_dl = sub["pei_delta"][si, cell_mask]
    samp = sub["sampling"][cell_mask]
    valid = (np.isfinite(tei_bs) & np.isfinite(tei_dl)
             & np.isfinite(pei_bs) & np.isfinite(pei_dl)
             & np.isfinite(samp))
    if valid.sum() == 0:
        return None
    return dict(
        n=int(valid.sum()),
        tei_bs=tei_bs[valid], tei_delta=tei_dl[valid],
        pei_bs=pei_bs[valid], pei_delta=pei_dl[valid],
        sampling=samp[valid],
    )


def species_re(sub, sp):
    if sp in sub["re_levels"]:
        return float(sub["re_means"][sub["re_levels"].index(sp)])
    return 0.0


def eta_full(sub, pred, scaling, drop_interactions=False):
    """Compute mean η across cells; optionally zero out interaction-term contributions."""
    def z(x, col):
        m, s = scaling[col]
        return (x - m) / s
    a = z(pred["tei_bs"], "TEI_bs")
    b = z(pred["tei_delta"], "TEI_delta")
    c = z(pred["pei_bs"], "PEI_bs")
    d = z(pred["pei_delta"], "PEI_delta")
    e = z(pred["sampling"], "sampling")
    cf = sub["coef"]
    eta = (
        cf["Intercept"]
        + cf["sc_sampling"] * e
        + cf["sc_TEI_bs"] * a
        + cf["sc_TEI_delta"] * b
        + cf["sc_PEI_bs"] * c
        + cf["sc_PEI_delta"] * d
    )
    if not drop_interactions:
        eta = eta + (
            cf["sc_TEI_bs:sc_TEI_delta"] * a * b
            + cf["sc_PEI_bs:sc_PEI_delta"] * c * d
            + cf["sc_TEI_bs:sc_PEI_bs"] * a * c
            + cf["sc_TEI_delta:sc_PEI_delta"] * b * d
        )
    return float(eta.mean())


def tei_exceed_fraction(pred, threshold=1.0):
    tei_future = pred["tei_bs"] + pred["tei_delta"]
    return float((tei_future > threshold).mean())


def tei_future_mean(pred):
    """Mean future TEI = TEI_bs + TEI_delta, the substrate-invariant physical signal."""
    return float((pred["tei_bs"] + pred["tei_delta"]).mean())


def rank_species(scores: dict[str, float], higher_is_riskier=True):
    items = [(sp, v) for sp, v in scores.items() if np.isfinite(v)]
    items.sort(key=lambda x: x[1], reverse=higher_is_riskier)
    return {sp: i + 1 for i, (sp, _) in enumerate(items)}


def run(horizon: str):
    sub64 = load_substrate(REF64, horizon)
    sub128 = load_substrate(REF128, horizon)
    cea_scaling = cea_reference_scaling()

    species_common = sorted(set(sub64["species"]) & set(sub128["species"]))

    rows = []
    for sp in species_common:
        p64 = per_cell_predictors(sub64, sp)
        p128 = per_cell_predictors(sub128, sp)
        if p64 is None or p128 is None:
            continue
        re64 = species_re(sub64, sp)
        re128 = species_re(sub128, sp)
        row = dict(species=sp, n_cells_64=p64["n"], n_cells_128=p128["n"])

        # (a) Full GLMM η, within-substrate scaling
        row["a_eta_64"] = eta_full(sub64, p64, sub64["scaling_local"]) + re64
        row["a_eta_128"] = eta_full(sub128, p128, sub128["scaling_local"]) + re128

        # (b) Main-effects only, within-substrate scaling
        row["b_eta_64"] = eta_full(sub64, p64, sub64["scaling_local"], drop_interactions=True) + re64
        row["b_eta_128"] = eta_full(sub128, p128, sub128["scaling_local"], drop_interactions=True) + re128

        # (c) Full GLMM η but predictors standardised against shared CEA reference
        row["c_eta_64"] = eta_full(sub64, p64, cea_scaling) + re64
        row["c_eta_128"] = eta_full(sub128, p128, cea_scaling) + re128

        # (d) Substrate-invariant: mean future TEI (continuous, in physical units)
        row["d_frac_64"] = tei_future_mean(p64)
        row["d_frac_128"] = tei_future_mean(p128)
        # (d2) Threshold metric: fraction of cells with future TEI > 0.5
        row["d2_frac_64"] = tei_exceed_fraction(p64, threshold=0.5)
        row["d2_frac_128"] = tei_exceed_fraction(p128, threshold=0.5)

        rows.append(row)

    df = pd.DataFrame(rows)

    # Compute rankings per substrate per variant
    for var in ("a", "b", "c"):
        for nside in ("64", "128"):
            df[f"{var}_rank_{nside}"] = df[f"{var}_eta_{nside}"].rank(method="min", ascending=False).astype(int)
    for var in ("d", "d2"):
        for nside in ("64", "128"):
            df[f"{var}_rank_{nside}"] = df[f"{var}_frac_{nside}"].rank(method="min", ascending=False).astype(int)

    # Spearman correlation across substrates per variant — under increasing
    # min-n-cells filters
    from scipy.stats import pearsonr
    variants = [
        ("a", "eta", "Full GLMM η (within-substrate)"),
        ("b", "eta", "Main-effects only η (within-substrate)"),
        ("c", "eta", "Full η, shared CEA reference"),
        ("d", "frac", "Mean future TEI (substrate-invariant)"),
        ("d2", "frac", "Frac TEI_future>0.5 (substrate-invariant)"),
    ]
    summary = {"horizon": horizon, "n_species_total": len(df), "by_min_n": {}}
    for min_n in (1, 5, 10, 20):
        sub = df[(df["n_cells_64"] >= min_n) & (df["n_cells_128"] >= min_n)]
        block = {"n_species": len(sub), "variants": {}}
        for var, col, label in variants:
            x = sub[f"{var}_{col}_64"].values
            y = sub[f"{var}_{col}_128"].values
            if len(sub) < 3 or x.std() == 0 or y.std() == 0:
                rho = float("nan"); pearson_r = float("nan")
            else:
                rho, _ = spearmanr(x, y)
                pearson_r, _ = pearsonr(x, y)
            top_k = min(10, len(sub))
            top64 = set(sub.nlargest(top_k, f"{var}_{col}_64")["species"])
            top128 = set(sub.nlargest(top_k, f"{var}_{col}_128")["species"])
            block["variants"][var] = {
                "label": label,
                "spearman_rho": float(rho) if rho == rho else None,
                "pearson_r": float(pearson_r) if pearson_r == pearson_r else None,
                f"top{top_k}_overlap": len(top64 & top128),
                "top_k": top_k,
            }
        summary["by_min_n"][str(min_n)] = block

    return df, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--horizon", default="2030_2039", choices=("2020_2029", "2030_2039"))
    ap.add_argument(
        "--out-dir",
        default=str(REF64 / "results"),
        help="Where to write the comparison artefacts",
    )
    args = ap.parse_args()
    df, summary = run(args.horizon)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / f"variant_comparison_{args.horizon}.csv"
    json_path = out / f"variant_comparison_{args.horizon}.json"
    df.to_csv(csv_path, index=False)
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nVariant comparison @ horizon {args.horizon}, total n={summary['n_species_total']} species")
    print(f"{'(higher Spearman = more substrate-stable; 1.00 means rankings agree perfectly)':<80}")
    for min_n, block in summary["by_min_n"].items():
        print()
        print("=" * 88)
        print(f"  Restricted to species with n_cells >= {min_n} at BOTH substrates  (n = {block['n_species']})")
        print("=" * 88)
        print(f"  {'Variant':<42} {'Spearman':>10} {'Pearson':>10} {'top∩':>10}")
        print("-" * 88)
        for var in ("a", "b", "c", "d", "d2"):
            v = block["variants"][var]
            rho = v["spearman_rho"]
            r = v["pearson_r"]
            rho_s = f"{rho:>+10.3f}" if rho is not None else f"{'nan':>10}"
            r_s = f"{r:>+10.3f}" if r is not None else f"{'nan':>10}"
            top = next(k for k in v if k.startswith("top") and k.endswith("_overlap"))
            print(f"  {v['label']:<42} {rho_s} {r_s} {v[top]:>4}/{v['top_k']}")
    print(f"\nFull table: {csv_path}")
    print(f"Summary:    {json_path}")


if __name__ == "__main__":
    main()
