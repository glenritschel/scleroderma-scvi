#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

def _symbol_series(a):
    for c in ["gene_symbol","gene_symbols","symbol","symbols","SYMBOL",
              "gene_name","gene_names","feature_name","features"]:
        if c in a.var.columns:
            return a.var[c].astype(str)
    return a.var_names.astype(str)

def compute_qc_from_full(ad_full):
    base = ad_full  # full gene counts in X
    names_u = _symbol_series(base).str.upper()

    # mito %
    base = base.copy()
    base.var["mt"] = names_u.str.startswith(("MT-","MT.","MT_"))
    sc.pp.calculate_qc_metrics(base, qc_vars=["mt"], layer=None, inplace=True)
    res = pd.DataFrame(index=base.obs_names)
    res["pct_counts_mt_recalc"] = base.obs["pct_counts_mt"].astype("float32")

    # ribo fraction
    ribo_mask = names_u.str.startswith(("RPS","RPL")).values
    if ribo_mask.any():
        X = base.X
        if sp.issparse(X):
            ribo = X[:, ribo_mask].sum(axis=1).A1
            total = X.sum(axis=1).A1
        else:
            ribo = X[:, ribo_mask].sum(axis=1)
            total = X.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total > 0, ribo / total, np.nan).astype("float32")
        res["ribo_score_from_raw"] = frac
    else:
        res["ribo_score_from_raw"] = np.nan

    return res

def safe_corr(a, b):
    df = pd.DataFrame({"a": a, "b": b}).dropna()
    return df["a"].corr(df["b"]) if len(df) > 1 else np.nan

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scvi", default="data/processed/ssc_skin_scvi.h5ad")
    ap.add_argument("--full", default="data/processed/ssc_skin_qc.h5ad")
    ap.add_argument("--save-into-scvi", action="store_true")
    args = ap.parse_args()

    ad_scvi = sc.read_h5ad(args.scvi)
    ad_full = sc.read_h5ad(args.full)
    print(f"Loaded scVI: {args.scvi} {ad_scvi.shape}  | full counts: {args.full} {ad_full.shape}")

    # compute QC from full counts
    qc = compute_qc_from_full(ad_full)

    # align to scVI cells
    qc_scvi = qc.reindex(ad_scvi.obs_names)
    ad_scvi.obs["pct_counts_mt_recalc"] = qc_scvi["pct_counts_mt_recalc"].values
    ad_scvi.obs["ribo_score_from_raw"]  = qc_scvi["ribo_score_from_raw"].values

    # compare with any existing pct_counts_mt in scVI object (may not exist)
    if "pct_counts_mt" in ad_scvi.obs:
        r = safe_corr(ad_scvi.obs["pct_counts_mt"], ad_scvi.obs["pct_counts_mt_recalc"])
        diff = ad_scvi.obs["pct_counts_mt_recalc"] - ad_scvi.obs["pct_counts_mt"]
    else:
        r = np.nan
        diff = ad_scvi.obs["pct_counts_mt_recalc"] * np.nan

    print("\n[QC] pct_counts_mt comparison (existing vs recomputed):")
    print(f"  Pearson r: {r:.3f}")
    with np.errstate(all="ignore"):
        print(f"  Δ (new-old): mean={np.nanmean(diff):.4f}, median={np.nanmedian(diff):.4f}, "
              f"p10={np.nanpercentile(diff,10):.4f}, p90={np.nanpercentile(diff,90):.4f}")

    if "leiden" in ad_scvi.obs.columns:
        print("\nTop clusters by pct_counts_mt_recalc:")
        print(ad_scvi.obs.groupby("leiden", observed=True)["pct_counts_mt_recalc"]
                     .median().sort_values(ascending=False).head(5))
        print("\nTop clusters by ribo_score_from_raw:")
        print(ad_scvi.obs.groupby("leiden", observed=True)["ribo_score_from_raw"]
                     .median().sort_values(ascending=False).head(5))

    if args.save_into_scvi:
        ad_scvi.write(args.scvi, compression="gzip")
        print(f"\n[save] wrote QC columns into {args.scvi}")

if __name__ == "__main__":
    main()
