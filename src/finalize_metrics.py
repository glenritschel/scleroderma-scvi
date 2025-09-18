#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

def pick_symbol_series(a: sc.AnnData) -> pd.Series:
    for c in ["gene_symbol","gene_symbols","symbol","symbols","SYMBOL",
              "gene_name","gene_names","feature_name","features"]:
        if c in a.var.columns:
            return a.var[c].astype(str)
    return a.var_names.astype(str)

def recompute_pct_mt_from_counts(base: sc.AnnData) -> pd.Series:
    names = pick_symbol_series(base)
    names_u = names.str.upper()
    base.var["mt"] = names_u.str.startswith(("MT-","MT.","MT_"))
    layer = "counts" if "counts" in base.layers else None
    sc.pp.calculate_qc_metrics(base, qc_vars=["mt"], layer=layer, inplace=True)
    return base.obs["pct_counts_mt"].astype("float32").copy()

def frac_ribo_from_counts(base: sc.AnnData) -> pd.Series:
    names = pick_symbol_series(base)
    names_u = names.str.upper()
    mask = names_u.str.startswith(("RPS","RPL")).values
    if not mask.any():
        return pd.Series(np.full(base.n_obs, np.nan, dtype="float32"), index=base.obs_names)
    X = base.layers["counts"] if "counts" in base.layers else base.X
    if sp.issparse(X):
        ribo = X[:, mask].sum(axis=1).A1
        total = X.sum(axis=1).A1
    else:
        ribo = X[:, mask].sum(axis=1)
        total = X.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(total > 0, ribo / total, np.nan).astype("float32")
    return pd.Series(frac, index=base.obs_names)

def main(path):
    p = Path(path)
    ad = sc.read_h5ad(p)
    assert ad.raw is not None, ".raw is missing"
    base = ad.raw.to_adata()
    # ensure counts in base.X for clarity
    if "counts" in ad.layers and base.var_names.equals(ad.var_names):
        base.layers["counts"] = ad.layers["counts"].copy()
    # recompute & persist metrics
    ad.obs["pct_counts_mt_counts"] = recompute_pct_mt_from_counts(base)  # keep old pct_counts_mt as-is
    ad.obs["frac_counts_ribo"] = frac_ribo_from_counts(base).reindex(ad.obs_names).values
    ad.obs["pct_counts_ribo"] = (100.0 * ad.obs["frac_counts_ribo"]).astype("float32")
    # provenance
    ad.uns.setdefault("qc_provenance", {})
    ad.uns["qc_provenance"].update({
        "pct_counts_mt_source": "counts, full-gene snapshot (.raw)",
        "ribo_metric": "frac_counts_ribo = sum(RPS/RPL) / total_counts (from counts)",
        "gene_symbol_col": next((c for c in ["gene_symbol","gene_symbols","symbol","symbols","SYMBOL",
                                             "gene_name","gene_names","feature_name","features"]
                                 if c in base.var.columns), "var_names"),
    })
    ad.write(p, compression="gzip")
    print(f"[done] wrote counts-based QC metrics back to {p}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/processed/ssc_skin_scvi.h5ad")
