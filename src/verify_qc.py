#!/usr/bin/env python3
"""
verify_qc.py

Loads an AnnData .h5ad file, verifies presence of counts/raw snapshots,
recomputes pct_counts_mt from counts using gene SYMBOLS when available,
computes a ribosomal score (RPS/RPL mean) from the full-gene snapshot,
and (optionally) writes results back to the file.

Usage:
  python verify_qc.py                       # defaults to data/processed/ssc_skin_scvi.h5ad
  python verify_qc.py path/to/file.h5ad
  python verify_qc.py path/to/file.h5ad --save
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import warnings

# Silence the "pkg_resources is deprecated as an API" warning
warnings.filterwarnings(
    "ignore",
    message=r".*pkg_resources is deprecated as an API.*",
    category=UserWarning,
)



def pick_symbol_series(a: sc.AnnData) -> pd.Series:
    """
    Prefer a gene symbol column if present; otherwise fall back to var_names.
    Handles common column namings found in various pipelines.
    """
    for c in [
        "gene_symbol", "gene_symbols", "symbol", "symbols", "SYMBOL",
        "gene_name", "gene_names", "feature_name", "features"
    ]:
        if c in a.var.columns:
            return a.var[c].astype(str)
    return a.var_names.astype(str)


def recompute_qc(base: sc.AnnData) -> pd.Series:
    """
    Flag mitochondrial genes by symbol and recompute pct_counts_mt on `base`.
    Uses layer='counts' if present; otherwise falls back to X.
    Returns a Series aligned to base.obs_names.
    """
    names = pick_symbol_series(base)
    names_u = names.str.upper()
    base.var["mt"] = names_u.str.startswith(("MT-", "MT.", "MT_"))
    layer = "counts" if "counts" in base.layers else None
    sc.pp.calculate_qc_metrics(base, qc_vars=["mt"], layer=layer, inplace=True)
    return base.obs["pct_counts_mt"].copy()


def ribo_score_from_counts(base: sc.AnnData) -> pd.Series:
    """
    Compute mean expression across ribosomal genes (RPS*/RPL*) on the base matrix.
    Uses X (or base.layers['counts'] if the caller sets it) and returns a Series.
    """
    names = pick_symbol_series(base)
    names_u = names.str.upper()
    ribo_mask = names_u.str.startswith(("RPS", "RPL"))

    if not ribo_mask.any():
        print("WARNING: No RPS/RPL genes found in the base matrix. Returning NaN.")
        return pd.Series(np.full(base.n_obs, np.nan), index=base.obs_names, name="ribo_score")

    Xr = base[:, ribo_mask].X
    m = Xr.mean(axis=1)
    m = m.A1 if hasattr(m, "A1") else np.asarray(m).ravel()
    return pd.Series(m, index=base.obs_names, name="ribo_score")


def main():
    p = argparse.ArgumentParser(description="Verify counts/raw QC in AnnData")
    p.add_argument(
        "h5ad",
        nargs="?",
        default="data/processed/ssc_skin_scvi.h5ad",
        help="Path to .h5ad (default: data/processed/ssc_skin_scvi.h5ad)",
    )
    p.add_argument(
        "--save", action="store_true",
        help="Write recomputed metrics back into the file"
    )
    args = p.parse_args()

    path = Path(args.h5ad)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")

    adata = sc.read_h5ad(path)
    print(f"Loaded {path}: {adata.shape}")
    print(
        f"layers={list(getattr(adata, 'layers', {}).keys())}, "
        f"raw={'yes' if adata.raw is not None else 'no'}"
    )

    # ---- Build a base object for QC (prefer full-gene snapshot) ----
    if adata.raw is not None:
        base = adata.raw.to_adata()
        # If counts layer exists on the main object *and* var matches, copy it to base
        if "counts" in adata.layers and base.var_names.equals(adata.var_names):
            base.layers["counts"] = adata.layers["counts"].copy()
    else:
        base = adata.copy()  # last resort; may be HVG-only or log1p

    # ---- Recompute QC (mito%) on base ----
    pct_mt_re = recompute_qc(base)
    adata.obs["pct_counts_mt_recalc"] = pct_mt_re.reindex(adata.obs_names).values

    # ---- Compare with existing pct_counts_mt if present ----
    if "pct_counts_mt" in adata.obs:
        old = adata.obs["pct_counts_mt"].values
        new = adata.obs["pct_counts_mt_recalc"].values
        mask = np.isfinite(old) & np.isfinite(new)
        if mask.sum() > 1:
            corr = np.corrcoef(old[mask], new[mask])[0, 1]
        else:
            corr = np.nan
        diff = new - old
        print("\n[QC] pct_counts_mt comparison (existing vs recomputed):")
        print(
            f"  Pearson r: {corr:.4f}\n"
            f"  Δ (new-old): mean={np.nanmean(diff):.4f}, median={np.nanmedian(diff):.4f}, "
            f"p10={np.nanpercentile(diff,10):.4f}, p90={np.nanpercentile(diff,90):.4f}"
        )
    else:
        print("\n[QC] pct_counts_mt was not present; only recomputed value is available.")

    # ---- Ribosomal score from counts/full genes on base ----
    ribo = ribo_score_from_counts(base)
    adata.obs["ribo_score_from_raw"] = ribo.reindex(adata.obs_names).values

    # ---- Cluster summaries if leiden exists ----
    if "leiden" in adata.obs:
        def top5(col):
            s = (adata.obs.groupby("leiden", observed=False)[col]
                 .median().sort_values(ascending=False).head(5))
            print(f"\nTop clusters by {col}:")
            print(s)

        top5("pct_counts_mt_recalc")
        top5("ribo_score_from_raw")
    else:
        print("\nℹINFO:  'leiden' not found; skipping cluster summaries.")

    if args.save:
        print(f"\n[save] Writing updated obs columns back to {path}")
        adata.write(path, compression="gzip")
    else:
        print("\n(run with --save to persist new obs columns)")


if __name__ == "__main__":
    main()
