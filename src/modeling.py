#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import scipy.sparse as sp

PROC = Path("data/processed")

def _symbol_series(a: sc.AnnData) -> pd.Series:
    """Prefer gene symbols if available; fall back to var_names."""
    for c in ["gene_symbol","gene_symbols","symbol","symbols","SYMBOL",
              "gene_name","gene_names","feature_name","features"]:
        if c in a.var.columns:
            return a.var[c].astype(str)
    return a.var_names.astype(str)

def _ensure_hvgs_from_counts(adata, n_top_genes=4000, batch_key=None):
    """Ensure HVGs exist, computed from COUNTS using Seurat v3."""
    if "highly_variable" in adata.var.columns and adata.var["highly_variable"].dtype == bool:
        return
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=n_top_genes,
        layer="counts",
        batch_key=batch_key if (batch_key and batch_key in adata.obs) else None,
    )

def _attach_full_gene_counts_raw(ad_hvg: sc.AnnData, ad_full: sc.AnnData) -> None:
    """
    Store a full-gene counts snapshot in .raw so QC can be recomputed later.
    Assumes ad_full.layers['counts'] holds raw counts; falls back to ad_full.X.
    """
    base = ad_full.copy()
    if "counts" in ad_full.layers:
        base.X = ad_full.layers["counts"].copy()
    # align obs to the HVG object
    base = base[ad_hvg.obs_names].copy()
    ad_hvg.raw = base  # full-gene counts snapshot

def _finalize_counts_qc(ad_hvg: sc.AnnData) -> None:
    """
    Recompute counts-based QC on the full-gene snapshot in .raw and
    write results into ad_hvg.obs.
    - pct_counts_mt  (percent of counts in mitochondrial genes)
    - frac_counts_ribo / pct_counts_ribo (fraction/percent of counts in RPS/RPL)
    """
    assert ad_hvg.raw is not None, ".raw missing; attach it before finalizing QC."
    base = ad_hvg.raw.to_adata()  # full genes, counts in X

    # Mito flag from symbols (robust to Ensembl var_names)
    names = _symbol_series(base)
    names_u = names.str.upper()
    base.var["mt"] = names_u.str.startswith(("MT-","MT.","MT_"))

    # pct_counts_mt from counts
    sc.pp.calculate_qc_metrics(base, qc_vars=["mt"], layer=None, inplace=True)
    ad_hvg.obs["pct_counts_mt"] = base.obs["pct_counts_mt"].astype("float32").values

    # Ribosomal fraction from counts (RPS/RPL)
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
        ad_hvg.obs["frac_counts_ribo"] = frac
        ad_hvg.obs["pct_counts_ribo"] = (100.0 * frac).astype("float32")
    else:
        ad_hvg.obs["frac_counts_ribo"] = np.nan
        ad_hvg.obs["pct_counts_ribo"] = np.nan

    # provenance
    ad_hvg.uns.setdefault("qc_provenance", {})
    ad_hvg.uns["qc_provenance"].update({
        "pct_counts_mt_source": "counts from .raw (full gene set)",
        "ribo_metric": "sum(RPS/RPL)/total_counts from counts",
        "gene_symbol_col": next((c for c in
                                 ["gene_symbol","gene_symbols","symbol","symbols","SYMBOL",
                                  "gene_name","gene_names","feature_name","features"]
                                 if c in base.var.columns), "var_names"),
    })

def train_scvi(
    input_h5ad: str = "ssc_skin_qc.h5ad",
    out_basename: str = "ssc_skin_scvi",
    n_top_hvgs: int = 4000,
    batch_key: str = "sample",  # change or set to None if not present
    n_latent: int = 30,
    max_epochs: int = 100,
    leiden_resolution: float = 0.3,
):
    in_path = PROC / input_h5ad
    ad_full = sc.read_h5ad(in_path)
    print(f"[info] loaded {in_path}: {ad_full.shape}")
    print(f"[info] layers={list(getattr(ad_full, 'layers', {}).keys())}, raw={'yes' if ad_full.raw is not None else 'no'}")
    assert "counts" in ad_full.layers, "layers['counts'] is required (create it in preprocess.py before normalize/log)."

    # Ensure HVGs from counts (Seurat v3)
    _ensure_hvgs_from_counts(ad_full, n_top_genes=n_top_hvgs, batch_key=batch_key)

    # Subset to HVGs for modeling
    hvgs = ad_full.var["highly_variable"].values
    ad = ad_full[:, hvgs].copy()
    print(f"[info] subset to HVGs: {ad_full.shape} → {ad.shape}")

    # scVI on counts (with optional batch)
    scvi.model.SCVI.setup_anndata(
        ad,
        layer="counts",
        batch_key=batch_key if (batch_key and batch_key in ad.obs) else None,
    )
    model = scvi.model.SCVI(ad, n_latent=n_latent)
    model.train(max_epochs=max_epochs)
    ad.obsm["X_scvi"] = model.get_latent_representation()

    # Neighbors / UMAP / Leiden
    sc.pp.neighbors(ad, use_rep="X_scvi")
    sc.tl.umap(ad)
    sc.tl.leiden(ad, resolution=leiden_resolution)

    # Attach full-gene counts snapshot, then finalize QC in-place
    _attach_full_gene_counts_raw(ad, ad_full)
    _finalize_counts_qc(ad)

    # (optional) keep a counts layer projected to HVGs on the final object
    if "counts" not in ad.layers:
        # Project counts from the full object onto HVGs if needed
        ad.layers["counts"] = ad_full[:, hvgs].layers["counts"].copy()

    out = PROC / f"{out_basename}.h5ad"
    ad.write(out, compression="gzip")
    print(f"[done] wrote {out}")

if __name__ == "__main__":
    train_scvi()
