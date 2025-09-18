#!/usr/bin/env python3
from pathlib import Path
import re
import scanpy as sc
import anndata as ad

import warnings
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=".*Variable names are not unique.*"
)

RAW = Path("data/raw")
PROC = Path("data/processed")
PROC.mkdir(parents=True, exist_ok=True)

###############################################################################
def _sample_name_from_path(p: Path) -> str:
    m = re.match(r"(GSM\d+)", p.stem)
    return m.group(1) if m else p.stem

###############################################################################
def load_any_merge(path: Path) -> ad.AnnData:
    import anndata as ad
    h5s = sorted(path.glob("**/*.h5"))
    if h5s:
        adatas = []
        batch_names = []
        for p in h5s:
            print(f"[info] loading 10x HDF5: {p}")
            a = sc.read_10x_h5(p)

            # Prefer stable gene IDs as var_names if present
            if "gene_ids" in a.var.columns:
                a.var["gene_symbol"] = a.var_names  # keep symbols
                a.var_names = a.var["gene_ids"].astype(str)

            # Ensure uniqueness
            a.var_names_make_unique()
            a.obs_names_make_unique()

            # Track sample name
            s = _sample_name_from_path(p)
            a.obs["sample"] = s
            adatas.append(a)
            batch_names.append(s)

        # Concatenate safely (prefix obs names by sample, keep outer union of genes)
        merged = ad.concat(
            adatas,
            axis=0,
            join="outer",
            label="sample",
            keys=batch_names,
            index_unique="-",   # make obs_names like <barcode>-<sample>
            merge="same"
        )
        merged.var_names_make_unique()
        merged.obs_names_make_unique()
        return merged

    # ---- SAME LOGIC for .h5ad and MTX as before, but apply the same hygiene ----
    h5ads = sorted(path.glob("**/*.h5ad"))
    if len(h5ads) > 1:
        adatas, batch_names = [], []
        for p in h5ads:
            a = sc.read_h5ad(p)
            if "gene_ids" in a.var.columns:
                a.var["gene_symbol"] = a.var_names
                a.var_names = a.var["gene_ids"].astype(str)
            a.var_names_make_unique()
            a.obs_names_make_unique()
            s = _sample_name_from_path(p)
            a.obs["sample"] = s
            adatas.append(a)
            batch_names.append(s)
        return ad.concat(adatas, axis=0, join="outer", label="sample", keys=batch_names, index_unique="-", merge="same")
    elif len(h5ads) == 1:
        a = sc.read_h5ad(h5ads[0])
        if "gene_ids" in a.var.columns:
            a.var["gene_symbol"] = a.var_names
            a.var_names = a.var["gene_ids"].astype(str)
        a.var_names_make_unique()
        a.obs_names_make_unique()
        return a

    mtx_dirs = [p for p in path.glob("**/*") if p.is_dir() and (p/"matrix.mtx").exists()]
    if mtx_dirs:
        a = sc.read_10x_mtx(mtx_dirs[0], var_names="gene_symbols", cache=True)
        # For MTX we only have symbols; just make unique
        a.var_names_make_unique()
        a.obs_names_make_unique()
        return a

    raise FileNotFoundError("No supported data found. Place .h5, .h5ad, or a 10x mtx folder under data/raw/.")

###############################################################################
def _symbol_series(adata):
    """
    Prefer gene symbols if available; fall back to var_names.
    Works whether var_names are gene_ids (e.g., Ensembl) or symbols.
    """
    if "gene_symbol" in adata.var.columns:
        return adata.var["gene_symbol"].astype(str)
    return adata.var_names.astype(str)

###############################################################################
def basic_qc(adata, min_genes=200, max_mito=0.10, min_cells_per_gene=3):
# delete below
#    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
#    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

    # Flag mito using SYMBOLS if present (robust to Ensembl var_names)
    names = _symbol_series(adata)
    adata.var["mt"] = names.str.upper().str.startswith(("MT-", "MT.", "MT_"))

    # Compute QC on raw counts (X is still counts at this stage)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["pct_counts_mt"] <= max_mito * 100)
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    return adata

###############################################################################
def run_qc(output="ssc_skin_qc.h5ad"):
    adata = load_any_merge(RAW)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    print(f"[info] raw shape: {adata.shape}")

    # Basic QC (assumes this flags mito via symbols & filters on counts)
    adata = basic_qc(adata)
    print(f"[info] post-QC (pre-freeze) shape: {adata.shape}")

    # Freeze counts + raw BEFORE any transform
    adata.layers["counts"] = adata.X.copy()        # preserve raw counts
    adata.raw = adata.copy()                       # snapshot with all genes (pre-HVG)

    # HVGs from counts (Seurat v3 expects counts)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=4000,
        layer="counts",
        batch_key="sample",   # remove if you don't have obs["sample"]
    )

    # Now normalize/log in X for downstream analyses
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Subset to HVGs
    adata = adata[:, adata.var.highly_variable].copy()

    # Sanity checks
    assert "counts" in adata.layers
    assert adata.raw is not None
    assert "highly_variable" in adata.var

    print(f"[info] post-QC (final) shape: {adata.shape}")

    out = PROC / output
    adata.write(out, compression="gzip")
    print(f"[done] wrote {out}")

###############################################################################
if __name__ == "__main__":
    run_qc()
