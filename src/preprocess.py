#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import gc, anndata as an
from scipy import sparse
import warnings
import re

# robust sparse import (handles environments where one import path fails)
try:
    import scipy.sparse as sp
except Exception:               # very defensive, but safe
    from scipy import sparse as sp

# silence specific, known-benign messages
warnings.filterwarnings("ignore", message=".*observed=False is deprecated.*")
warnings.filterwarnings("ignore", message=".*Variable names are not unique.*")
warnings.filterwarnings("ignore", message=".*Observation names are not unique.*")

RAW_DIR   = Path("data/raw/GSE138669")
PROC_DIR  = Path("data/processed")
PER_DIR   = PROC_DIR / "per_sample"
PROC_DIR.mkdir(parents=True, exist_ok=True)
PER_DIR.mkdir(parents=True, exist_ok=True)

# --- helpers for robust symbols/QC ---
def _strip_ens_version(idx):
    return pd.Index(idx).astype(str).str.replace(r"\.\d+$", "", regex=True)

def _build_id_to_symbol_from_per_sample(per_dir: Path) -> pd.Series:
    """
    Build a map: Ensembl(without version) -> gene_symbol
    by scanning the per-sample *_qc.h5ad files.
    """
    maps = []
    for fp in sorted(per_dir.glob("*_qc.h5ad")):
        a = sc.read_h5ad(fp, backed=None)
        ens = _strip_ens_version(a.var_names)
        sym = (a.var.get("gene_symbol", a.var_names)).astype(str)
        m = pd.DataFrame({"ens": ens, "sym": sym})
        # keep real symbols first: drop rows where sym still looks like ENSG*
        m["is_symbol"] = ~m["sym"].str.upper().str.startswith("ENSG")
        m = m.sort_values("is_symbol", ascending=False).drop_duplicates("ens")
        maps.append(m[["ens","sym"]])
        del a
    if not maps:
        return pd.Series(dtype=str)
    cat = pd.concat(maps, axis=0, ignore_index=True).drop_duplicates("ens")
    return cat.set_index("ens")["sym"]

def _symbol_series(var, fallback_index):
    for c in ("gene_symbol","gene_symbols","symbol","SYMBOL","gene_name","features","gene_names"):
        if c in var.columns:
            return var[c].astype(str)
    return pd.Index(fallback_index).astype(str)

def annotate_and_qc(ad):
    # Keep raw counts
    if "counts" not in ad.layers:
        ad.layers["counts"] = ad.X.copy()

    # Gene symbols
    if "gene_symbol" not in ad.var.columns:
        ad.var["gene_symbol"] = _symbol_series(ad.var, ad.var_names)

    # Mito/ribo flags + QC
    symu = ad.var["gene_symbol"].astype(str).str.upper()
    ad.var["mt"]   = symu.str.startswith(("MT-","MT.","MT_"))
    ad.var["ribo"] = symu.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt","ribo"], layer="counts", inplace=True)

def load_one_10x(h5_path: Path, sample_id: str) -> sc.AnnData:
    ad = sc.read_10x_h5(str(h5_path))

    # stable IDs if available
    if "gene_ids" in ad.var.columns:
        ad.var_names = ad.var["gene_ids"].astype(str)
    else:
        ad.var_names = ad.var_names.astype(str)
    ad.var_names_make_unique()

    # sparse + dtype
    ad.X = ad.X.tocsr()
    ad.layers["counts"] = ad.X.copy()
    ad.X = ad.X.astype(np.float32)

    # QC metrics
    annotate_and_qc(ad)

    # ensure gene_symbol
    if "gene_symbol" in ad.var.columns:
        ad.var["gene_symbol"] = ad.var["gene_symbol"].astype(str)
    elif "gene_names" in ad.var.columns:
        ad.var["gene_symbol"] = ad.var["gene_names"].astype(str)
    else:
        for c in ("features","feature_name","SYMBOL","symbol","gene_name"):
            if c in ad.var.columns:
                ad.var["gene_symbol"] = ad.var[c].astype(str)
                break
        else:
            ad.var["gene_symbol"] = ad.var_names.astype(str)

    # filter empty/near-empty barcodes
    keep = (ad.obs["total_counts"] > 500) & (ad.obs["n_genes_by_counts"] > 200)
    if keep.sum() == 0:
        keep = ad.obs["total_counts"] > 0
    ad = ad[keep].copy()

    # unique obs names, sample tag
    ad.obs_names = sample_id + "_" + ad.obs_names.astype(str)
    ad.obs_names_make_unique()
    ad.obs["sample"] = sample_id
    return ad

def main():
    files = sorted(RAW_DIR.glob("GSM*_*feature_bc_matrix.h5"))
    if not files:
        raise FileNotFoundError(f"No 10x HDF5 files under {RAW_DIR} (expected GSM*_*feature_bc_matrix.h5)")

    # per-sample processing, save & free
    for fp in files:
        sample_id = fp.stem.split("_")[0]
        print(f"[info] reading {fp.name} → sample {sample_id}")
        ad = load_one_10x(fp, sample_id)

        # (optional) light per-sample prune
        sc.pp.filter_genes(ad, min_cells=3)

        outp = PER_DIR / f"{sample_id}_qc.h5ad"
        ad.write(outp, compression="gzip")
        print(f"[saved] {outp}  | cells: {ad.n_obs:,} genes: {ad.n_vars:,}")
        del ad; gc.collect()

    # streamed concat from disk
    print("[info] concatenating per-sample objects…")
    ad_all = None
    for fp in sorted(PER_DIR.glob("*_qc.h5ad")):
        a = sc.read_h5ad(fp, backed=None)
        a.X = a.X.tocsr()
        if ad_all is None:
            ad_all = a
        else:
            ad_all = an.concat(
                [ad_all, a],
                join="outer",
                label="sample",
                index_unique="-",
                merge="unique",
                fill_value=0,
            )
        del a; gc.collect()

    # ---- POST-CONCAT HARDENING ----
    ad_all.obs_names_make_unique()
    ad_all.var_names_make_unique()

    # Ensure CSR and counts layer
    if not sp.issparse(ad_all.X):
        ad_all.X = sp.csr_matrix(ad_all.X)
    if "counts" not in ad_all.layers:
        ad_all.layers["counts"] = ad_all.X.copy()

    # Rebuild gene_symbol from per-sample files if needed
    id2sym = _build_id_to_symbol_from_per_sample(PER_DIR)

    ens_full = _strip_ens_version(ad_all.var_names)
    gene_symbol = id2sym.reindex(ens_full).astype(object)

    # Fallbacks: keep any existing gene_symbol; otherwise var_names
    if "gene_symbol" in ad_all.var.columns:
        existing = ad_all.var["gene_symbol"].astype(str)
    else:
        existing = ad_all.var_names.astype(str)

    mask = gene_symbol.isna() | (gene_symbol.astype(str) == "") | (gene_symbol.str.upper().str.startswith("ENSG"))
    gene_symbol[mask] = existing[mask].values

    ad_all.var["gene_symbol"] = gene_symbol.astype(str)
    ad_all.var["gene_symbol_upper"] = ad_all.var["gene_symbol"].str.upper()

    # Robust mito/ribo flags from symbol column
    symu = ad_all.var["gene_symbol_upper"]
    ad_all.var["mt"]   = symu.str.startswith(("MT-","MT.","MT_"))
    ad_all.var["ribo"] = symu.str.startswith(("RPS","RPL"))

    # Counts-based QC (ensures pct_counts_mt / pct_counts_ribo exist)
    sc.pp.calculate_qc_metrics(
        ad_all,
        qc_vars=["mt","ribo"],
        layer="counts",
        percent_top=None,
        inplace=True,
    )

    # --- repair gene_symbol on merged object using one per-sample mapping ---
    per_files = sorted(PER_DIR.glob("*_qc.h5ad"))
    assert per_files, "No per-sample QC files found in data/processed/per_sample/"

    a0 = sc.read_h5ad(per_files[0])

    # Build a robust mapping from gene IDs -> symbols (strip version suffixes like .7)
    def strip_ver(s: pd.Series | pd.Index) -> pd.Series:
        return pd.Index(s.astype(str)).str.replace(r"\.\d+$", "", regex=True)

    # map from first per-sample object
    id_keys_src  = strip_ver(a0.var_names)
    id_keys_full = strip_ver(ad_all.var_names)

    # if there are duplicate Ensembl IDs, keep the first symbol
    map_sym = (
        pd.Series(a0.var["gene_symbol"].astype(str).values, index=id_keys_src)
        .groupby(level=0)
        .first()
    )

    # mapped symbols for the merged object
    mapped   = id_keys_full.to_series(index=ad_all.var.index).map(map_sym)   # Series, not Index

    # Fallback MUST be a Series aligned to ad_all.var.index (not an Index)
    fallback = pd.Series(ad_all.var_names.astype(str), index=ad_all.var.index)

    # Final gene_symbol
    ad_all.var["gene_symbol"] = mapped.combine_first(fallback).astype(str)

    print("NaNs in gene_symbol:", int(ad_all.var["gene_symbol"].isna().sum()))


    # Now set QC flags
    symu = ad_all.var["gene_symbol"].astype(str).str.upper()
    ad_all.var["mt"]   = symu.str.startswith(("MT-","MT.","MT_"))
    ad_all.var["ribo"] = symu.str.startswith(("RPS","RPL"))

    print("NaN gene_symbol fraction:", ad_all.var["gene_symbol"].isna().mean())

    # Save combined QC object
    combined_out = PROC_DIR / "ssc_skin_qc.h5ad"
    ad_all.write(combined_out, compression="gzip")
    print(f"[saved] {combined_out}  | cells: {ad_all.n_obs:,} genes: {ad_all.n_vars:,}")

if __name__ == "__main__":
    main()
