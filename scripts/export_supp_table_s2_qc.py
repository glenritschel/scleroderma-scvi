#!/usr/bin/env python3
# Supplementary Table S2: per-cluster QC medians + sizes
# Joins QC metrics from data/processed/ssc_skin_qc.h5ad with cluster labels
# from your scVI object (e.g., ssc_skin_scvi.withraw_umap.h5ad) via obs_names.

import pandas as pd, numpy as np
from pathlib import Path

QC_ANNDATA = Path("data/processed/ssc_skin_qc.h5ad")
CLUSTER_ANNDATA_CANDIDATES = [
    Path("data/processed/ssc_skin_scvi.withraw_umap.h5ad"),
    Path("data/processed/ssc_skin_scvi.withraw.h5ad"),
    Path("data/processed/ssc_skin_scvi.h5ad"),
]
OUTCSV = Path("results/supplementary/supplementary_table_s2_qc_by_cluster.csv")
OUTCSV.parent.mkdir(parents=True, exist_ok=True)

def pick_cluster_file():
    for p in CLUSTER_ANNDATA_CANDIDATES:
        if p.exists():
            try:
                import anndata as ad
                a = ad.read_h5ad(str(p))
                col = find_cluster_column(a.obs.columns)
                if col:
                    return p, col
            except Exception:
                pass
    raise SystemExit("Could not find a scVI AnnData with a cluster column among: "
                     + ", ".join(map(str, CLUSTER_ANNDATA_CANDIDATES)))

def find_cluster_column(cols):
    # Prefer standard Leiden names; fall back to common variants
    cols = list(cols)
    for c in cols:
        if c == "leiden": return c
    for c in cols:
        if c.startswith("leiden"): return c
    for c in ["clusters","cluster","louvain","labels","label","celltype","cell_type"]:
        if c in cols: return c
    return None

def main():
    import anndata as ad

    if not QC_ANNDATA.exists():
        raise SystemExit(f"Missing QC AnnData: {QC_ANNDATA}")
    qc_ad = ad.read_h5ad(str(QC_ANNDATA))

    # QC columns (be forgiving about names)
    obs = qc_ad.obs.copy()
    def pick_col(cands):
        for c in cands:
            if c in obs.columns: return c
        return None

    genes_col = pick_col(["n_genes_by_counts","n_genes","genes"])
    total_col = pick_col(["total_counts","umi","counts","n_counts"])
    mt_col    = pick_col(["pct_counts_mt","mt","percent_mt","pct_mt"])
    ribo_col  = pick_col(["pct_counts_ribo","ribo","percent_ribo","pct_ribo"])

    need = [genes_col, total_col, mt_col, ribo_col]
    if any(x is None for x in need):
        raise SystemExit("QC AnnData is missing one or more needed columns "
                         f"(found: genes={genes_col}, total={total_col}, mt={mt_col}, ribo={ribo_col}).")

    qc_df = obs[[genes_col, total_col, mt_col, ribo_col]].copy()
    qc_df.columns = ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"]
    qc_df.index = qc_ad.obs_names.astype(str)

    # Load cluster labels from scVI AnnData
    cl_file, cl_col = pick_cluster_file()
    scvi_ad = ad.read_h5ad(str(cl_file))
    if cl_col not in scvi_ad.obs.columns:
        raise SystemExit(f"Cluster column '{cl_col}' not in {cl_file}")
    cl = scvi_ad.obs[cl_col].astype(str).str.strip()
    cl.index = scvi_ad.obs_names.astype(str)

    # Intersect cells
    shared = qc_df.index.intersection(cl.index)
    if len(shared) == 0:
        raise SystemExit("No shared cells between QC AnnData and scVI AnnData. "
                         "Check that obs_names/barcodes align.")

    # Join and compute per-cluster medians
    df = pd.DataFrame({"cluster_id": cl.loc[shared]}).join(qc_df.loc[shared])
    g = df.groupby("cluster_id", observed=False)

    out = pd.DataFrame({
        "cluster_id": g.size().index.astype(str),
        "median_n_genes_by_counts": g["n_genes_by_counts"].median().values,
        "median_total_counts":      g["total_counts"].median().values,
        "median_pct_counts_mt":     g["pct_counts_mt"].median().values,
        "median_pct_counts_ribo":   g["pct_counts_ribo"].median().values,
        "n_cells":                  g.size().values
    })

    # Sort clusters numerically if possible
    try:
        out["_n"] = out["cluster_id"].astype(int)
        out = out.sort_values("_n").drop(columns="_n")
    except Exception:
        out = out.sort_values("cluster_id")

    out.to_csv(OUTCSV, index=False)

    # Status print
    print(f"[OK] Wrote {OUTCSV}")
    print(f"Used cluster file: {cl_file} (column: {cl_col})")
    print(f"Cells: QC={qc_df.shape[0]}  scVI={cl.shape[0]}  shared={len(shared)}")
    print(out.head(12).to_string(index=False))

if __name__ == "__main__":
    main()
