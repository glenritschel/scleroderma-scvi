#!/usr/bin/env python3
# Export Supplementary Table S3: top markers per cluster (robust to column names)

import pandas as pd, numpy as np
from pathlib import Path

CANDIDATE_DE = [
    Path("results/tables/rank_genes_groups_leiden_wilcoxon.csv"),
    Path("results/tables/de_leiden_wilcoxon.csv"),
]
MAP_CANDIDATES = [
    Path("results/supplementary/supplementary_table_s1_ensembl_to_hgnc.csv"),
    Path("results/supplementary/supp_table_s1_ensembl_to_hgnc.csv"),
]
TOP_K = 20
PADJ_MAX = 0.05
OUT = Path("results/supplementary/supplementary_table_s3_top_markers_by_cluster.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

def pick_first(paths):
    for p in paths:
        if p.exists(): return p
    return None

def normalize_de_columns(df):
    df = df.copy()
    # lower-case a copy of the names for matching, but keep originals as columns
    lower = {c: c.lower() for c in df.columns}
    inv = {v:k for k,v in lower.items()}  # map lower->original

    group_col = inv.get("group") or inv.get("cluster") or inv.get("leiden") or inv.get("cluster_id")
    gene_col  = inv.get("names") or inv.get("gene") or inv.get("feature") or inv.get("gene_id")
    logfc_col = inv.get("logfoldchanges") or inv.get("logfc") or inv.get("log2fc")
    padj_col  = inv.get("pvals_adj") or inv.get("padj") or inv.get("qval") or inv.get("qvalue") or inv.get("adj_pval") or inv.get("adj_p")

    if not all([group_col, gene_col, logfc_col, padj_col]):
        raise SystemExit(f"Could not locate required DE columns. Found -> "
                         f"group:{group_col}, gene:{gene_col}, logFC:{logfc_col}, padj:{padj_col}")

    out = df[[group_col, gene_col, logfc_col, padj_col]].copy()
    out.columns = ["cluster_id","gene_id","logFC","padj"]
    out["cluster_id"] = out["cluster_id"].astype(str).str.strip()
    out["gene_id"] = out["gene_id"].astype(str).str.strip()
    out["logFC"] = pd.to_numeric(out["logFC"], errors="coerce")
    out["padj"]  = pd.to_numeric(out["padj"],  errors="coerce")
    out = out.dropna(subset=["logFC","padj"])
    return out

def load_symbol_map():
    mpath = pick_first(MAP_CANDIDATES)
    if not mpath: return None
    m = pd.read_csv(mpath)
    m.columns = [c.lower() for c in m.columns]
    ens_col = next((c for c in ["ensembl","ensembl_id","ensembl_id_stripped","ensg","gene_id"] if c in m.columns), None)
    sym_col = next((c for c in ["symbol","hgnc_symbol","gene_symbol","symbol_upper"] if c in m.columns), None)
    if not ens_col or not sym_col: return None
    m["ensembl_stripped"] = m[ens_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    m["gene_symbol"] = m[sym_col].astype(str)
    return m[["ensembl_stripped","gene_symbol"]].drop_duplicates()

def main():
    path = pick_first(CANDIDATE_DE)
    if not path:
        raise SystemExit("No DE table found at expected locations.")
    de = pd.read_csv(path)
    de = normalize_de_columns(de)

    mapping = load_symbol_map()
    if mapping is not None:
        de["ensembl_stripped"] = de["gene_id"].str.replace(r"\.\d+$", "", regex=True)
        de = de.merge(mapping, on="ensembl_stripped", how="left")
        de["gene_display"] = de["gene_symbol"].fillna(de["gene_id"])
    else:
        de["gene_display"] = de["gene_id"]

    de["abs_logFC"] = de["logFC"].abs()
    if PADJ_MAX is not None:
        de = de[de["padj"] <= PADJ_MAX].copy()

    de = de.sort_values(["cluster_id","padj","abs_logFC"], ascending=[True, True, False])
    de["rank_in_cluster"] = de.groupby("cluster_id").cumcount() + 1
    top = de[de["rank_in_cluster"] <= TOP_K].copy()
    top["direction"] = np.where(top["logFC"] >= 0, "up", "down")

    top = top[["cluster_id","rank_in_cluster","gene_display","gene_id","logFC","padj","direction"]]
    try:
        top["_n"] = top["cluster_id"].astype(int)
        top = top.sort_values(["_n","rank_in_cluster"]).drop(columns="_n")
    except Exception:
        top = top.sort_values(["cluster_id","rank_in_cluster"])

    top.to_csv(OUT, index=False)
    print(f"[OK] Wrote {OUT}")
    print(f"Source: {path}")
    print(f"Clusters: {top['cluster_id'].nunique()}   Rows: {len(top)}")
    print(top.head(12).to_string(index=False))

if __name__ == "__main__":
    main()
