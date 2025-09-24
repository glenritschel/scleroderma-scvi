#!/usr/bin/env python3
# Supplementary Table S4: Top markers per lineage (aggregated across clusters)

import pandas as pd, numpy as np
from pathlib import Path

# ---- Inputs (robust fallbacks) ----
DE_CANDIDATES = [
    Path("results/tables/rank_genes_groups_leiden_wilcoxon.csv"),
    Path("results/tables/de_leiden_wilcoxon.csv"),
]
MAP_CANDIDATES = [
    Path("results/tables/leiden_to_celltype.csv"),
    Path("results/tables/cell_type_map.csv"),
]
SYMBOL_MAP_CANDIDATES = [
    Path("results/supplementary/supplementary_table_s1_ensembl_to_hgnc.csv"),
    Path("results/supplementary/supp_table_s1_ensembl_to_hgnc.csv"),
]

# ---- Params ----
TOP_K = 30           # markers per lineage
PADJ_MAX = 0.05      # only count support when padj <= PADJ_MAX; set None to disable
OUT = Path("results/supplementary/supplementary_table_s4_top_markers_by_lineage.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

def pick_first(paths):
    for p in paths:
        if p.exists(): return p
    return None

def normalize_de_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # map lower->original
    lower = {c: c.lower() for c in df.columns}
    inv = {v:k for k,v in lower.items()}
    group_col = inv.get("group") or inv.get("cluster") or inv.get("leiden") or inv.get("cluster_id")
    gene_col  = inv.get("names") or inv.get("gene") or inv.get("feature") or inv.get("gene_id")
    logfc_col = inv.get("logfoldchanges") or inv.get("logfc") or inv.get("log2fc")
    padj_col  = (inv.get("pvals_adj") or inv.get("padj") or inv.get("qval") or inv.get("qvalue")
                 or inv.get("adj_pval") or inv.get("adj_p"))
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

def normalize_cluster_map(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    # cluster id column
    cl = None
    for c in ["cluster_id","leiden","cluster","group"]:
        if c in df.columns:
            cl = c; break
    if cl is None: raise SystemExit("Cluster map missing a cluster id column (cluster_id/leiden/cluster/group).")

    # cell type column (a curated / final label if available)
    ct = None
    for c in ["cell_type_final","celltype_final","cell_type_curated","celltype_curated",
              "cell_type","celltype","label","labels","annotation","ann"]:
        if c in df.columns:
            ct = c; break
    if ct is None: raise SystemExit("Cluster map missing a cell type column (cell_type* / label / annotation).")

    out = df[[cl, ct]].copy()
    out.columns = ["cluster_id","cell_type"]
    out["cluster_id"] = out["cluster_id"].astype(str).str.strip()
    out["cell_type"]  = out["cell_type"].astype(str).str.strip()
    # derive lineage
    out["lineage"] = out["cell_type"].apply(to_lineage)
    return out

def to_lineage(ct: str) -> str:
    s = ct.lower()
    if "fibro" in s: return "Fibroblast"
    if "kerat" in s: return "Keratinocyte"
    if "endothel" in s: return "Endothelial"
    if "pericyte" in s or "smc" in s or "smooth" in s: return "Pericyte/SMC"
    if "myeloid" in s: return "Myeloid"
    if "mast" in s: return "Mast"
    if "plasma" in s: return "Plasma cell"
    if s.startswith("b") or " b_" in s or "b cell" in s: return "B cell"
    if s.startswith("t") or " t_" in s or "t cell" in s or "cd8" in s or "naive" in s or "cytotoxic" in s: return "T cell"
    if "cdc" in s or "pdc" in s or "dendritic" in s: return "Dendritic cell"
    return "Other"

def load_symbol_map():
    mpath = pick_first(SYMBOL_MAP_CANDIDATES)
    if not mpath: return None
    m = pd.read_csv(mpath)
    m.columns = [c.lower() for c in m.columns]
    ens = next((c for c in ["ensembl_stripped","ensembl","ensembl_id","ensg","gene_id"] if c in m.columns), None)
    sym = next((c for c in ["gene_symbol","symbol","hgnc_symbol","symbol_upper"] if c in m.columns), None)
    if ens is None or sym is None: return None
    mm = m.copy()
    mm["ensembl_stripped"] = mm[ens].astype(str).str.replace(r"\.\d+$","", regex=True)
    mm = mm[["ensembl_stripped", sym]].drop_duplicates()
    mm.columns = ["ensembl_stripped","gene_symbol"]
    return mm

def main():
    de_path = pick_first(DE_CANDIDATES)
    if not de_path:
        raise SystemExit("No DE table found.")
    de = pd.read_csv(de_path)
    de = normalize_de_columns(de)

    cmap_path = pick_first(MAP_CANDIDATES)
    if not cmap_path:
        raise SystemExit("No cluster→celltype map found (expected leiden_to_celltype.csv or cell_type_map.csv).")
    cmap = normalize_cluster_map(pd.read_csv(cmap_path))

    # attach lineage to DE via cluster_id
    de = de.merge(cmap[["cluster_id","lineage"]], on="cluster_id", how="left")
    missing_lineage = de["lineage"].isna().mean()
    if missing_lineage > 0:
        print(f"[WARN] {missing_lineage:.1%} of DE rows lack lineage (unmapped clusters). They will be dropped.")
        de = de.dropna(subset=["lineage"])

    # only consider significant rows for support (optional)
    de["abs_logFC"] = de["logFC"].abs()
    if PADJ_MAX is not None:
        de_sig = de[de["padj"] <= PADJ_MAX].copy()
    else:
        de_sig = de.copy()

    # aggregate per lineage + gene
    def agg_func(g):
        # clusters supporting (significant subset)
        sig = g[g["padj"] <= PADJ_MAX] if PADJ_MAX is not None else g
        clusters = sig["cluster_id"].astype(str).unique().tolist()
        return pd.Series({
            "best_padj": g["padj"].min(),
            "median_abs_logFC": sig["abs_logFC"].median() if len(sig) else g["abs_logFC"].median(),
            "max_abs_logFC": sig["abs_logFC"].max() if len(sig) else g["abs_logFC"].max(),
            "n_support_clusters": int(len(clusters)),
            "clusters_supporting": ",".join(clusters),
            "frac_up": float((sig["logFC"] >= 0).mean()) if len(sig) else float((g["logFC"] >= 0).mean())
        })

    agg = de_sig.groupby(["lineage","gene_id"], observed=False).apply(agg_func).reset_index()

    # attach symbols if mapping exists
    symmap = load_symbol_map()
    if symmap is not None:
        agg["ensembl_stripped"] = agg["gene_id"].str.replace(r"\.\d+$","", regex=True)
        agg = agg.merge(symmap, on="ensembl_stripped", how="left")
        agg["gene_display"] = agg["gene_symbol"].fillna(agg["gene_id"])
    else:
        agg["gene_display"] = agg["gene_id"]

    # rank per lineage: best_padj asc, then median_abs_logFC desc, then n_support desc
    agg = agg.sort_values(["lineage","best_padj","median_abs_logFC","n_support_clusters"],
                          ascending=[True, True, False, False])
    agg["rank_in_lineage"] = agg.groupby("lineage").cumcount() + 1
    top = agg[agg["rank_in_lineage"] <= TOP_K].copy()

    # tidy columns
    top = top[[
        "lineage","rank_in_lineage","gene_display","gene_id",
        "best_padj","median_abs_logFC","max_abs_logFC",
        "n_support_clusters","frac_up","clusters_supporting"
    ]]

    # sort lineages alphabetically, then rank
    top = top.sort_values(["lineage","rank_in_lineage"])

    top.to_csv(OUT, index=False)
    print(f"[OK] Wrote {OUT}")
    print(f"Source DE: {de_path}")
    print(f"Cluster map: {cmap_path}")
    print(f"Lineages: {top['lineage'].nunique()}   Rows: {len(top)}")
    print(top.head(12).to_string(index=False))

if __name__ == "__main__":
    main()
