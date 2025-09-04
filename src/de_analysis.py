#!/usr/bin/env python3
from pathlib import Path
import scanpy as sc
import pandas as pd

PROC = Path("data/processed")
RES = Path("results/tables")
RES.mkdir(parents=True, exist_ok=True)

def de_by_cluster(h5ad='ssc_skin_scvi.h5ad'):
    adata = sc.read_h5ad(PROC / h5ad)
    if 'leiden' not in adata.obs:
        raise RuntimeError("No 'leiden' clusters found. Run modeling.py first.")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    df = sc.get.rank_genes_groups_df(adata, group=None)
    out = RES / 'de_leiden_wilcoxon.csv'
    df.to_csv(out, index=False)
    print(f"[done] wrote {out}")

if __name__ == "__main__":
    de_by_cluster()
