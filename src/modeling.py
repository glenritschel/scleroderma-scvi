#!/usr/bin/env python3
from pathlib import Path
import scanpy as sc
import scvi

PROC = Path("data/processed")

def train_scvi(input_h5ad='ssc_skin_qc.h5ad', out_basename='ssc_skin_scvi'):
    adata = sc.read_h5ad(PROC / input_h5ad)
    if 'highly_variable' in adata.var.columns:
        adata = adata[:, adata.var['highly_variable']].copy()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata, n_latent=30)
    model.train(max_epochs=100)
    adata.obsm['X_scvi'] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep='X_scvi')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.3)
    out = PROC / f"{out_basename}.h5ad"
    adata.write(out)
    print(f"[done] wrote {out}")

if __name__ == "__main__":
    train_scvi()
