#!/usr/bin/env python3
from pathlib import Path
import re
import scanpy as sc
import anndata as ad

RAW = Path("data/raw")
PROC = Path("data/processed")
PROC.mkdir(parents=True, exist_ok=True)

def _sample_name_from_path(p: Path) -> str:
    m = re.match(r"(GSM\d+)", p.stem)
    return m.group(1) if m else p.stem

def load_any_merge(path: Path) -> ad.AnnData:
    h5s = sorted(path.glob("**/*.h5"))
    if h5s:
        adatas = []
        for p in h5s:
            adata_i = sc.read_10x_h5(p)
            adata_i.obs["sample"] = _sample_name_from_path(p)
            adatas.append(adata_i)
        if len(adatas) == 1:
            return adatas[0]
        return adatas[0].concatenate(*adatas[1:], batch_key="sample",
                                     batch_categories=[a.obs["sample"][0] for a in adatas])
    h5ads = sorted(path.glob("**/*.h5ad"))
    if len(h5ads) > 1:
        adatas = []
        for p in h5ads:
            adata_i = sc.read_h5ad(p)
            adata_i.obs["sample"] = _sample_name_from_path(p)
            adatas.append(adata_i)
        return adatas[0].concatenate(*adatas[1:], batch_key="sample",
                                     batch_categories=[a.obs["sample"][0] for a in adatas])
    elif len(h5ads) == 1:
        return sc.read_h5ad(h5ads[0])
    mtx_dirs = [p for p in path.glob("**/*") if p.is_dir() and (p/'matrix.mtx').exists()]
    if mtx_dirs:
        return sc.read_10x_mtx(mtx_dirs[0], var_names='gene_symbols', cache=True)
    raise FileNotFoundError("No supported data found. Place .h5, .h5ad, or a 10x mtx folder under data/raw/.")

def basic_qc(adata, min_genes=200, max_mito=0.10, min_cells_per_gene=3):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["pct_counts_mt"] <= max_mito * 100)
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    return adata

def run_qc(output="ssc_skin_qc.h5ad"):
    adata = load_any_merge(RAW)
    print(f"[info] raw shape: {adata.shape}")
    adata = basic_qc(adata)
    print(f"[info] post-QC shape: {adata.shape}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3")
    out = PROC / output
    adata.write(out)
    print(f"[done] wrote {out}")

if __name__ == "__main__":
    run_qc()
