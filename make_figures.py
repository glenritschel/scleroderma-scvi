#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import scanpy as sc
except Exception as e:
    raise SystemExit(
        "This script requires scanpy. Please install it first, e.g.:\n"
        "  pip install scanpy==1.10.2 anndata==0.10.7\n"
        f"Underlying error: {e}"
    )

try:
    from adjustText import adjust_text
except Exception:
    adjust_text = None

PRIMARY = Path("data/processed/ssc_skin_scvi.withraw_umap.h5ad")
FALLBACK = Path("data/processed/ssc_skin_qc.h5ad")
OUTDIR = Path("results/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

GENE_SYMBOL_COLS = ["SYMBOL","symbol","gene_symbol","GeneSymbol","Gene","gene","genesymbol"]

def load_adata():
    path = PRIMARY if PRIMARY.exists() else FALLBACK
    print(f"[i] Loading {path} ...")
    adata = sc.read_h5ad(path)
    return adata

def ensure_umap(adata):
    if "X_umap" not in adata.obsm:
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata)
        sc.tl.umap(adata)

def get_label_key(adata) -> str:
    if "cell_type" in adata.obs:
        return "cell_type"
    if "celltype" in adata.obs:
        return "celltype"
    if "leiden" in adata.obs:
        print("[i] Using 'leiden' as labels (cell_type not present).")
        return "leiden"
    adata.obs["cell_type"] = "cells"
    return "cell_type"

def find_symbol_series(df: pd.DataFrame) -> Optional[pd.Series]:
    for col in GENE_SYMBOL_COLS:
        if col in df.columns:
            s = df[col].astype(str)
            return s
    return None

def build_symbol_table(adata) -> Tuple[pd.DataFrame, bool]:
    using_raw = False
    rows = []
    if adata.raw is not None:
        s = find_symbol_series(adata.raw.var)
        if s is not None:
            using_raw = True
            for vn, sym in zip(adata.raw.var_names, s):
                if sym and sym.strip() and sym != "nan":
                    rows.append((sym.upper(), vn))
    if not rows:
        s = find_symbol_series(adata.var)
        if s is not None:
            for vn, sym in zip(adata.var_names, s):
                if sym and sym.strip() and sym != "nan":
                    rows.append((sym.upper(), vn))
    if not rows:
        for vn in adata.var_names:
            rows.append((str(vn).upper(), vn))
    table = pd.DataFrame(rows, columns=["SYMBOL_U","var_name"]).drop_duplicates("SYMBOL_U")
    return table, using_raw

def map_symbols(adata, symbols: List[str]) -> Tuple[List[str], List[str], bool]:
    tab, using_raw = build_symbol_table(adata)
    symbols_u = [s.upper() for s in symbols]
    df = pd.DataFrame({"SYMBOL_U": symbols_u})
    df = df.merge(tab, on="SYMBOL_U", how="left")
    missing = [s for s, vn in zip(symbols, df["var_name"].tolist()) if pd.isna(vn)]
    mapped = [vn for vn in df["var_name"].tolist() if isinstance(vn, str)]
    return mapped, missing, using_raw

def _figsize_for_categories(n: int, base_w: float = 10.0, base_h: float = 6.0) -> tuple:
    h = base_h + max(0, (n - 10)) * 0.2
    return (base_w, h)

def plot_umap_with_centroids(adata, label_key: str):
    ensure_umap(adata)
    sc.set_figure_params(figsize=(8, 7), dpi=150, fontsize=12)
    ax = sc.pl.umap(adata, color=None, show=False, frameon=False, size=5, alpha=0.8)
    ax.figure.axes[0].set_xlabel("UMAP-1", fontsize=12)
    ax.figure.axes[0].set_ylabel("UMAP-2", fontsize=12)
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    groups = adata.obs[label_key].astype("category")
    texts = []
    for ct in groups.cat.categories:
        idx = groups[groups == ct].index
        if len(idx) == 0: 
            continue
        xy = umap.loc[idx].mean(0).values
        t = ax.text(
            xy[0], xy[1], str(ct),
            ha="center", va="center",
            fontsize=12, weight="bold",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7)
        )
        texts.append(t)
    if adjust_text is not None and len(texts) > 1:
        try:
            adjust_text(texts, ax=ax, only_move={'points':'y','texts':'y'})
        except Exception:
            pass
    for ext in ("png","pdf"):
        ax.figure.savefig(OUTDIR / f"fig_umap_celltype_centroid.{ext}", bbox_inches="tight")
    plt.close(ax.figure)

def _dotplot(adata, gene_symbols: List[str], groupby: str, title: str, basename: str):
    varnames, missing, using_raw = map_symbols(adata, gene_symbols)
    if missing:
        print(f"[!] Missing symbols (not found in {'raw' if using_raw else 'var'}): {', '.join(missing)}")
    if len(varnames) == 0:
        raise SystemExit(f"No genes could be mapped for {basename}. Please check symbol names or provide mapping.")

    use_raw = using_raw and (adata.raw is not None)
    sc.set_figure_params(figsize=_figsize_for_categories(adata.obs[groupby].nunique(), 10, 6), dpi=150, fontsize=12)
    dp = sc.pl.dotplot(
        adata,
        varnames,
        groupby=groupby,
        dendrogram=False,
        standard_scale=None,
        swap_axes=False,
        show=False,
        use_raw=use_raw,
        color_map="RdBu_r",
    )

    # Handle both modern DotPlot object and legacy dict-of-axes return types
    fig = None
    main_ax = None
    if hasattr(dp, "ax_dict"):
        main_ax = dp.ax_dict.get("mainplot_ax", None)
        fig = getattr(dp, "fig", None) or (main_ax.figure if main_ax is not None else plt.gcf())
        if main_ax is not None:
            main_ax.set_title(title, fontsize=14, pad=10)
            main_ax.tick_params(axis='x', labelrotation=45, labelsize=10)
        saver = getattr(dp, "savefig", None)
        for ext in ("png","pdf"):
            out = OUTDIR / f"{basename}.{ext}"
            if callable(saver):
                dp.savefig(out, bbox_inches="tight")
            else:
                fig.savefig(out, bbox_inches="tight")
            print(f"[✓] Wrote {out}")
    elif isinstance(dp, dict):
        main_ax = dp.get("mainplot_ax", None)
        fig = main_ax.figure if main_ax is not None else plt.gcf()
        if main_ax is not None:
            main_ax.set_title(title, fontsize=14, pad=10)
            main_ax.tick_params(axis='x', labelrotation=45, labelsize=10)
        for ext in ("png","pdf"):
            out = OUTDIR / f"{basename}.{ext}"
            fig.savefig(out, bbox_inches="tight")
            print(f"[✓] Wrote {out}")
    else:
        # Fallback: just save current figure
        fig = plt.gcf()
        for ext in ("png","pdf"):
            out = OUTDIR / f"{basename}.{ext}"
            fig.savefig(out, bbox_inches="tight")
            print(f"[✓] Wrote {out}")
    plt.close(fig)

def make_dotplots(adata, label_key: str):
    fibro_syms = None
    for k in ["fibro_core","Fibro_core","gene_sets","dotplot_fibro_core"]:
        if k in adata.uns and isinstance(adata.uns[k], (list, tuple)):
            fibro_syms = list(adata.uns[k]); break
        if k in adata.uns and isinstance(adata.uns[k], dict):
            for v in adata.uns[k].values():
                if isinstance(v, (list, tuple)):
                    fibro_syms = list(v); break
        if fibro_syms: break
    if fibro_syms is None:
        fibro_syms = ["COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","THY1","SFRP2","COL5A1","TAGLN"]
        print("[i] Using default fibro_core symbols:", ", ".join(fibro_syms))

    ker_syms = None
    for k in ["keratin_basal","Keratinocyte_basal","gene_sets","dotplot_keratin_basal"]:
        if k in adata.uns and isinstance(adata.uns[k], (list, tuple)):
            ker_syms = list(adata.uns[k]); break
        if k in adata.uns and isinstance(adata.uns[k], dict):
            for v in adata.uns[k].values():
                if isinstance(v, (list, tuple)):
                    ker_syms = list(v); break
        if ker_syms: break
    if ker_syms is None:
        ker_syms = ["KRT14","KRT5","TP63","ITGA6","KRT6A","KRT15","KRT16","DSC3","ITGB4","COL17A1"]
        print("[i] Using default keratin_basal symbols:", ", ".join(ker_syms))

    _dotplot(adata, fibro_syms, label_key, "Fibroblast core markers", "dotplot__fibro_core")
    _dotplot(adata, ker_syms, label_key, "Keratinocyte (basal) markers", "dotplot__keratin_basal")

def make_overview(adata, label_key: str):
    ncols = 3
    fig, axes = plt.subplots(1, ncols, figsize=(18, 5), constrained_layout=True)

    ax = axes[0]
    if "priority_final" in adata.obs.columns:
        dfp = adata.obs[[label_key, "priority_final"]].copy()
        g = dfp.groupby(label_key)["priority_final"].mean().sort_values(ascending=False)
        g.plot(kind="barh", ax=ax)
        ax.invert_yaxis()
        ax.set_title("Priority (mean by cell type)")
        ax.set_xlabel("Mean priority_final")
        ax.set_ylabel("")
    else:
        ax.text(0.5, 0.5, "priority_final not found", ha="center", va="center")
        ax.set_axis_off()

    ax = axes[1]
    rev = None
    if "rev_score_max" in adata.obsm:
        rev = adata.obsm["rev_score_max"]
        if not isinstance(rev, pd.DataFrame):
            if hasattr(rev, "toarray"):
                rev = pd.DataFrame(rev.toarray(), index=adata.obs_names)
            else:
                rev = pd.DataFrame(np.asarray(rev), index=adata.obs_names)
        if rev.shape[1] and (rev.columns is None or any([c is None for c in rev.columns])):
            rev.columns = [f"CMPD_{i}" for i in range(rev.shape[1])]
        dfh = pd.concat([adata.obs[label_key], rev], axis=1)
        agg = dfh.groupby(label_key).max(numeric_only=True)
        im = ax.imshow(agg.values, aspect="auto")
        ax.set_title("Compound × Cell type (max rev_score)")
        ax.set_yticks(range(agg.shape[0]))
        ax.set_yticklabels(agg.index, fontsize=10)
        ax.set_xticks(range(agg.shape[1]))
        ax.set_xticklabels(agg.columns, fontsize=8, rotation=90)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    else:
        ax.text(0.5, 0.5, "rev_score_max not found", ha="center", va="center")
        ax.set_axis_off()

    ax = axes[2]
    if "rev_score_max" in adata.obsm:
        dfh = pd.concat([adata.obs[label_key], rev], axis=1)
        agg = dfh.groupby(label_key).max(numeric_only=True)
        thr = np.nanpercentile(agg.values, 80)
        hits = (agg >= thr).astype(int)
        breadth = hits.sum(axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            selectivity = (1 / breadth.replace(0, np.nan)).fillna(0)
        x = breadth.values
        y = selectivity.values
        ax.scatter(x, y, s=30)
        for i, name in enumerate(agg.columns):
            if x[i] <= np.nanmedian(x) or y[i] >= np.nanmedian(y):
                ax.text(x[i], y[i], str(name), fontsize=8, ha="left", va="bottom")
        ax.set_xlabel("Breadth (# cell types hit)")
        ax.set_ylabel("Selectivity (1 / breadth)")
        ax.set_title("Selectivity vs Breadth")
    else:
        ax.text(0.5, 0.5, "rev_score_max not found", ha="center", va="center")
        ax.set_axis_off()

    for ext in ("png","pdf"):
        fig.savefig(OUTDIR / f"fig_overview.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def main():
    adata = load_adata()
    label_key = get_label_key(adata)
    adata.obs[label_key] = adata.obs[label_key].astype("category")
    plot_umap_with_centroids(adata, label_key)
    make_dotplots(adata, label_key)
    make_overview(adata, label_key)

if __name__ == "__main__":
    main()
