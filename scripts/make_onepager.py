#!/usr/bin/env python3
import argparse, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_first_existing(paths):
    for p in paths:
        if Path(p).exists():
            return pd.read_csv(p)
    raise FileNotFoundError(f"None of the files exist: {paths}")

def best_p_col(df):
    for c in ["best_padj","padj","p_adj","pval_adj","FDR","padj_min","p"]:
        if c in df.columns:
            return c
    # fallbacks used by some dumps (take smallest present)
    for c in ["best_p","pval","pvalue"]:
        if c in df.columns:
            return c
    raise KeyError("Could not find an adjusted p-value column")

def main(compound):
    base = Path("results/drug_repurposing")
    fb = base / "final_bundle"

    # Load candidate metadata (MOA, dose window) if available
    meta = None
    meta_path = fb / "final_shortlist_ranked.csv"
    if meta_path.exists():
        meta = pd.read_csv(meta_path)
        # be liberal about matching the compound column
        for c in ["compound", "base_compound", "name"]:
            if c in meta.columns:
                name_col = c
                break
        else:
            name_col = None
        if name_col:
            meta_row = meta[meta[name_col].str.lower()==compound.lower()]
            meta_row = meta_row.iloc[0] if not meta_row.empty else None
        else:
            meta_row = None
    else:
        meta_row = None

    # Load top-by-cluster and/or aggregate reversal
    top_paths = [
        fb / "lincs_reversal_top15_by_cluster.csv",
        base / "lincs_reversal_top15_by_cluster.csv",
    ]
    agg_paths = [
        fb / "lincs_reversal_aggregate.csv",
        base / "lincs_reversal_aggregate.csv",
    ]
    top = load_first_existing(top_paths)
    # normalize column names
    top.columns = [c.strip().lower() for c in top.columns]
    # standardize expected columns
    if "compound" not in top.columns:
        # try 'term' containing compound name
        if "term" in top.columns:
            top["compound"] = top["term"]
        else:
            raise KeyError("No 'compound' (or 'term') column in top15 file")
    pcol = best_p_col(top)
    # try to coerce p-values numeric
    top[pcol] = pd.to_numeric(top[pcol], errors="coerce")

    # Filter to requested compound
    m = top["compound"].str.contains(compound, case=False, na=False)
    sub = top[m].copy()
    if sub.empty:
        raise SystemExit(f"No rows found for compound containing: {compound}")

    # Pick columns for display
    # unify time and dose if necessary
    time_col = "time_h" if "time_h" in sub.columns else ("time" if "time" in sub.columns else None)
    dose_col = "dose" if "dose" in sub.columns else None
    if "cluster" not in sub.columns:
        # sometimes 'group' or 'leiden'
        for alt in ["group","leiden","cluster_id"]:
            if alt in sub.columns:
                sub["cluster"] = sub[alt].astype(str)
                break
    sub["neglog10p"] = -np.log10(sub[pcol].clip(lower=1e-300))

    # Top 10 by significance
    sub = sub.sort_values(pcol).head(10).copy()
    clusters = sub["cluster"].astype(str).tolist()
    neglog10p = sub["neglog10p"].values

    # Title bits
    title = f"One-pager: {compound}"
    subtitle = ""
    if meta_row is not None:
        moa_col = "moa" if "moa" in meta_row.index else None
        if moa_col:
            subtitle += f"MOA: {meta_row[moa_col]}; "
        for dc in ["dose_min","dose_max"]:
            if dc not in meta_row.index: meta_row[dc] = np.nan
        if pd.notnull(meta_row["dose_min"]) or pd.notnull(meta_row["dose_max"]):
            subtitle += f"dose window: {meta_row['dose_min']}–{meta_row['dose_max']}"
        subtitle = subtitle.strip(" ;")

    # --------- Plot ----------
    plt.figure(figsize=(8.5, 11))  # US Letter-ish
    # Header
    y = 1.0
    plt.figtext(0.5, y-0.03, title, ha='center', va='top', fontsize=16, weight='bold')
    if subtitle:
        plt.figtext(0.5, y-0.07, subtitle, ha='center', va='top', fontsize=11)

    # Bar plot
    ax = plt.axes([0.12, 0.55, 0.76, 0.35])
    ax.barh(range(len(sub)), neglog10p)
    ax.set_yticks(range(len(sub)))
    labels = [f"cluster {c}" for c in clusters]
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(r"$- \log_{10}$(adjusted p)")
    ax.set_title("Top reversal hits by cluster")

    # Annotate dose/time
    txt_lines = []
    for _, row in sub.iterrows():
        line = f"cluster {row['cluster']}: p={row[pcol]:.2e}"
        if dose_col and not pd.isna(row.get(dose_col, np.nan)):
            line += f", dose={row[dose_col]}"
        if time_col and not pd.isna(row.get(time_col, np.nan)):
            line += f", time={row[time_col]}"
        txt_lines.append(line)
    ax2 = plt.axes([0.10, 0.24, 0.80, 0.26])
    ax2.axis('off')
    txt = "Details (top 10):\n" + "\n".join(txt_lines)
    ax2.text(0.0, 1.0, txt, ha='left', va='top', family='monospace', fontsize=9)

    # Footer
    plt.figtext(0.5, 0.06,
                "Note: LINCS reversal in immortalized cell lines is hypothesis-generating; "
                "not evidence of clinical efficacy/safety.",
                ha='center', va='center', fontsize=8)

    outdir = fb
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"onepager_{compound.replace('/','-')}.png"
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print("Wrote", outpath)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--compound", required=True, help='Compound name (case-insensitive substring match). e.g., "dasatinib"')
    args = ap.parse_args()
    main(args.compound)
