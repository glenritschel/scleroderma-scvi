#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/processed/ssc_skin_qc.h5ad")
    ap.add_argument("--meta", default="data/metadata/GSE138669_sample_metadata.tsv")
    ap.add_argument("--cols", default="condition,chemistry", help="Comma-separated metadata columns to map")
    args = ap.parse_args()

    h5 = Path(args.h5ad)
    meta_path = Path(args.meta)
    cols = [c.strip() for c in args.cols.split(",") if c.strip()]

    ad = sc.read_h5ad(h5)
    if "sample" not in ad.obs.columns:
        raise SystemExit("obs['sample'] is missing; run merge_per_sample_qc.py first")

    meta = pd.read_csv(meta_path, sep="\t")
    if "sample" not in meta.columns:
        raise SystemExit(f"Metadata file missing required 'sample' column: {meta_path}")

    meta = meta.set_index("sample")

    for col in cols:
        if col not in meta.columns:
            print(f"[warn] column '{col}' not found in metadata; skipping")
            continue
        ad.obs[col] = ad.obs["sample"].map(meta[col])

    print("Injected columns:", [c for c in cols if c in ad.obs.columns])
    for c in cols:
        if c in ad.obs.columns:
            print(c, "nunique:", ad.obs[c].nunique(dropna=True))

    ad.write_h5ad(h5)
    print("Updated:", h5.resolve())


if __name__ == "__main__":
    main()

