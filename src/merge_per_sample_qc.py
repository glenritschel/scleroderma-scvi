#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc
import anndata as ad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/processed/per_sample", help="Directory containing GSM*_qc.h5ad")
    ap.add_argument("--pattern", default="GSM*_qc.h5ad")
    ap.add_argument("--out", default="data/processed/ssc_skin_qc.h5ad")
    ap.add_argument("--join", default="outer", choices=["outer", "inner"], help="Gene join mode")
    args = ap.parse_args()

    indir = Path(args.indir)
    paths = sorted(indir.glob(args.pattern))
    if not paths:
        raise SystemExit(f"No files matched {indir}/{args.pattern}")

    ads = []
    keys = []
    for p in paths:
        gsm = p.name.replace("_qc.h5ad", "")
        a = sc.read_h5ad(p)
        a.obs["sample"] = gsm
        ads.append(a)
        keys.append(gsm)

    ad_full = ad.concat(
        ads,
        join=args.join,
        label="sample",
        keys=keys,
        index_unique="-",
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ad_full.write_h5ad(out)

    print("Wrote:", out.resolve())
    print("Cells:", ad_full.n_obs, "Genes:", ad_full.n_vars)
    print("obs columns (first 25):", list(ad_full.obs.columns)[:25])


if __name__ == "__main__":
    main()

