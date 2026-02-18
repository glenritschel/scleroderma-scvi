#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    ad = sc.read_h5ad(args.h5ad)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # by sample
    sc.pl.umap(ad, color="sample", show=False)
    p1 = outdir / "umap_by_sample.png"
    plt.savefig(p1, dpi=args.dpi, bbox_inches="tight")
    plt.close()

    # by condition
    sc.pl.umap(ad, color="condition", show=False)
    p2 = outdir / "umap_by_condition.png"
    plt.savefig(p2, dpi=args.dpi, bbox_inches="tight")
    plt.close()

    print("Wrote:", p1.resolve())
    print("Wrote:", p2.resolve())


if __name__ == "__main__":
    main()

