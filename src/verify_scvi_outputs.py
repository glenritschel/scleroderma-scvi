#!/usr/bin/env python3
import argparse
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    args = ap.parse_args()

    ad = sc.read_h5ad(args.h5ad)
    print("shape:", ad.shape)
    print("obs cols include:", [c for c in ["sample", "condition", "chemistry"] if c in ad.obs.columns])
    print("obsm keys:", list(ad.obsm.keys()))
    print("X_scVI present:", "X_scVI" in ad.obsm)
    print("X_scvi present:", "X_scvi" in ad.obsm)
    print("gene_symbol in var:", "gene_symbol" in ad.var.columns)
    print("X_umap present:", "X_umap" in ad.obsm)


if __name__ == "__main__":
    main()

