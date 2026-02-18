#!/usr/bin/env python3
import argparse
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    args = ap.parse_args()

    ad = sc.read_h5ad(args.h5ad)

    if "X_scVI" not in ad.obsm and "X_scvi" in ad.obsm:
        ad.obsm["X_scVI"] = ad.obsm["X_scvi"]
        ad.write_h5ad(args.h5ad)
        print("Aliased X_scvi -> X_scVI in", args.h5ad)
    else:
        print("No alias needed")


if __name__ == "__main__":
    main()

