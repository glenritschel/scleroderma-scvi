#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--use-rep", default="", help="Override rep key; otherwise auto-detect X_scVI/X_scvi")
    ap.add_argument("--resolution", type=float, default=0.3)
    ap.add_argument("--leiden-flavor", default="leidenalg", choices=["leidenalg", "igraph"])
    args = ap.parse_args()

    ad = sc.read_h5ad(args.in_h5ad)

    if args.use_rep:
        rep = args.use_rep
    else:
        rep = "X_scVI" if "X_scVI" in ad.obsm else ("X_scvi" if "X_scvi" in ad.obsm else "")
    if not rep:
        raise SystemExit("No scVI latent found in obsm. Expected X_scVI or X_scvi.")

    sc.pp.neighbors(ad, use_rep=rep)
    sc.tl.umap(ad)

    if args.leiden_flavor == "igraph":
        sc.tl.leiden(ad, resolution=args.resolution, flavor="igraph", n_iterations=2, directed=False)
    else:
        sc.tl.leiden(ad, resolution=args.resolution)

    out = Path(args.out_h5ad)
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.write_h5ad(out)
    print("Wrote:", out.resolve())


if __name__ == "__main__":
    main()

