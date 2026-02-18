#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/processed/ssc_skin_qc.h5ad")
    ap.add_argument("--species", default="human")
    ap.add_argument("--chunk", type=int, default=1000)
    args = ap.parse_args()

    try:
        import mygene
    except ImportError as e:
        raise SystemExit(
            "mygene is not installed. Install with: pip install mygene\n"
            f"ImportError: {e}"
        )

    h5 = Path(args.h5ad)
    ad = sc.read_h5ad(h5)

    ens = ad.var_names.astype(str).tolist()
    mg = mygene.MyGeneInfo()

    symbols = {}
    for i in range(0, len(ens), args.chunk):
        sub = ens[i : i + args.chunk]
        res = mg.querymany(
            sub,
            scopes="ensembl.gene",
            fields="symbol",
            species=args.species,
            as_dataframe=False,
        )
        for r in res:
            q = r.get("query")
            sym = r.get("symbol")
            if q and sym:
                symbols[q] = sym

    ad.var["gene_symbol"] = [symbols.get(e, e) for e in ad.var_names.astype(str)]
    mapped = int((ad.var["gene_symbol"].astype(str) != ad.var_names.astype(str)).sum())

    ad.write_h5ad(h5)
    print("Mapped gene_symbol for", mapped, "genes")
    print("Updated:", h5.resolve())


if __name__ == "__main__":
    main()

