#!/usr/bin/env python3
import argparse
from pathlib import Path
import GEOparse

RAW = Path("data/raw")

def download_geo(geo_accession: str):
    RAW.mkdir(parents=True, exist_ok=True)
    print(f"[info] downloading {geo_accession} into {RAW.resolve()}")
    _ = GEOparse.get_GEO(geo=geo_accession, destdir=str(RAW))
    print("[done] Download complete. Check data/raw for supplementary files (.h5 / .h5ad / mtx).")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--geo", type=str, required=True)
    args = ap.parse_args()
    download_geo(args.geo)
