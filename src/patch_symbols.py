#!/usr/bin/env python3
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np

PROC = Path("data/processed")
PER  = PROC / "per_sample"

full_p = PROC / "ssc_skin_qc.h5ad"
ad = sc.read_h5ad(full_p)
print("[loaded full]", full_p, ad.shape)

# Build a mapping var_name → gene_symbol by scanning per-sample files
maps = []
for fp in sorted(PER.glob("*_qc.h5ad")):
    a = sc.read_h5ad(fp, backed=None)
    symcol = next((c for c in [
        "gene_symbol","gene_names","SYMBOL","symbol","gene_name","features","feature_name"
    ] if c in a.var.columns), None)
    if symcol is None:
        del a
        continue
    m = (pd.DataFrame({
            "var": a.var_names.astype(str),
            "sym": a.var[symcol].astype(str)
        })
         .drop_duplicates(subset="var"))
    maps.append(m)
    del a

if not maps:
    raise SystemExit("No per-sample gene symbol columns found; nothing to map.")

map_df = (pd.concat(maps, ignore_index=True)
            .drop_duplicates(subset="var", keep="first"))
print("[map] entries:", len(map_df))

# IMPORTANT: use bracket access; .var is a method on DataFrame!
map_sym = pd.Series(map_df["sym"].to_numpy(), index=map_df["var"].to_numpy())

# Apply to full object
ids = ad.var_names.astype(str)
mapped = ids.to_series(index=ad.var.index).map(map_sym)

# Fallback to existing symbol column or to IDs
fallback = (ad.var["gene_symbol"].astype(str)
            if "gene_symbol" in ad.var.columns
            else pd.Series(ids, index=ad.var.index))
ad.var["gene_symbol"] = mapped.combine_first(fallback).astype(str)

# Set mt/ribo flags from symbols
names_u = ad.var["gene_symbol"].str.upper()
ad.var["mt"]   = names_u.str.startswith(("MT-","MT.","MT_"))
ad.var["ribo"] = names_u.str.startswith(("RPS","RPL"))

# Save back
ad.write_h5ad(full_p, compression="gzip")
print("[patched & saved]", full_p)
print("MT- genes detected:", int(ad.var["mt"].sum()))
print("Ribo genes detected:", int(ad.var["ribo"].sum()))
