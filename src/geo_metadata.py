#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

import GEOparse
import pandas as pd

TECH_PATTERNS = re.compile(
    r"(batch|lane|pool|run|seq|flow|plate|library|prep|date|center|instrument|kit)",
    re.IGNORECASE,
)

def p(*args):
    print(*args, flush=True)

def clean_value(x):
    if pd.isna(x):
        return x
    if isinstance(x, str):
        return x.strip()
    return x

def flatten_characteristics(char_list):
    out = {}
    if not isinstance(char_list, list):
        return out
    for item in char_list:
        if ":" in item:
            key, val = item.split(":", 1)
            key = key.strip().lower().replace(" ", "_")
            out[key] = val.strip()
    return out

def build_metadata(gse_id: str) -> pd.DataFrame:
    p(f"[geo_metadata] start; python={sys.version.split()[0]}")
    p(f"[geo_metadata] downloading/parsing GEO: {gse_id}")

    gse = GEOparse.get_GEO(geo=gse_id, silent=False)

    p(f"[geo_metadata] parsed GSE; GSM count={len(gse.gsms)}")

    rows = []
    for gsm_name, gsm in gse.gsms.items():
        row = {"sample": gsm_name}
        meta = gsm.metadata

        # Basic fields
        for field in ["title", "source_name_ch1", "organism_ch1", "platform_id"]:
            if field in meta:
                v = meta[field]
                if isinstance(v, list):
                    row[field] = "; ".join(v)
                else:
                    row[field] = str(v)

        # characteristics -> columns
        if "characteristics_ch1" in meta:
            row.update(flatten_characteristics(meta["characteristics_ch1"]))

        # technical-ish keys
        for k, v in meta.items():
            if TECH_PATTERNS.search(k):
                key = k.lower()
                if isinstance(v, list):
                    row[key] = "; ".join(v)
                else:
                    row[key] = str(v)

        rows.append(row)

    df = pd.DataFrame(rows)

    for col in df.columns:
        df[col] = df[col].apply(clean_value)

    p(f"[geo_metadata] dataframe built; shape={df.shape}")
    return df

def summarize(df: pd.DataFrame):
    p("\n==============================")
    p("Column uniqueness summary:")
    p("==============================")

    uniq = sorted(
        [(c, df[c].nunique(dropna=True)) for c in df.columns],
        key=lambda x: x[1],
        reverse=True,
    )
    for c, n in uniq:
        p(f"{c:30s}  unique_values={n}")

    p("\n==============================")
    p("Candidate technical/batch-like columns:")
    p("==============================")

    for c in df.columns:
        n = df[c].nunique(dropna=True)
        if n > 1 and TECH_PATTERNS.search(c):
            vals = df[c].dropna().unique()
            vals_preview = vals if len(vals) <= 25 else list(vals[:25]) + ["..."]
            p(f"\n{c} (unique={n})")
            p(vals_preview)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse", required=True)
    args = ap.parse_args()

    df = build_metadata(args.gse)

    outdir = Path("data/metadata")
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"{args.gse}_sample_metadata.tsv"
    df.to_csv(outfile, sep="\t", index=False)

    p(f"\nWrote metadata to: {outfile.resolve()}")
    summarize(df)

if __name__ == "__main__":
    main()

