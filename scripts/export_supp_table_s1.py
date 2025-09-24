# save as scripts/export_supp_table_s1.py and run: python scripts/export_supp_table_s1.py
import anndata as ad
import pandas as pd
from pathlib import Path

HVG = ad.read_h5ad("data/processed/ssc_skin_scvi.withraw.h5ad")  # or your current .withraw object
RAW = HVG.raw  # must be present

# Build symbol column in RAW if missing (strip Ensembl version and reuse HVG var map if available)
raw_var = pd.DataFrame(index=RAW.var_names.astype(str))
raw_var["ensembl_id_stripped"] = raw_var.index.str.split(".").str[0]
if "gene_symbol" in RAW.var.columns:
    raw_var["gene_symbol"] = RAW.var["gene_symbol"].astype(str)
elif "gene_symbol" in HVG.var.columns:
    # map from HVG.var onto RAW by stripped Ensembl
    hvg_map = (
        pd.DataFrame({"ens": HVG.var_names.astype(str).str.split(".").str[0],
                      "sym": HVG.var.get("gene_symbol", HVG.var_names.astype(str)).astype(str)})
        .drop_duplicates("ens")
        .set_index("ens")["sym"]
    )
    raw_var["gene_symbol"] = raw_var["ensembl_id_stripped"].map(hvg_map).fillna(raw_var["ensembl_id_stripped"])
else:
    raw_var["gene_symbol"] = raw_var["ensembl_id_stripped"]

# Presence flags
hvg_stripped = pd.Index(HVG.var_names.astype(str).str.split(".").str[0])
raw_var["present_in_hvg"] = raw_var["ensembl_id_stripped"].isin(hvg_stripped)
raw_var["present_in_raw"] = True

# Duplicated symbol flag within RAW
sym_counts = raw_var["gene_symbol"].value_counts()
raw_var["duplicated_symbol_in_raw"] = raw_var["gene_symbol"].map(sym_counts).gt(1)

# Source: best-effort label
raw_var["source_map"] = "10x_features/GTF"  # customize if you tracked exact source
raw_var["notes"] = ""

cols = ["ensembl_id_stripped","gene_symbol","source_map","present_in_hvg","present_in_raw","duplicated_symbol_in_raw","notes"]
out = Path("results") / "supplementary" / "supplementary_table_s1_ensembl_to_hgnc.csv"
out.parent.mkdir(parents=True, exist_ok=True)
raw_var[cols].reset_index(drop=True).to_csv(out, index=False)
print("Wrote", out)
