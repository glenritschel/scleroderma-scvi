#!/usr/bin/env bash
set -euo pipefail

GSE="${GSE:-GSE138669}"

# Inputs / outputs
PER_SAMPLE_DIR="${PER_SAMPLE_DIR:-data/processed/per_sample}"
MERGED_H5AD="${MERGED_H5AD:-data/processed/ssc_skin_qc.h5ad}"
META_TSV="${META_TSV:-data/metadata/${GSE}_sample_metadata.tsv}"

SCVI_OUT_BASE="${SCVI_OUT_BASE:-ssc_skin_scvi_batch_sample}"
SCVI_OUT_H5AD="${SCVI_OUT_H5AD:-data/processed/${SCVI_OUT_BASE}.h5ad}"
UMAP_OUT_H5AD="${UMAP_OUT_H5AD:-data/processed/${SCVI_OUT_BASE}_umap.h5ad}"

FIG_DIR="${FIG_DIR:-figures}"

echo "[run_all] 1) GEO metadata export"
python -u src/geo_metadata.py --gse "${GSE}"

echo "[run_all] 2) Merge per-sample QC into merged h5ad"
python -u src/merge_per_sample_qc.py \
  --indir "${PER_SAMPLE_DIR}" \
  --out "${MERGED_H5AD}" \
  --join outer

echo "[run_all] 3) Inject GEO metadata into merged h5ad (condition, chemistry)"
python -u src/inject_geo_metadata.py \
  --h5ad "${MERGED_H5AD}" \
  --meta "${META_TSV}" \
  --cols "condition,chemistry"

echo "[run_all] 4) Add gene_symbol mapping (requires mygene)"
python -u src/add_gene_symbols_mygene.py \
  --h5ad "${MERGED_H5AD}" \
  --species human

echo "[run_all] 5) Train scVI (batch_key=sample, covariates=condition,chemistry, epochs=200)"
python -u src/modeling.py \
  --input "${MERGED_H5AD}" \
  --out-basename "${SCVI_OUT_BASE}" \
  --batch-key sample \
  --cat-cov condition,chemistry \
  --epochs 200

echo "[run_all] 6) Ensure X_scVI key exists (alias if needed)"
python -u src/alias_scvi_key.py --h5ad "${SCVI_OUT_H5AD}"

echo "[run_all] 7) Compute neighbors/UMAP/Leiden on scVI latent and write UMAP h5ad"
python -u src/make_umap_leiden.py \
  --in-h5ad "${SCVI_OUT_H5AD}" \
  --out-h5ad "${UMAP_OUT_H5AD}" \
  --resolution 0.3

echo "[run_all] 8) Export proof figures"
python -u src/export_umap_figures.py \
  --h5ad "${UMAP_OUT_H5AD}" \
  --outdir "${FIG_DIR}" \
  --dpi 300

echo "[run_all] 9) Verify outputs"
python -u src/verify_scvi_outputs.py --h5ad "${SCVI_OUT_H5AD}"
python -u src/verify_scvi_outputs.py --h5ad "${UMAP_OUT_H5AD}"

echo "[run_all] DONE"
echo "  scVI:  ${SCVI_OUT_H5AD}"
echo "  UMAP:  ${UMAP_OUT_H5AD}"
echo "  figs:  ${FIG_DIR}/umap_by_sample.png, ${FIG_DIR}/umap_by_condition.png"

