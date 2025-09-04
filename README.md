# Systemic Sclerosis Single-Cell Atlas (Scanpy + scVI)

# Systemic Sclerosis Single-Cell Atlas (Scanpy + scVI)

Reproducible pipeline to analyze open single-cell datasets in **systemic sclerosis (SSc/CREST)** and prioritize therapeutic targets with **Open Targets**.

- Handles **multiple 10x HDF5 (`.h5`) files** (e.g., from `GSE138669_RAW.tar`) or `.h5ad` or 10x MTX.
- QC → scVI (UMAP/Leiden) → DE → target integration.
- WSL2/Ubuntu friendly. Uses **Miniforge3** (includes `mamba`).

## Quickstart
```bash
wget -O Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p ~/miniforge3
source ~/miniforge3/bin/activate && conda init
mamba env create -f environment.yml
conda activate ssc-scvi
python src/preprocess.py
python src/modeling.py
python src/de_analysis.py

