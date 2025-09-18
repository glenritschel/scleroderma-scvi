# Systemic Sclerosis Single-Cell Atlas with scvi-tools

This repository contains a reproducible pipeline to analyze open single-cell datasets in **systemic sclerosis (SSc/CREST syndrome)** and prioritize therapeutic targets.

## Features
- Preprocessing & QC of GSE138669 (SSc skin biopsies).
- Latent embedding with scVI/totalVI (scvi-tools).
- Semi-supervised annotation with scANVI.
- Differential expression & pathway analysis.
- Target prioritization with Open Targets Platform.

## Quickstart
# 1) Clone repo
```bash
 git clone https://github.com/glenritschel/scleroderma-scvi.git
 cd scleroderma-scvi
```

# 2) Install Miniforge3 (includes mamba)
```bash
wget -O Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p ~/miniforge3
source ~/miniforge3/bin/activate
conda init     # <-- then close & reopen your terminal
```

# 3) Create the environment (adds scikit-misc so HVGs work)
```bash
mamba env create -f environment.yml
```

# 4) Activate
```bash
conda activate ssc-scvi
```

# 5) Put data
# Download raw data for systemic sclerosis (GSE138669) from GEO:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138669
# Copy *.h5 from GSE138669_RAW.tar into:
#   ~/scleroderma-scvi/data/raw/GSE138669/

# 6) Run the pipeline
```bash
python src/preprocess.py
python src/modeling.py
python src/de_analysis.py
```

# Maintenance
## 100) Reproducible env TODO
- When environment stabilizes:
  - Run: conda-lock -f environment.yml -p linux-64 -p osx-64 -p win-64
  - Commit: conda-lock.yml and platform lockfiles
- Interim: keep key versions pinned in environment.yml (scanpy, scvi-tools, anndata, gseapy).

### Current Versions
- python: 3.10.18
- scanpy: 1.11.4
- scvi-tools: 1.3.3
- anndata: 0.11.4
- gseapy: 1.1.10