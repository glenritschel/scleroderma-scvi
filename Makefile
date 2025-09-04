CONDA_ENV ?= ssc-scvi
PY ?= python
GEO ?= GSE138669
VERSION ?= 1.0.0
TAG ?= v$(VERSION)
GITHUB_REPO ?= YOURUSERNAME/scleroderma-scvi

.PHONY: help env data qc scvi de targets all clean release zip

help:
	@echo "Targets: env | data | qc | scvi | de | targets | all | release | clean"

env:
	conda env update -f environment.yml -n $(CONDA_ENV) || conda env create -f environment.yml

data:
	$(PY) src/data_download.py --geo $(GEO)

qc:
	$(PY) src/preprocess.py

scvi:
	$(PY) src/modeling.py

de:
	$(PY) src/de_analysis.py

targets:
	$(PY) src/open_targets.py

all: qc scvi de targets

zip:
	@rm -f scleroderma-scvi-$(VERSION).zip
	zip -r scleroderma-scvi-$(VERSION).zip . -x "*.git*"

release: zip
	@git add -A
	@git commit -m "Release $(TAG)" || true
	@git tag -a $(TAG) -m "Release $(TAG)" || true
	@git push --follow-tags
	@echo "Create GitHub release for $(TAG) (or use gh release create)"
