PY=python
CFG=conf/config.yaml

.PHONY: all prep qc map umap de enrichr reversal shortlist provenance bundle clean tidy

all: prep qc map umap de enrichr reversal shortlist provenance

prep:
	@echo ">> Preprocess / modeling"
	$(PY) src/preprocess.py --config $(CFG) || true
	$(PY) src/modeling.py   --config $(CFG) || true

qc:
	@echo ">> Verify QC & cluster medians"
	$(PY) src/verify_qc.py  --config $(CFG) || true

map:
	@echo ">> Attach .raw / symbol patch (idempotent)"
	$(PY) src/patch_symbols.py --config $(CFG) || true

umap:
	@echo ">> Neighbors/UMAP on X_scvi"
	$(PY) src/modeling.py --config $(CFG) --umap || true

de:
	@echo ">> DE (wilcoxon)"
	$(PY) src/de_analysis.py --config $(CFG) --out $(shell yq '.paths.de_csv' $(CFG) 2>/dev/null || echo results/tables/rank_genes_groups_leiden_wilcoxon.csv) || true

enrichr:
	@echo ">> Summarize Enrichr outputs (GO/Reactome/Drug)"
	$(PY) src/open_targets.py --config $(CFG) || true

reversal:
	@echo ">> Build/refresh reversal tables from Enrichr drug outputs"
	$(PY) src/finalize_metrics.py --config $(CFG) --write-reversal || true

shortlist:
	@echo ">> Day-7 scoring snapshot"
	$(PY) src/finalize_metrics.py --config $(CFG) --write-shortlist || true

provenance:
	@echo ">> Save env + provenance"
	@mkdir -p results/metadata
	conda env export > results/metadata/environment.yml || true
	pip freeze > results/metadata/pip_freeze.txt || true
	$(PY) - <<'PY'
import json, os, time, subprocess
def sh(cmd): 
  try: return subprocess.check_output(cmd, shell=True, text=True).strip()
  except: return "NA"
prov = {
 "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
 "git_sha": sh("git rev-parse HEAD"),
 "inputs": ["data/processed/ssc_skin_qc.h5ad","data/processed/ssc_skin_scvi.h5ad"],
 "outputs": ["data/processed/ssc_skin_scvi.withraw_umap.h5ad","results/tables/rank_genes_groups_leiden_wilcoxon.csv","results/shortlist_day7.csv"],
 "env": {"conda_env": "results/metadata/environment.yml","pip_freeze": "results/metadata/pip_freeze.txt"},
 "seeds": {"numpy": 0, "torch": 0, "scvi": 0}
}
open("results/metadata/provenance.json","w").write(json.dumps(prov, indent=2))
open("results/metadata/session.txt","w").write(f"{prov['generated_at']} | git {prov['git_sha']}\n")
PY

bundle:
	@echo ">> Build/share bundle (optional)"
	$(PY) src/finalize_metrics.py --config $(CFG) --build-bundle || true

clean:
	@echo ">> Non-destructive: remove transient CSVs/images (keeps final_bundle/)"
	@rm -f results/shortlist_day7.csv

tidy:
	scripts/cleanup_repo.sh
