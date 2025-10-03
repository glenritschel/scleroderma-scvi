#!/usr/bin/env bash
set -euo pipefail

# ---------------------- config you can tweak ----------------------
# Tag (default: vYYYY.MM.DD.N). Override with: ./tools/freeze_and_release.sh v0.5
TAG="${1:-v$(date +%Y.%m.%d).1}"
REPO_REMOTE="${REPO_REMOTE:-origin}"        # which git remote to push to
NBS=(
  "notebooks/09_paper_cleanup.ipynb"                  # add more notebooks here if needed
)
# Expected *minimum* artifacts (add/remove if your pipeline differs)
FIGS_EXPECT=(
  "results/figures/fig_umap_celltype_centroid.png"
  "results/figures/fig_umap_celltype_centroid.pdf"
  "results/figures/dotplot__fibro_core.png"
  "results/figures/dotplot__fibro_core.pdf"
  "results/figures/dotplot__keratin_basal.png"
  "results/figures/dotplot__keratin_basal.pdf"
  "results/figures/fig_overview.png"
  "results/figures/fig_overview.pdf"
)
TABLES_EXPECT=(
  "results/tables/supp_table_s1_gene_map.csv"
  "results/tables/supp_table_s2_cluster_qc_medians.csv"
  "results/tables/supp_table_s3_de_markers_excerpt.csv"
  "results/tables/supp_table_s4_lineage_markers.csv"
)
# What to bundle into the release zip (directories are recursed)
BUNDLE_PATHS=(
  "results/figures"
  "results/tables"
  "notebooks"
  "environment.yml"
  "README.md"
)
# ------------------------------------------------------------------

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "$ROOT" ]]; then
  echo "[!] Not inside a git repo. Run from repository root."; exit 1
fi
cd "$ROOT"

echo "==> Freeze & release for tag: $TAG"
echo "==> Repo: $ROOT"

# 0) Pre-flight
command -v git >/dev/null || { echo "[!] git not found"; exit 1; }
command -v python >/dev/null || { echo "[!] python not found"; exit 1; }
if command -v jupyter >/dev/null; then EXEC_NB="jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=0 --execute"
elif python -c "import nbconvert" 2>/dev/null; then EXEC_NB="python -m nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=0 --execute"
elif python -c "import papermill" 2>/dev/null; then EXEC_NB="papermill"  # will use: papermill in out
else
  echo "[!] Need jupyter nbconvert or papermill installed."; exit 1
fi

# 1) Re-run notebooks
echo "==> Executing notebooks"
for nb in "${NBS[@]}"; do
  [[ -f "$nb" ]] || { echo "[!] Notebook missing: $nb"; exit 1; }
  echo "   -> $nb"
  if [[ "$EXEC_NB" == "papermill" ]]; then
    papermill "$nb" "$nb"
  else
    $EXEC_NB "$nb"
  fi
done

# 2) Export environment snapshots
mkdir -p env_freeze
if command -v conda >/dev/null 2>&1; then
  conda env export --no-builds > env_freeze/conda_environment.yml || true
fi
python - <<'PY'
import sys, subprocess, pathlib, json
outdir = pathlib.Path("env_freeze"); outdir.mkdir(exist_ok=True)
# pip freeze
try:
    req = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (outdir/"requirements.txt").write_text(req)
except Exception as e:
    (outdir/"requirements.txt").write_text(f"# pip freeze failed: {e}\n")
# python & package versions
import platform
meta = {
    "python": sys.version,
    "implementation": platform.python_implementation(),
}
(outdir/"python_meta.json").write_text(json.dumps(meta, indent=2))
PY

# 3) Validate key artifacts exist
echo "==> Validating expected artifacts"
missing=0
for f in "${FIGS_EXPECT[@]}" "${TABLES_EXPECT[@]}"; do
  if [[ ! -f "$f" ]]; then echo "[miss] $f"; missing=$((missing+1)); fi
done
if (( missing > 0 )); then
  echo "[!] Missing $missing expected artifact(s). Fix before releasing, or edit the expect lists in this script."
  exit 1
fi

# 4) Freeze directory and manifest
FREEZE_DIR="freeze/${TAG}"
echo "==> Assembling freeze bundle in $FREEZE_DIR"
rm -rf "$FREEZE_DIR"
mkdir -p "$FREEZE_DIR"

# Copy artifacts
for p in "${BUNDLE_PATHS[@]}"; do
  if [[ -e "$p" ]]; then
    rsync -a --exclude='*.ipynb_checkpoints' "$p" "$FREEZE_DIR"/
  fi
done
# Always add environment snapshot
rsync -a env_freeze "$FREEZE_DIR"/

# Create MANIFEST with checksums and git info
echo "==> Writing manifest & checksums"
(
  echo "# Freeze manifest for $TAG"
  echo
  echo "## Git"
  echo "- Commit: $(git rev-parse HEAD)"
  echo "- Branch: $(git rev-parse --abbrev-ref HEAD)"
  echo "- Remote: $REPO_REMOTE ($(git remote get-url $REPO_REMOTE 2>/dev/null || echo n/a))"
  echo
  echo "## Files (SHA256)"
  find "$FREEZE_DIR" -type f | sort | while read -r f; do
    sha256sum "$f"
  done
) > "$FREEZE_DIR/MANIFEST.sha256.txt"

# 5) Zip the freeze
ZIP="release_${TAG}.zip"
echo "==> Creating zip: $ZIP"
rm -f "$ZIP"
zip -r -q "$ZIP" "$FREEZE_DIR"

# 6) Commit (optional), tag, and push tag
echo "==> Creating git tag $TAG"
if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo "[i] Tag $TAG already exists locally; skipping create."
else
  git tag -a "$TAG" -m "Freeze & release $TAG"
fi
git push "$REPO_REMOTE" --tags

# 7) Create GitHub release and upload asset
CHANGELOG="$(git log --oneline --decorate --no-merges "$(git describe --tags --abbrev=0 --tags --match 'v*' --always 2>/dev/null || echo '')"..HEAD || true)"

create_release_with_gh() {
  if gh release view "$TAG" >/dev/null 2>&1; then
    echo "[i] Release $TAG already exists on GitHub; uploading asset (may overwrite)."
    gh release upload "$TAG" "$ZIP" --clobber
  else
    gh release create "$TAG" "$ZIP" \
      --title "$TAG" \
      --notes "Automated freeze and release.

**Changelog since last tag:**
$(echo "$CHANGELOG" | sed 's/^/- /')" \
      --verify-tag
  fi
}

create_release_with_api() {
  : "${GITHUB_TOKEN:?Set GITHUB_TOKEN or use GitHub CLI 'gh' to create releases.}"
  owner_repo="$(git remote get-url "$REPO_REMOTE" | sed -E 's#.*github.com[:/](.+/[^/]+?)(\.git)?$#\1#')"
  api="https://api.github.com/repos/${owner_repo}"
  # Create release
  rel_json=$(curl -sS -H "Authorization: token $GITHUB_TOKEN" \
    -d "{\"tag_name\":\"$TAG\",\"name\":\"$TAG\",\"body\":\"Automated freeze and release.\n\n**Changelog since last tag:**\n$(printf '%s' "$CHANGELOG" | sed 's/"/\\"/g')\",\"draft\":false,\"prerelease\":false}" \
    "$api/releases")
  upload_url=$(echo "$rel_json" | python - <<'PY'
import sys, json, re
j=json.load(sys.stdin)
print(re.sub(r'\{.*\}$','',j.get("upload_url","")))
PY
)
  curl -sS -H "Authorization: token $GITHUB_TOKEN" \
       -H "Content-Type: application/zip" \
       --data-binary @"$ZIP" \
       "${upload_url}?name=$(basename "$ZIP")" >/dev/null
}

echo "==> Publishing GitHub release"
if command -v gh >/dev/null 2>&1; then
  create_release_with_gh
else
  create_release_with_api
fi

echo "==> Done. Release $TAG with asset $(basename "$ZIP") ready."
