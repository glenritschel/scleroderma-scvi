#!/usr/bin/env bash
# Prepare a local medRxiv submission bundle.
# Usage:
#   scripts/prepare_submission.sh [PATH_TO_MANUSCRIPT.docx]
#
# Creates:
#   medrxiv_submission/
#     ├── manuscript.docx              (if provided)
#     ├── figures/                     (main figure PNGs incl. one-pager)
#     ├── supplement/                  (Supplementary tables S1..Sn)
#     └── admin/                       (text blobs for the form + license)
#
# After running, you can zip it (optional):
#   (cd medrxiv_submission && zip -r ../medrxiv_submission.zip .)

set -euo pipefail
shopt -s nullglob

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

log()    { echo "[INFO $(timestamp)] $*"; }
warn()   { echo "[WARN $(timestamp)] $*" >&2; }
err()    { echo "[ERROR $(timestamp)] $*" >&2; }

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$ROOT"

OUTDIR="medrxiv_submission"
FIGDIR="$OUTDIR/figures"
SUPPDIR="$OUTDIR/supplement"
ADMINDIR="$OUTDIR/admin"

mkdir -p "$FIGDIR" "$SUPPDIR" "$ADMINDIR"

# 1) Manuscript (optional first arg)
MANU="${1:-}"
if [[ -n "${MANU}" ]]; then
  if [[ -f "$MANU" ]]; then
    cp -f "$MANU" "$OUTDIR/manuscript.docx"
    log "Copied manuscript -> $OUTDIR/manuscript.docx"
  else
    warn "Manuscript not found: $MANU"
  fi
else
  warn "No manuscript DOCX supplied. Pass the path as the first argument."
fi

# Small helper to copy if a file exists
copy_if() {
  local src="$1" destdir="$2" destname="${3:-}"
  if [[ -f "$src" ]]; then
    if [[ -n "$destname" ]]; then
      cp -f "$src" "$destdir/$destname"
      log "Copied '$src' -> '$destdir/$destname'"
    else
      cp -f "$src" "$destdir/"
      log "Copied '$src' -> '$destdir/'"
    fi
    return 0
  fi
  return 1
}

# Helper to copy all files matching a glob (silently skip if none)
copy_glob() {
  local pattern="$1" destdir="$2"
  local matched=0
  for f in $pattern; do
    [[ -e "$f" ]] || continue
    cp -f "$f" "$destdir/"
    log "Copied '$f' -> '$destdir/'"
    matched=1
  done
  if [[ $matched -eq 0 ]]; then
    warn "No matches for $pattern"
  fi
}

# 2) Figures (add or trim to taste)
#    Keep your previously curated set; add dasatinib one-pager (two spellings, multiple locations).
FIGS=(
  "results/drug_repurposing/fig_overview.png"
  "results/figures/umap_leiden_final.png"
  "results/figures/umap_pct_counts_mt_final.png"
  "results/figures/umap_pct_counts_ribo_final.png"
  "results/figures/umap_celltype_final.png"
  "results/figures/dotplot__Fibroblast_core_dot.png"
  "results/figures/dotplot__Kerat_basal_dot.png"
  "results/drug_repurposing/fig_compound_vs_celltype_heatmap.png"
  "results/drug_repurposing/fig_final_shortlist_bar.png"
  "results/drug_repurposing/dose_time_plots/dose_resp__PD-0325901.png"
)

for f in "${FIGS[@]}"; do
  if ! copy_if "$f" "$FIGDIR"; then
    warn "Figure not found: $f"
  fi
done

# Try to include the dasatinib one-pager (various common names/locations)
DASA_CANDIDATES=(
  "results/drug_repurposing/final_bundle/onepager_dasatinib.png"
  "results/drug_repurposing/final_bundle/one_pager_dasatinib.png"
  "results/drug_repurposing/onepager_dasatinib.png"
  "results/drug_repurposing/one_pager_dasatinib.png"
)
added_dasa=0
for c in "${DASA_CANDIDATES[@]}"; do
  if copy_if "$c" "$FIGDIR" "onepager_dasatinib.png"; then
    added_dasa=1
    break
  fi
done
if [[ $added_dasa -eq 0 ]]; then
  warn "Could not find dasatinib one-pager (tried: ${DASA_CANDIDATES[*]})."
fi

# 3) Supplementary tables (now include ALL supplementary_table_*)
#    This picks up S1 mapping, S2 QC, S3 markers, S4 lineage markers, etc.
copy_glob "results/supplementary/supplementary_table_*" "$SUPPDIR"

# 4) Admin text blocks (if they exist next to repo root)
copy_if "submission_form_text.txt" "$ADMINDIR"
copy_if "cover_letter.txt"         "$ADMINDIR"
copy_if "figure_captions.txt"      "$ADMINDIR"
copy_if "medrxiv_checklist.txt"    "$ADMINDIR"
copy_if "license_CC_BY_4.0.txt"    "$ADMINDIR"

# 5) Friendly summary + zip hint
echo
log "Prepared $OUTDIR/"
if command -v tree >/dev/null 2>&1; then
  tree -a -I '.git' "$OUTDIR" || true
else
  find "$OUTDIR" -type f | sed "s|^|  |"
fi

echo
log "You can zip it if desired:"
echo "  (cd $OUTDIR && zip -r ../medrxiv_submission.zip .)"
