#!/usr/bin/env bash
set -euo pipefail

# usage: scripts/cleanup_repo.sh [--dry-run] [--archive-dir DIR] [--keep-latest-deliverable|--no-keep-latest-deliverable]
DRY=0
ARCHIVE_DIR="archive"
KEEP_LATEST=1

# ---- args ----
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --archive-dir) ARCHIVE_DIR="$2"; shift 2 ;;
    --keep-latest-deliverable) KEEP_LATEST=1; shift ;;
    --no-keep-latest-deliverable) KEEP_LATEST=0; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

say() { echo "[$(date +%H:%M:%S)] $*"; }
do_mv() { if [ "$DRY" -eq 1 ]; then say "DRY mv $1 -> $2"; else mv "$1" "$2"; fi; }
do_mkdir() { if [ ! -d "$1" ]; then if [ "$DRY" -eq 1 ]; then say "DRY mkdir -p $1"; else mkdir -p "$1"; fi; fi; }
do_rm() { if [ "$DRY" -eq 1 ]; then say "DRY rm -rf $*"; else rm -rf "$@"; fi; }
do_ln() { if [ "$DRY" -eq 1 ]; then say "DRY ln -sf $1 $2"; else ln -sf "$1" "$2"; fi; }

# 1) Ensure archive dir
do_mkdir "$ARCHIVE_DIR"

# 2) Move older deliverables under results/drug_repurposing/ into archive/
DELIV_ROOT="results/drug_repurposing"
if [ -d "$DELIV_ROOT" ]; then
  # collect deliverable dirs
  DELIVS=$(find "$DELIV_ROOT" -maxdepth 1 -type d -name "deliverable_*" | sort)
  if [ -n "$DELIVS" ]; then
    if [ "$KEEP_LATEST" -eq 1 ]; then
      LATEST=$(printf "%s\n" $DELIVS | tail -n 1)
      say "Keeping latest deliverable: $LATEST"
      while IFS= read -r d; do
        [ "$d" = "$LATEST" ] && continue
        tgt="$ARCHIVE_DIR/$(basename "$d")"
        say "Archiving $d -> $tgt"
        do_mv "$d" "$tgt"
      done <<EOF
$DELIVS
EOF
    else
      while IFS= read -r d; do
        tgt="$ARCHIVE_DIR/$(basename "$d")"
        say "Archiving $d -> $tgt"
        do_mv "$d" "$tgt"
      done <<EOF
$DELIVS
EOF
    fi
  fi

  # ZIPs too
  ZIPS=$(find "$DELIV_ROOT" -maxdepth 1 -type f -name "deliverable_*.zip" | sort)
  if [ -n "$ZIPS" ]; then
    while IFS= read -r z; do
      tgt="$ARCHIVE_DIR/$(basename "$z")"
      say "Archiving $z -> $tgt"
      do_mv "$z" "$tgt"
    done <<EOF
$ZIPS
EOF
  fi
fi

# 3) Consolidate figures: move top-level ./figures/*.png into results/figures/
if [ -d "figures" ]; then
  do_mkdir "results/figures"
  found_any=0
  for f in figures/*.png; do
    [ -e "$f" ] || continue
    found_any=1
    say "Moving $f -> results/figures/"
    do_mv "$f" "results/figures/"
  done
  if [ "$found_any" -eq 1 ]; then
    rmdir figures 2>/dev/null || true
  fi
fi

# 4) Remove notebooks/results (duplicates of results/) -> archive
if [ -d "notebooks/results" ]; then
  tgt="$ARCHIVE_DIR/notebooks_results_$(date +%Y%m%d_%H%M%S)"
  say "Archiving notebooks/results -> $tgt"
  do_mv "notebooks/results" "$tgt"
fi

# 5) Ensure DE canonical filename exists (create symlink if only legacy present)
CANON="results/tables/rank_genes_groups_leiden_wilcoxon.csv"
LEGACY="results/tables/de_leiden_wilcoxon.csv"
if [ -f "$LEGACY" ] && [ ! -e "$CANON" ]; then
  say "Creating symlink $CANON -> $LEGACY"
  do_mkdir "$(dirname "$CANON")"
  do_ln "$LEGACY" "$CANON"
fi

# 6) Optional: prune sensitivity plots (large, reproducible)
SENS_DIR="results/drug_repurposing/sensitivity"
if [ -d "$SENS_DIR" ]; then
  say "Pruning sensitivity plots (keeping summary CSVs)"
  PNGS=$(find "$SENS_DIR" -type f -name "*.png")
  if [ -n "$PNGS" ]; then
    do_rm $PNGS
  fi
fi

say "Cleanup complete. (Dry-run: $DRY)"
