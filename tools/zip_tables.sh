mkdir -p tables_zip_each
shopt -s nullglob
for f in results/tables/*.csv tables/*.xlsx; do
  base="$(basename "$f")"
  zip -j -q "tables_zip_each/${base}.zip" "$f"
  echo "[✓] Created tables_zip_each/${base}.zip"
done
