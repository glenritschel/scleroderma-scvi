#!/usr/bin/env python3
import argparse, json, datetime, pathlib
p = argparse.ArgumentParser()
p.add_argument("--version", required=True)
p.add_argument("--file", default="zenodo.json")
a = p.parse_args()
path = pathlib.Path(a.file)
meta = json.loads(path.read_text(encoding="utf-8"))
meta["version"] = a.version
meta["publication_date"] = datetime.date.today().isoformat()
path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
print(f"[ok] {a.file} updated → version={a.version}, date={meta['publication_date']}")
