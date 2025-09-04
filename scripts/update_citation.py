#!/usr/bin/env python3
import argparse, re, pathlib
p = argparse.ArgumentParser()
p.add_argument("--version", required=True)
p.add_argument("--doi", default=None)
p.add_argument("--file", default="CITATION.cff")
a = p.parse_args()
path = pathlib.Path(a.file)
txt = path.read_text(encoding="utf-8")
txt = re.sub(r'(?m)^version:\s*".*"$', f'version: "{a.version}"', txt)
if a.doi:
    if re.search(r'(?m)^doi:\s*".*"$', txt):
        txt = re.sub(r'(?m)^doi:\s*".*"$', f'doi: "{a.doi}"', txt)
    else:
        txt = txt.rstrip() + f'\ndoi: "{a.doi}"\n'
path.write_text(txt, encoding="utf-8")
print(f"[ok] {a.file} updated → version={a.version}, doi={a.doi or 'unchanged'}")
