import os, json, time, urllib.parse, urllib.request, sys, re, subprocess, pathlib

ZENODO_BASE = os.environ.get("ZENODO_BASE","https://zenodo.org").rstrip("/")
OWNER = os.environ.get("OWNER") or (os.environ.get("GITHUB_REPOSITORY","/").split("/")[0])
REPO  = os.environ.get("REPO")  or (os.environ.get("GITHUB_REPOSITORY","/").split("/")[1] if "/" in os.environ.get("GITHUB_REPOSITORY","/") else "")
TAG   = os.environ.get("TAG") or os.environ.get("GITHUB_REF_NAME","")

if not TAG:
    print("::warning ::No TAG env; exiting neutral.")
    sys.exit(78)

release_url = f"https://github.com/{OWNER}/{REPO}/releases/tag/{TAG}"
q = f'metadata.related_identifiers.identifier:"{release_url}"'
url = f'{ZENODO_BASE}/api/records/?q={urllib.parse.quote(q)}&size=5&sort=mostrecent'

def fetch(u):
    with urllib.request.urlopen(u) as r:
        return json.load(r)

hits = []
for i in range(20):  # up to ~5 min
    try:
        data = fetch(url)
        hits = data.get("hits", {}).get("hits", [])
        if hits:
            break
    except Exception as e:
        print("Fetch error:", e)
    print(f"Zenodo record not found yet (attempt {i+1}/20); sleeping 15s...")
    time.sleep(15)

if not hits:
    print("::warning ::No Zenodo record found; exiting neutral.")
    sys.exit(78)

rec = hits[0]
doi = rec.get("doi")
conceptdoi = rec.get("conceptdoi") or rec.get("metadata", {}).get("conceptdoi")
rels = rec.get("metadata", {}).get("related_identifiers", []) or []
if not any(release_url == (ri.get("identifier") or "") for ri in rels):
    print("::error ::Found a record but it doesn't reference this release URL; aborting.")
    sys.exit(1)

print(f"Found DOI: {doi}")
print(f"Concept DOI: {conceptdoi}")

def patch_file(path, pattern, repl, flags=re.M):
    p = pathlib.Path(path)
    if not p.exists(): return False
    txt = p.read_text(encoding='utf-8')
    new = re.sub(pattern, repl, txt, flags=flags)
    if new != txt:
        p.write_text(new, encoding='utf-8')
        return True
    return False

changed = False

# CITATION.cff
p_cff = pathlib.Path("CITATION.cff")
if p_cff.exists():
    txt = p_cff.read_text(encoding='utf-8')
    if "identifiers:" in txt:
        if re.search(r'(?ms)-\s*type:\s*doi\s*\n\s*value:\s*".*?"', txt):
            changed |= patch_file("CITATION.cff",
                                  r'(?ms)(-\s*type:\s*doi\s*\n\s*value:\s*).*',
                                  r'\1"'+doi+'"')
        else:
            changed |= patch_file("CITATION.cff",
                                  r'(?ms)(identifiers:\s*(?:-.*\n)*)$',
                                  r'\1- type: doi\n  value: "'+doi+'"\n')
    else:
        p_cff.write_text(txt + f'\nidentifiers:\n- type: doi\n  value: "{doi}"\n', encoding='utf-8')
        changed = True

# README.md
p_md = pathlib.Path("README.md")
if p_md.exists():
    txt = p_md.read_text(encoding='utf-8')
    if "DOI:" in txt:
        changed |= patch_file("README.md", r'(DOI:\s*)(\S+)', r'\1'+doi)
    else:
        p_md.write_text(txt + f'\n\n**DOI:** {doi}\n', encoding='utf-8')
        changed = True

if not changed:
    print("::notice ::No file changes needed.")
    sys.exit(0)

subprocess.check_call(["git","config","user.name","doi-bot"])
subprocess.check_call(["git","config","user.email","doi-bot@users.noreply.github.com"])
branch = f"update-doi-{TAG}"
subprocess.check_call(["git","checkout","-b", branch])
subprocess.check_call(["git","add","CITATION.cff","README.md"])
subprocess.check_call(["git","commit","-m", f"Add Zenodo DOI for {TAG}: {doi}"])
subprocess.check_call(["git","push","origin","HEAD"])
print("::notice ::Pushed branch with DOI update.")
