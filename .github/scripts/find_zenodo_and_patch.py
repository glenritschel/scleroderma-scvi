#!/usr/bin/env python3
import os, re, sys, time, json, urllib.parse, urllib.request, subprocess, pathlib

def zenodo_get(url):
    with urllib.request.urlopen(url) as r:
        return json.load(r)

def search_by_related_url(release_url):
    q = f'metadata.related_identifiers.identifier:"{release_url}"'
    u = f'https://zenodo.org/api/records/?q={urllib.parse.quote(q)}&size=3&sort=mostrecent'
    return zenodo_get(u).get("hits", {}).get("hits", [])

def search_by_version(tag):
    # try with and without leading 'v'
    versions = {tag, tag.lstrip('v')}
    hits = []
    for v in versions:
        if not v: 
            continue
        q = f'metadata.version:"{v}"'
        u = f'https://zenodo.org/api/records/?q={urllib.parse.quote(q)}&size=3&sort=mostrecent'
        hits.extend(zenodo_get(u).get("hits", {}).get("hits", []))
    return hits

def patch(path, pattern, repl, flags=re.M):
    p = pathlib.Path(path)
    if not p.exists():
        return False
    txt = p.read_text(encoding="utf-8")
    new = re.sub(pattern, repl, txt, flags=flags)
    if new != txt:
        p.write_text(new, encoding="utf-8")
        return True
    return False

owner, repo = (os.environ.get("GITHUB_REPOSITORY") or "/").split("/", 1)
owner = owner or os.environ.get("OWNER", "")
repo  = repo  or os.environ.get("REPO", "")
tag   = os.environ.get("TAG") or os.environ.get("INPUT_TAG") or os.environ.get("GITHUB_REF_NAME")

if not owner or not repo:
    print("::error ::Could not determine owner/repo.")
    sys.exit(1)
if not tag:
    print("::error ::No tag provided or detected.")
    sys.exit(1)

release_url = f"https://github.com/{owner}/{repo}/releases/tag/{tag}"
print(f"Looking for Zenodo record referencing: {release_url}")

attempts = 40  # up to ~20 min with gentle backoff
for i in range(attempts):
    # Prefer exact relation by release URL
    hits = search_by_related_url(release_url)
    # Fallback search by version string
    if not hits:
        hits = search_by_version(tag)

    if hits:
        rec = hits[0]
        doi = rec.get("doi") or rec.get("metadata", {}).get("doi")
        conceptdoi = rec.get("conceptdoi") or rec.get("metadata", {}).get("conceptdoi")
        if not doi:
            print("::error ::Zenodo record found but DOI missing in payload.")
            sys.exit(1)

        print(f"Found DOI: {doi}")
        if conceptdoi:
            print(f"Concept DOI: {conceptdoi}")

        changed = False

        # CITATION.cff: add/replace a version DOI entry
        p_cff = pathlib.Path("CITATION.cff")
        if p_cff.exists():
            txt = p_cff.read_text(encoding="utf-8")
            if "identifiers:" in txt:
                if re.search(r'(?ms)-\s*type:\s*doi\s*\n\s*value:\s*".*?"', txt):
                    changed |= patch("CITATION.cff",
                                     r'(?ms)(-\s*type:\s*doi\s*\n\s*value:\s*).*',
                                     r'\1"'+doi+'"')
                else:
                    changed |= patch("CITATION.cff",
                                     r'(?ms)(identifiers:\s*(?:-.*\n)*)$',
                                     r'\1- type: doi\n  value: "'+doi+'"\n')
            else:
                p_cff.write_text(txt + f'\nidentifiers:\n- type: doi\n  value: "{doi}"\n', encoding="utf-8")
                changed = True

        # README.md: update or append a DOI line
        p_md = pathlib.Path("README.md")
        if p_md.exists():
            txt = p_md.read_text(encoding="utf-8")
            if "DOI:" in txt:
                changed |= patch("README.md", r'(DOI:\s*)(\S+)', r'\1'+doi)
            else:
                p_md.write_text(txt + f'\n\n**DOI:** {doi}\n', encoding="utf-8")
                changed = True

        if not changed:
            print("::notice ::No file changes needed.")
            sys.exit(0)

        subprocess.check_call(["git","config","user.name","doi-bot"])
        subprocess.check_call(["git","config","user.email","doi-bot@users.noreply.github.com"])
        branch = f"update-doi-{tag}"
        subprocess.check_call(["git","checkout","-b", branch])
        subprocess.check_call(["git","add","CITATION.cff","README.md"])
        subprocess.check_call(["git","commit","-m", f"Add Zenodo DOI for {tag}: {doi}"])
        subprocess.check_call(["git","push","origin","HEAD"])
        print("::notice ::Pushed branch with DOI update.")
        sys.exit(0)

    wait = 15 + min(i, 30)  # start 15s, creep to ~45s
    print(f"Zenodo record not found yet (attempt {i+1}/{attempts}); sleeping {wait}s...")
    time.sleep(wait)

print("::warning ::No Zenodo record found after extended wait; exiting neutral.")
sys.exit(78)
