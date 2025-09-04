#!/usr/bin/env python3
import requests, time, sys, csv
from pathlib import Path

OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

def fetch_associations(ensembl_ids, disease_efo="EFO_0009671", sleep=0.2):
    q = """
    query($ensemblId: String!, $efoId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        associatedDiseases(efoId: $efoId) {
          rows { disease { id name } score }
        }
        tractability {
          smallmolecule { topCategory }
          antibody { topCategory }
        }
      }
    }
    """
    rows = []
    for gid in ensembl_ids:
        r = requests.post(OT_URL, json={"query": q, "variables": {"ensemblId": gid, "efoId": disease_efo}})
        if not r.ok:
            print(f"[warn] Open Targets request failed for {gid}: {r.status_code}", file=sys.stderr)
            continue
        data = r.json()
        t = data.get("data", {}).get("target", {})
        if not t:
            continue
        score = None
        assoc = t.get("associatedDiseases", {}).get("rows", [])
        if assoc:
            score = assoc[0].get("score")
        sm = (t.get("tractability", {}).get("smallmolecule", {}) or {}).get("topCategory")
        ab = (t.get("tractability", {}).get("antibody", {}) or {}).get("topCategory")
        rows.append([t.get("id"), t.get("approvedSymbol"), score, sm, ab])
        time.sleep(sleep)
    return rows

def main():
    ensembl_ids = ["ENSG00000141510"]  # Example: TP53 placeholder; replace with your DE genes
    rows = fetch_associations(ensembl_ids)
    out = OUT / "open_targets_example.csv"
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ensembl_id", "symbol", "score_systemic_sclerosis", "tractability_small_molecule", "tractability_antibody"])
        w.writerows(rows)
    print(f"[done] wrote {out}")

if __name__ == "__main__":
    main()
