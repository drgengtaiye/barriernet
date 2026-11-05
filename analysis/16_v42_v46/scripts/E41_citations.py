"""Full citations of the 19 trial papers in results/E41_adjudicated.csv, from the PubMed esummary record saved in
results/E41_pubmed_esummary.json (NCBI E-utilities, fetched 2026-11-01). Writes results/E41_citations.csv."""
import json
from pathlib import Path
import pandas as pd

R = Path(__file__).resolve().parent / "results"
J = json.load(open(R / "E41_pubmed_esummary.json"))["result"]
A = pd.read_csv(R / "E41_adjudicated.csv", dtype=str)
rows = []
for t, pmid in zip(A.trial, A.pmid.str.strip()):
    r = J[pmid]; au = [a["name"] for a in r["authors"] if a.get("authtype") == "Author"]
    who = ", ".join(au[:6]) + (", et al." if len(au) > 6 else ".")
    vol = r.get("volume", ""); iss = f"({r['issue']})" if r.get("issue") else ""; pg = f":{r['pages']}" if r.get("pages") else ""
    rows.append(dict(trial=t, pmid=pmid, citation=f"{who} {r['title'].rstrip('.')}. {r['source']}. {r['pubdate'][:4]}{';' + vol + iss if vol else ''}{pg}. "
                     f"PMID {pmid}."))
pd.DataFrame(rows).to_csv(R / "E41_citations.csv", index=False)
print(len(rows), "citations")
