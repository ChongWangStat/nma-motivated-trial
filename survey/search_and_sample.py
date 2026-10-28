"""Europe PMC searches and random samples for the survey (manuscript Section 2.1). Searches were run on 28 October 2026.
Counts change as Europe PMC is updated; the sampled identifiers used in the paper are in sample_ids.json."""
import json, random, urllib.request, urllib.parse, time
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
RCT = 'PUB_TYPE:"Randomized Controlled Trial" AND PUB_YEAR:2024'
MA = '(REF:"meta-analysis" OR REF:"meta-analyses" OR REF:"meta analysis" OR REF:"meta analyses")'
NMA = ('(REF:"network meta-analysis" OR REF:"network meta-analyses" OR REF:"network meta analysis" OR REF:"mixed treatment comparison" '
       'OR REF:"mixed treatment comparisons" OR REF:"multiple treatments meta-analysis" OR REF:"indirect treatment comparison")')
QUERIES = {  # counts reported in Section 3.1 use the open-access frame
    "oa_all": f"{RCT} AND OPEN_ACCESS:y",
    "oa_cite_ma": f"{RCT} AND OPEN_ACCESS:y AND IN_EPMC:y AND {MA}",
    "oa_cite_nma": f"{RCT} AND OPEN_ACCESS:y AND IN_EPMC:y AND {NMA}",
    "ft_cite_nma": f"{RCT} AND IN_EPMC:y AND {NMA}",                     # sampling frame, stratum 1
    "ft_cite_ma_only": f"{RCT} AND IN_EPMC:y AND {MA} AND NOT {NMA}",    # sampling frame, stratum 2
}
def get(params):
    for k in range(6):
        try: return json.load(urllib.request.urlopen(API + urllib.parse.urlencode(params), timeout=180))
        except Exception: time.sleep(10 * (k + 1))
def count(q): return get({"query": q, "format": "json", "pageSize": 1, "resultType": "idlist"})["hitCount"]
def ids(q):
    out, cur = [], "*"
    while True:
        d = get({"query": q, "format": "json", "pageSize": 1000, "resultType": "idlist", "cursorMark": cur})
        rs = d["resultList"]["result"]; out += [r["pmcid"] for r in rs if r.get("pmcid")]
        if not rs or d.get("nextCursorMark") == cur: return sorted(set(out))
        cur = d["nextCursorMark"]; time.sleep(1)
if __name__ == "__main__":
    for k, q in QUERIES.items(): print(k, count(q))
    random.seed(2024); nma = random.sample(ids(QUERIES["ft_cite_nma"]), 200)
    random.seed(2025); ma = random.sample(ids(QUERIES["ft_cite_ma_only"]), 100)
    print(len(nma), len(ma))  # full texts: https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML (open access only)
