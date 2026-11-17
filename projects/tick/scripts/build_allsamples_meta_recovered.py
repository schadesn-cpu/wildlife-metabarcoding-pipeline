#!/usr/bin/env python3
"""
Build all-samples 16S metadata with host RECOVERED from names + cytb bloodmeal.
Run on ron:  python scripts/build_allsamples_meta_recovered.py <output_dir>
Writes <output_dir>/meta_all_recovered.tsv for use with 09_plot_diversity.py.

Recovery logic (host_source column records provenance):
  1. collection : sheet__Host Spp is populated (real collection record)
  2. name       : sample-id contains 'Deer'->WTD or 'Bear'->Bear (record in name)
  3. bloodmeal  : cytb dominant_host_family, ONLY if dominant_pct >= 75 and a
                  clean family (Cervidae->WTD*, Mustelidae->Mustelidae,
                  Ursidae->Bear*).  *inferred, flagged in host_source.
  4. unknown    : nothing resolved (bare c__Mammalia, no record, etc.)
Bloodmeal-derived hosts are labelled inferred in host_source so the paper can
keep them distinct from collection records (a tick's bloodmeal is what it fed
on, not necessarily where it was collected).
"""
import sys, re
import pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
META = "reads_staged/16S_metadata.tsv"
BM   = "results/cytb/bloodmeal_per_sample_flagged.tsv"

FAMILY_TO_HOST = {   # cytb family -> figure host label (starred = inferred)
    "Cervidae":    "WTD",
    "Ursidae":     "Bear",
    "Mustelidae":  "Mustelidae",
}
BM_MIN_PCT = 75.0    # only trust bloodmeal host above this dominance

m = pd.read_csv(META, sep="\t", dtype=str)
m = m[m['is_control'].astype(str).str.lower().isin(['no','false','nan',''])].copy()

bm = pd.read_csv(BM, sep="\t")
def norm(s): return re.sub(r'-16s.*|-18s.*|-ITS.*|\..*','',str(s),flags=re.I).strip().lower()
bm['key'] = bm['sample_id'].map(norm)
bm_dom = bm.set_index('key')['dominant_host_family'].to_dict()
bm_pct = bm.set_index('key')['dominant_pct_of_host'].to_dict()

def resolve(row):
    raw = row['sheet__Host Spp']
    sid = str(row['sample-id'])
    # 1. collection record
    if pd.notna(raw) and str(raw) not in ('0','not recorded',''):
        h = str(raw)
        if h == 'Deer': h = 'WTD'
        return h, 'collection'
    # 2. name-based
    if re.search(r'Deer', sid, re.I):  return 'WTD',  'name'
    if re.search(r'Bear', sid, re.I):  return 'Bear', 'name'
    # 3. bloodmeal
    k = norm(sid)
    fam = bm_dom.get(k); pct = bm_pct.get(k)
    if fam in FAMILY_TO_HOST and pd.notna(pct) and float(pct) >= BM_MIN_PCT:
        return FAMILY_TO_HOST[fam], 'bloodmeal_inferred'
    return 'unknown', 'unresolved'

res = m.apply(resolve, axis=1, result_type='expand')
m['host'] = res[0]; m['host_source'] = res[1]

def scrub(x):
    if x is None: return 'unknown'
    s = str(x)
    if s.strip()=='' or s.lower()=='nan': return 'unknown'
    return re.sub(r'\s+',' ',s).strip()

out = pd.DataFrame({
    'sample-id': m['sample-id'].map(scrub),
    'host': m['host'].map(scrub),
    'host_source': m['host_source'].map(scrub),
    'batch': m['batch_short'].map(scrub),
    'state': m['sheet__State'].fillna('unknown').map(scrub),
})
path = f"{OUT}/meta_all_recovered.tsv"
with open(path,'w') as f:
    f.write("sample-id\thost\thost_source\tbatch\tstate\n")
    f.write("#q2:types\tcategorical\tcategorical\tcategorical\tcategorical\n")
    for _,r in out.iterrows():
        f.write("\t".join([r['sample-id'],r['host'],r['host_source'],r['batch'],r['state']])+"\n")

print("wrote", path, "rows:", len(out))
print("\nhost after recovery:")
print(out['host'].value_counts().to_string())
print("\nhost_source (provenance):")
print(out['host_source'].value_counts().to_string())
