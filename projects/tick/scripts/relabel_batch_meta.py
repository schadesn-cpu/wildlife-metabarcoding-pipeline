#!/usr/bin/env python3
"""
Add clean, human-readable batch and host labels to a diversity metadata TSV,
so 09_plot_diversity.py legends read nicely for slides/publication.
Usage: python scripts/relabel_batch_meta.py <meta.tsv>
Writes <meta>_clean.tsv with added 'batch_clean' column.
Maps the long run IDs -> Batch 1..N ordered by date; keeps a lookup printed.
"""
import sys, re
import pandas as pd

path = sys.argv[1]
# read, preserving the #q2:types row
with open(path) as f:
    lines = f.readlines()
header = lines[0].rstrip("\n").split("\t")
types  = lines[1].rstrip("\n").split("\t") if lines[1].startswith("#q2:types") else None
data_start = 2 if types else 1
rows = [l.rstrip("\n").split("\t") for l in lines[data_start:] if l.strip()]
df = pd.DataFrame(rows, columns=header)

# batch column present?
bcol = 'batch' if 'batch' in df.columns else None
if bcol is None:
    print("no 'batch' column; nothing to do"); sys.exit(0)

# order batches by the leading date (YYMMDD) in the run id, map to Batch N
def datekey(b):
    m = re.match(r'(\d{6})', str(b))
    return m.group(1) if m else '999999'
uniq = sorted(df[bcol].unique(), key=datekey)
# human label: 'Batch NN (mon-yy)' from the date
def human(b, i):
    m = re.match(r'(\d{2})(\d{2})(\d{2})', str(b))
    if m:
        yy, mm, dd = m.groups()
        mon = ['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][int(mm)]
        return f"Batch {i} ({mon} '{yy})"
    return f"Batch {i}"
mapping = {b: human(b, i+1) for i, b in enumerate(uniq)}
df['batch_clean'] = df[bcol].map(mapping)

# host clean (title-case, already fine but ensure consistency)
out = path.replace('.tsv', '_clean.tsv')
cols = list(df.columns)
with open(out, 'w') as f:
    f.write("\t".join(cols) + "\n")
    f.write("#q2:types\t" + "\t".join(["categorical"]*(len(cols)-1)) + "\n")
    for _, r in df.iterrows():
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

print("wrote", out)
print("\nbatch label mapping (date-ordered):")
for b, h in mapping.items():
    print(f"  {h:22s} <- {b}")
