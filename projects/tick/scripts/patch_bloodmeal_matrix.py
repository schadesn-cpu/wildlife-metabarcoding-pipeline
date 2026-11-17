import sys
p='scripts/tick_bloodmeal_matrix.py'
s=open(p).read()
n0=len(s); log=[]
def sub(old,new,label):
    global s
    assert s.count(old)==1, f'{label}: found {s.count(old)}'
    s=s.replace(old,new); log.append(label)

# 1. Ked is NOT a blood meal. Pull it out of the matrix and report it separately.
sub('''# which cytb bloodmeal columns collapse into "Cervidae"
CERVID_COLS = {"cytb_Deer", "cytb_Moose", "cytb_Cervidae_unresolved"}
# display order / labels for bloodmeal columns after collapsing
BM_ORDER = ["Cervidae", "Deer_ked", "Bird", "Bear", "Opossum", "Raccoon", "Porcupine", "Rabbit"]
BM_LABEL = {"Cervidae":"Cervidae\\n(deer/moose)", "Deer_ked":"Deer ked", "Bird":"Bird",
            "Bear":"Bear", "Opossum":"Opossum", "Raccoon":"Raccoon",
            "Porcupine":"Porcupine", "Rabbit":"Rabbit"}''',
'''# which cytb bloodmeal columns collapse into "Cervidae"
CERVID_COLS = {"cytb_Deer", "cytb_Moose", "cytb_Cervidae_unresolved"}

# I removed Deer_ked from the blood-meal columns on 2026-08-11. A deer ked is an
# ectoparasite whose DNA co-amplifies with the tick's blood meal, not a host the
# tick fed on, and this figure reads off-diagonal cells as prior-host meals. It is
# reported as a separate co-detection line beneath the matrix instead. The ked is
# confirmed on three markers (18S 100%, cytb accession-verified, 16S endosymbionts)
# and is a real finding; it just is not a blood meal.
NOT_A_BLOODMEAL = {"Deer_ked"}

# display order / labels for bloodmeal columns after collapsing
BM_ORDER = ["Cervidae", "Bird", "Bear", "Opossum", "Raccoon", "Porcupine", "Rabbit"]
BM_LABEL = {"Cervidae":"Cervidae\\n(deer/moose)", "Bird":"Bird",
            "Bear":"Bear", "Opossum":"Opossum", "Raccoon":"Raccoon",
            "Porcupine":"Porcupine", "Rabbit":"Rabbit"}''','ked out of BM_ORDER')

# 2. keep the ked counts so they can be reported
sub('''    for c in d.columns:
        if c.startswith("cytb_") and c not in CERVID_COLS:
            name=c.replace("cytb_","")
            bm[name]=(d[c]==1).astype(int)
    return bm''',
'''    for c in d.columns:
        if c.startswith("cytb_") and c not in CERVID_COLS:
            name=c.replace("cytb_","")
            bm[name]=(d[c]==1).astype(int)
    return bm

def ked_by_host(d, host_col, bm):
    """Ked co-detections per collection host. Reported under the matrix, not in it."""
    if "Deer_ked" not in bm.columns:
        return {}
    out={}
    for host in d[host_col].dropna().unique():
        idx=d.index[d[host_col]==host]
        n_ked=int(bm.loc[idx,"Deer_ked"].sum())
        if n_ked:
            out[host]=(n_ked, len(idx))
    return out''','ked_by_host helper')

# 3. drop the ked column before the matrix is built
sub('''    bm_cols=[c for c in BM_ORDER if c in bm.columns]''',
'''    bm_cols=[c for c in BM_ORDER if c in bm.columns and c not in NOT_A_BLOODMEAL]''','exclude ked from cols')

# 4. footnote
sub('''    fig.text(0.01,-0.03,"Cervidae merges Deer/Moose/Cervidae (cytb can't resolve cervid species). "
             "Diagonal = meal matches collection host; off-diagonal = prior-host meals.",
             fontsize=7,style='italic',color='#666')''',
'''    note=("Cervidae merges Deer/Moose/Cervidae (cytb can't resolve cervid species). "
          "Diagonal = meal matches collection host; off-diagonal = prior-host meals.")
    if KED_NOTE:
        note += "\\n" + KED_NOTE
    fig.text(0.01,-0.03,note,fontsize=7,style='italic',color='#666')''','footnote with ked')

# 5. main: compute and print the ked line
sub('''    d=pd.read_csv(a.input,sep='\\t')
    bm=collapse_bloodmeal(d)''',
'''    d=pd.read_csv(a.input,sep='\\t')
    bm=collapse_bloodmeal(d)

    global KED_NOTE
    keds=ked_by_host(d,a.host_col,bm)
    if keds:
        parts=[f"{h} {k}/{n}" for h,(k,n) in sorted(keds.items(), key=lambda kv:-kv[1][0])]
        total=sum(k for k,_ in keds.values())
        KED_NOTE=("Deer ked (Lipoptena cervi) DNA co-detected in %d tick(s): %s. "
                  "A ked is an ectoparasite sharing the host, NOT a blood meal, so it is "
                  "excluded from the matrix." % (total, "; ".join(parts)))
        print("\\nKED CO-DETECTION (not a blood meal): " + KED_NOTE)
    else:
        KED_NOTE=""''','main computes ked note')

sub("""CERVID_COLS = {"cytb_Deer", "cytb_Moose", "cytb_Cervidae_unresolved"}""",
"""CERVID_COLS = {"cytb_Deer", "cytb_Moose", "cytb_Cervidae_unresolved"}
KED_NOTE = ""   # filled in main(); appended to the figure footnote""",'KED_NOTE global')

open(p,'w').write(s)
print(f'{n0} -> {len(s)}')
for l in log: print('  +',l)
