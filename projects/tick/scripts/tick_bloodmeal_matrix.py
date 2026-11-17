#!/usr/bin/env python3
"""
tick_bloodmeal_matrix.py
Build a collection-host x blood-meal matrix (heatmap) for ticks, from the locked
single-tick file. Structured to ALSO support a 3-way (tick species x host x blood meal)
once a tick-species column is joined in (pending crosswalk / genome-center tick ID).

USAGE (single-tick, works now):
  python tick_bloodmeal_matrix.py \
      --input results/single_tick_v3/single_tick_v3_full.tsv \
      --out results/figures/bloodmeal_matrix

OPTIONAL 3-way (when you have tick species per sample):
  --tick-species results/18S/all/tick_genus_per_sample.tsv \
  --tick-join-col sample_id --tick-value-col dominant
  (requires the input to have a matching key; see JOIN note below)

DESIGN / HONESTY NOTES (baked in):
  - Cervidae columns (Deer, Moose, Cervidae_unresolved) are COLLAPSED into one "Cervidae"
    column, because cytb cannot resolve deer vs moose. Showing them separately would falsely
    imply a single tick fed on multiple cervids.
  - Cell = tick count; color = % of that collection host's ticks (row-normalized).
  - Strong diagonal = blood meal matches collection host (validates the cytb pipeline).
  - Off-diagonal = prior-host meals (three-host tick biology) or questing (no meal).
  - For SINGLE ticks, tick species is ~uniformly Ixodes; the 3-way is only informative on the
    FULL dataset (where Dermacentor/winter-tick appears). Script warns if species is constant.
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# which cytb bloodmeal columns collapse into "Cervidae"
CERVID_COLS = {"cytb_Deer", "cytb_Moose", "cytb_Cervidae_unresolved"}
# display order / labels for bloodmeal columns after collapsing
# I removed Deer_ked from the blood-meal columns on 2026-08-11. A deer ked is an
# ectoparasite whose DNA co-amplifies with the tick's blood meal, not a host the
# tick fed on, and this figure reads off-diagonal cells as prior-host meals. It is
# reported as a co-detection line beneath the matrix instead. The ked is confirmed
# on three markers and is a real finding; it just is not a blood meal.
NOT_A_BLOODMEAL = {"Deer_ked"}
# 'WTD' and 'Deer' are the same animal recorded two ways. Grouping on the raw
# string splits white-tailed deer into two rows, which broke this figure on
# 2026-08-11 (74 WTD + 4 Deer instead of one row of 78). build_single_tick_set.py
# carries host_bin for this reason, but host_bin collapses all six wildlife hosts
# into 'Other wildlife' and this figure needs the species detail off-diagonal.
# So I normalise the labels here and keep the species-level column.
HOST_LABEL_ALIASES = {
    "wtd": "White-tailed deer",
    "deer": "White-tailed deer",
    "white-tailed deer": "White-tailed deer",
    "odocoileus virginianus": "White-tailed deer",
    # Abbreviation expansions only. These add no information the collection
    # record does not already give; they make this figure agree with
    # fig_tick_species_by_host, which spells the names out.
    "na porcupine": "North American porcupine",
    "n.a. porcupine": "North American porcupine",
    "porcupine": "North American porcupine",
    "virginia opossum": "Virginia opossum",
    "eastern gray squirrel": "Eastern gray squirrel",
    "moose": "Moose",
    # NOT folded: 'Bear'. The collection record says Bear;
    # fig_tick_species_by_host says Black bear, which asserts a species the
    # record does not give. Resolve that deliberately rather than here.
}
KED_NOTE = ""   # filled in main(); appended to the figure footnote
# Human is deliberately NOT a blood-meal column. Human cytb is in 258 of ~830
# libraries, the highest count in the entire run is a POSITIVE CONTROL (124,701
# reads) and an extraction blank carries 8,938, and counts run continuously from
# 1 to 124,701 with no separation between positive and negative samples. That is
# reagent and handling contamination tracking library depth, not feeding. It is
# reported in the footnote instead, where it can be stated accurately.
BM_ORDER = ["Cervidae", "Bird", "Bear", "Opossum", "Raccoon", "Porcupine", "Rabbit"]
BM_LABEL = {"Cervidae":"Cervidae\n(deer/moose)", "Bird":"Bird",
            "Bear":"Bear", "Opossum":"Opossum", "Raccoon":"Raccoon",
            "Porcupine":"Porcupine", "Rabbit":"Rabbit"}

def collapse_bloodmeal(d):
    """return DataFrame with one 0/1 col per collapsed bloodmeal category, per tick."""
    bm=pd.DataFrame(index=d.index)
    # Cervidae = any of the cervid cols
    cerv=[c for c in CERVID_COLS if c in d.columns]
    bm["Cervidae"]=(d[cerv].sum(axis=1)>0).astype(int) if cerv else 0
    for c in d.columns:
        if c.startswith("cytb_") and c not in CERVID_COLS:
            name=c.replace("cytb_","")
            bm[name]=(d[c]==1).astype(int)
    return bm

def build_matrix(d, host_col, bm, species=None):
    """
    Returns (row_labels, col_labels, count_matrix).
    If species is a Series aligned to d, rows become (species, host) pairs (3-way);
    otherwise rows are collection hosts (2-way).
    """
    bm_cols=[c for c in BM_ORDER if c in bm.columns and c not in NOT_A_BLOODMEAL]
    # add a 'No meal' column
    has_meal=bm[bm_cols].sum(axis=1)>0
    if species is not None:
        d=d.assign(_sp=species.values)
        groups=d.groupby(["_sp", host_col])
        rows=[]; labels=[]
        for (sp,h),idx in groups.groups.items():
            sub=bm.loc[idx]
            counts=[int(sub[c].sum()) for c in bm_cols]
            nomeal=int((sub[bm_cols].sum(axis=1)==0).sum())
            rows.append(counts+[nomeal]); labels.append(f"{sp}\n{h} (n={len(idx)})")
        return labels, [BM_LABEL.get(c,c) for c in bm_cols]+["No meal\ndetected"], np.array(rows,float)
    else:
        hosts=d[host_col].value_counts().index.tolist()
        rows=[]; labels=[]
        for h in hosts:
            idx=d.index[d[host_col]==h]
            sub=bm.loc[idx]
            counts=[int(sub[c].sum()) for c in bm_cols]
            nomeal=int((sub[bm_cols].sum(axis=1)==0).sum())
            rows.append(counts+[nomeal]); labels.append(f"{h} (n={len(idx)})")
        return labels, [BM_LABEL.get(c,c) for c in bm_cols]+["No meal\ndetected"], np.array(rows,float)

def plot_matrix(rowlab, collab, M, out, title):
    from matplotlib.colors import LinearSegmentedColormap
    # Sequential colormap built from the Wong/Okabe-Ito blue (#0072B2): white -> Wong blue.
    # (A heatmap needs a sequential gradient; the Wong palette is categorical, so we derive a
    #  gradient anchored on the project's blue to stay palette-consistent + colorblind-safe.)
    wong_blues = LinearSegmentedColormap.from_list("wong_blues", ["#FFFFFF", "#56B4E9", "#0072B2"])
    # Colour is the percentage of THAT HOST'S specimens, not of the row. A
    # specimen can carry two meals, so a row can exceed 100% and normalising by
    # the row total would silently rescale that away.
    import re as _re
    _n=[]
    for _lab in rowlab:
        _m=_re.search(r"n=(\d+)", _lab)
        _n.append(int(_m.group(1)) if _m else 1)
    n_spec=np.array(_n, dtype=float).reshape(-1,1)
    n_spec[n_spec==0]=1
    pct=100*M/n_spec
    fig,ax=plt.subplots(figsize=(max(8,0.9*len(collab)+3), max(4.5,0.6*len(rowlab)+2)))
    im=ax.imshow(pct,cmap=wong_blues,vmin=0,vmax=100,aspect="auto")
    ax.set_xticks(range(len(collab))); ax.set_xticklabels(collab,fontsize=8.5)
    ax.set_yticks(range(len(rowlab))); ax.set_yticklabels(rowlab,fontsize=8.5)
    ax.set_xlabel("Blood-meal source (cytb)",fontsize=10.5)
    ax.set_ylabel("Collection host",fontsize=10.5)
    ax.set_title(title,fontsize=11,fontweight='bold',loc='left')
    for i in range(len(rowlab)):
        for j in range(len(collab)):
            c=int(M[i,j])
            if c>0:
                ax.text(j,i,str(c),ha='center',va='center',fontsize=9,
                        color='white' if pct[i,j]>55 else '#222',fontweight='bold')
    cbar=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.02); cbar.set_label("% of that host's specimens",fontsize=8.5)
    note=("Cervidae merges Deer/Moose/Cervidae (cytb can't resolve cervid species). "
          "Diagonal = meal matches collection host; off-diagonal = prior-host meals. "
          "A tick can return more than one blood meal, so the counts in a row can "
          "sum to more than that row's specimen count; colour is the percentage of "
          "that host's specimens, not of the row.")
    if KED_NOTE: note += "\n" + KED_NOTE
    note += ("\nHuman cytb was detected in 258 of ~830 libraries and is NOT shown as a "
             "blood meal. The highest human read count in the run is a positive control "
             "(124,701 reads) and an extraction blank carries 8,938; counts run "
             "continuously from 1 to 124,701 with no separation between positive and "
             "negative samples, consistent with reagent and handling contamination rather "
             "than feeding. 13 single ticks carry human reads, median 529.")
    # Wrap the footnote. Unwrapped it is one line the full width of the figure,
    # which forced a ~3.5:1 aspect ratio and would be unreadable at column width.
    import textwrap as _tw
    wrapped = "\n".join("\n".join(_tw.wrap(line, 150)) if line.strip() else line
                        for line in note.split("\n"))
    fig.text(0.01,-0.03,wrapped,fontsize=7,style='italic',color='#666',
             va='top')
    fig.tight_layout()
    for e in ('png','pdf','svg'):
        fig.savefig(f'{out}.{e}',dpi=300,bbox_inches='tight',facecolor='white')
    print("wrote", out+".{png,pdf,svg}")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--host-col',default='Host Spp')
    ap.add_argument('--out',default='results/figures/bloodmeal_matrix')
    ap.add_argument('--tick-species',default=None,help='optional tsv with per-sample tick species')
    ap.add_argument('--tick-join-col',default='sample_id')
    ap.add_argument('--tick-value-col',default='dominant')
    ap.add_argument('--input-join-col',default='Barcode',help='column in --input to join tick species on')
    a=ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or '.',exist_ok=True)

    d=pd.read_csv(a.input,sep='\t')
    bm=collapse_bloodmeal(d)

    # Normalise host labels before anything counts or groups on them.
    _raw = d[a.host_col].astype(str).str.strip()
    d[a.host_col] = _raw.where(~_raw.str.lower().isin(HOST_LABEL_ALIASES),
                               _raw.str.lower().map(HOST_LABEL_ALIASES))
    d.loc[_raw.str.lower().isin(["nan", "", "none"]), a.host_col] = None
    _changed = {}
    for _before, _after in zip(_raw, d[a.host_col]):
        if _after is not None and str(_before) != str(_after):
            _changed[str(_before)] = str(_after)
    if _changed:
        print("normalised host labels:")
        for _before in sorted(_changed):
            _n = int((_raw == _before).sum())
            print("    %-28s -> %-28s (%d specimen(s))"
                  % (_before, _changed[_before], _n))
    # Any label NOT normalised is printed too. A new spelling should be visible
    # here on the run that introduces it, not as a surprise extra row in the
    # figure six weeks later.
    _untouched = sorted(set(str(v) for v in _raw
                            if str(v).lower() not in HOST_LABEL_ALIASES
                            and str(v).lower() not in ("nan", "", "none")))
    if _untouched:
        print("host labels used as-is (add to HOST_LABEL_ALIASES if any of "
              "these are the same animal under two names):")
        for _u in _untouched:
            print("    %-28s (%d specimen(s))" % (_u, int((_raw == _u).sum())))

    global KED_NOTE
    if "Deer_ked" in bm.columns:
        per_host={}
        for host in d[a.host_col].dropna().unique():
            idx=d.index[d[a.host_col]==host]
            n_ked=int(bm.loc[idx,"Deer_ked"].sum())
            if n_ked: per_host[host]=(n_ked,len(idx))
        if per_host:
            parts=[f"{h} {k}/{n}" for h,(k,n) in sorted(per_host.items(),key=lambda kv:-kv[1][0])]
            total=sum(k for k,_ in per_host.values())
            KED_NOTE=("Deer ked (Lipoptena cervi) DNA co-detected in %d tick(s): %s. A ked is an "
                      "ectoparasite sharing the host, NOT a blood meal, so it is excluded from the "
                      "matrix." % (total, "; ".join(parts)))
            KED_NOTE += (" This is NOT a detection rate: keds were removed from tubes "
                         "after the first processing batch (Kellen Sawyer, 2026-09-14), and ked "
                         "co-detection ranges from 32.7% to 0% across sequencing batches, so it "
                         "reflects handling rather than host ecology.")
            print("\nKED CO-DETECTION (not a blood meal): "+KED_NOTE)

    # The matrix rows come from value_counts() on the host column, which drops
    # NaN. Titling with len(d) therefore claimed n=93 while plotting 89. I title
    # with the number actually drawn and name the dropped specimens in the note.
    n_plotted = int(d[a.host_col].notna().sum())
    n_nohost = len(d) - n_plotted
    if n_nohost:
        KED_NOTE += (" %d specimen(s) with no recorded collection host are not shown; "
                     "the matrix plots %d of %d." % (n_nohost, n_plotted, len(d)))
        # Two figures in this manuscript carry n=89 and they are NOT the same
        # ticks. This one is the sequenced set minus the specimens with no
        # recorded collection host. The qPCR concordance figure is the same
        # sequenced set minus the specimens whose only qPCR record is a pool.
        # Zero overlap between the two excluded sets. Say which, or a reader
        # reasonably assumes one cohort.
        KED_NOTE += (" This is not the same subset as the qPCR concordance figure, "
                     "which also plots %d: that one excludes specimens whose only "
                     "qPCR record is a pooled subsample, while this one excludes "
                     "specimens with no recorded collection host. No specimen is "
                     "in both excluded groups." % n_plotted)
    species=None; title="Blood meal vs collection host (n=%d)"%n_plotted
    if a.tick_species and os.path.exists(a.tick_species):
        ts=pd.read_csv(a.tick_species,sep='\t')
        # NOTE: join keys often differ (run-suffixes/replicates) — see Tick_Crosswalk_Task.md.
        # naive join; strip run-suffix on the species side.
        ts['_key']=ts[a.tick_join_col].astype(str).str.split('.').str[0]
        keymap=dict(zip(ts['_key'], ts[a.tick_value_col]))
        species=d[a.input_join_col].astype(str).map(keymap)
        nmatched=species.notna().sum()
        print("tick-species matched: %d / %d" % (nmatched, len(d)))
        if nmatched < 0.5*len(d):
            print("** WARNING: <50%% matched — join key mismatch (crosswalk needed). Falling back to 2-way. **")
            species=None
        elif species.dropna().nunique()<=1:
            print("** NOTE: tick species is constant (%s) — 3-way uninformative. Using 2-way. **"
                  % species.dropna().unique())
            species=None
        else:
            title="Tick species x collection host vs blood meal (n=%d)"%n_plotted

    rowlab,collab,M=build_matrix(d,a.host_col,bm,species)
    # print the matrix as text too
    print("\n"+"\t".join(["ROW"]+[c.replace('\n',' ') for c in collab]))
    for lab,row in zip(rowlab,M):
        print("\t".join([lab.replace('\n',' ')]+[str(int(x)) for x in row]))
    plot_matrix(rowlab,collab,M,a.out,title)

if __name__=='__main__':
    main()
