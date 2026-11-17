#!/usr/bin/env python3
"""
build_single_tick_set.py  —  run ON RON, from the project root:
    cd ~/meedlab/tick_project
    python3 scripts/build_single_tick_set.py

WHY THIS EXISTS
---------------
`results/single_tick_v2/single_tick_v2_full.tsv` had NO build script. It was made
interactively, and it is the source of every n=103 number in the manuscript: the
single-tick prevalences, the 21.4% co-infection figure, fig2, the bloodmeal matrix,
and the metabarcoding side of the qPCR concordance. The figure regenerated from the
table; the table regenerated from nothing. This closes that gap.

WHAT IT DOES
------------
Applies the single-tick inclusion rule, which now lives in ONE place
(`is_single_tick` in crosswalk_4marker.tsv):

    A specimen is a single tick iff (a) its collection record gives a tick count of
    exactly 1, (b) no other record contradicts that, and (c) it is a wild-caught
    animal from the study region.

Excluded on that basis (2026-08-07): 8 Batch-12 fungal-treatment ticks (experimental,
not wild-caught) and 4 specimens whose extraction record contradicts the collection
record (USGSTICKS-6503/6508/6510/9216).

Reads:
    results/metadata/crosswalk_4marker.tsv          (the rule, as is_single_tick)
    results/single_tick_v2/single_tick_v2_full.tsv  (the per-specimen calls)
Writes:
    results/single_tick_v3/single_tick_v3_full.tsv  (the corrected set)
    results/single_tick_v3/single_tick_v3_report.txt
    results/figures/fig2_single_tick_prevalence_coinfection.png/.pdf

NOTE ON CO-INFECTION: counted at GENUS level. Babesia odocoilei + B. microti in one
tick is ONE genus, not two. This is what the published fig2 does, and it is why the
old set gave 22 ticks with >=2 genera rather than 23 with >=2 species.

Bartonella is NOT counted, as of 2026-09-18. See CO_DETECTIONS below for the
reasoning and the numbers. Counting it gave >=1 55/93 (59.1%) and >=2 19/93
(20.4%); excluding it gives >=1 51/93 (54.8%) and >=2 14/93 (15.1%). If you are
reading a document that quotes 59.1% or 20.4%, it predates this change.
"""
import os
import sys
import pandas as pd

CW   = "results/metadata/crosswalk_4marker.tsv"
# Regenerated into a temp file on every run rather than read from disk. The
# stored copy was dated 2026-06-13 and predated every correction since, so v3
# carried a fresh timestamp and June's calls.
V2_BUILDER = "scripts/build_single_tick_v2.py"
MIN_READS = 1   # unchanged from the v2 default; see the note printed at the end
ODIR = "results/single_tick_v3"
FDIR = "results/figures"

# pathogen columns in the v2 table -> the GENUS they belong to.
# Rickettsia and Francisella are endosymbionts and are NOT pathogens (locked rule).
PATHOGEN_GENUS = {
    '16S_Anaplasma':          'Anaplasma',
    '16S_Borreliella':        'Borreliella',
    '16S_Borrelia':           'Borrelia',
    '18S_Babesia_sp':         'Babesia',
    '18S_Babesia_odocoilei':  'Babesia',
    '18S_Babesia_microti':    'Babesia',
}
# Wolbachia is an endosymbiont here, not a pathogen: these ASVs are deer ked
# (Lipoptena cervi) symbionts that SILVA had filed under Anaplasma.
ENDOSYMBIONTS = ['16S_Rickettsia', '16S_Francisella', '16S_Wolbachia']

# Bartonella is a CO-DETECTION, not a pathogen, and is therefore not in
# PATHOGEN_GENUS above. Section 3.7 of the manuscript already says so: these
# are ruminant-associated Bartonella, not tick-borne human pathogens, reported
# outside the pathogen results. The co-infection count used to contradict that
# by counting them as a pathogen genus, which put 'this tick had a ked on it'
# into a co-infection statistic.
#
# Measured on these 93 single ticks, 2026-09-18: Bartonella 10/93, deer ked
# 11/93, and 6 of the 10 Bartonella-positive ticks carry ked reads (odds ratio
# 23.4, Fisher p = 0.0001). A seventh carries Wolbachia, the other ked
# endosymbiont. Keds were removed from tubes after the first processing batch
# (Kellen Sawyer, 2026-09-12), so the rate is partly a handling artifact.
#
# Reclassified, NOT deleted: the prevalence is still reported below.
CO_DETECTIONS = ['16S_Bartonella']

rep = []
def say(m=""):
    print(m); rep.append(str(m))

os.makedirs(ODIR, exist_ok=True)
os.makedirs(FDIR, exist_ok=True)

# ---- 1. the rule -----------------------------------------------------------
cw = pd.read_csv(CW, sep="\t")
if 'is_single_tick' not in cw.columns:
    sys.exit(f"ERROR: {CW} has no is_single_tick column. Re-run build_4marker_crosswalk.py.")
cw['is_single_tick'] = cw['is_single_tick'].astype(str).str.lower().isin(['true','1'])
keep = cw[cw['is_single_tick']].copy()
say(f"crosswalk: {len(cw)} specimens, {len(keep)} flagged is_single_tick")

# ---- 2. the old set --------------------------------------------------------
import subprocess, tempfile
_tmp = tempfile.NamedTemporaryFile(suffix="_v2.tsv", delete=False); _tmp.close()
say(f"regenerating pathogen calls from current taxonomy (min_reads={MIN_READS}) ...")
_r = subprocess.run([sys.executable, V2_BUILDER, "--min-reads", str(MIN_READS),
                     "--out", _tmp.name], capture_output=True, text=True)
if _r.returncode != 0:
    sys.stderr.write(_r.stdout + _r.stderr); sys.exit("v2 builder failed; nothing written")
old = pd.read_csv(_tmp.name, sep="\t")
say(f"regenerated calls: {len(old)} specimens")

old_ids  = set(old['Barcode'].astype(str).str.strip())
keep_ids = set(keep['sample_id'].astype(str).str.strip())

dropped = sorted(old_ids - keep_ids)
added   = sorted(keep_ids - old_ids)
shared  = old_ids & keep_ids

say(f"\n  in both:            {len(shared)}")
say(f"  dropped by rule:    {len(dropped)}")
for b in dropped: say(f"      - {b}")
say(f"  newly eligible:     {len(added)}")
for b in added: say(f"      + {b}")

# ---- 3. the corrected set --------------------------------------------------
new = old[old['Barcode'].astype(str).str.strip().isin(keep_ids)].copy()

missing = [b for b in added if b not in old_ids]
if missing:
    say(f"\n  WARNING: {len(missing)} specimen(s) meet the rule but have no row in the")
    say(f"  v2 table, so they have NO pathogen calls and cannot be included yet.")
    say(f"  They need their per-marker calls computed before the set is complete:")
    for b in missing: say(f"      ? {b}")

say(f"\ncorrected set: {len(new)} specimens")
if len(new) != len(keep):
    say(f"  (rule says {len(keep)}; the gap is the {len(missing)} without calls above)")

# refresh host/state from the corrected crosswalk rather than trusting the old file
meta = keep.set_index(keep['sample_id'].astype(str).str.strip())
# NOTE: do NOT map host_bin into 'Town'. An earlier version did, which overwrote
# the town with the host group. The crosswalk has no town column; leave Town alone.
for col_new, col_cw in [('Host Spp','host_spp'), ('State','state')]:
    if col_cw in meta.columns and col_new in new.columns:
        new[col_new] = new['Barcode'].astype(str).str.strip().map(meta[col_cw])

# I carry host_bin through as its own column. 'Host Spp' holds the raw collection
# string, so 'WTD' and 'Deer' are the same animal under two labels and split into
# two rows in any figure that groups on it (this broke the bloodmeal matrix on
# 2026-08-11: 74 WTD + 4 Deer instead of one row of 78). host_bin is the
# normalised group and is what figures should group on.
if 'host_bin' in meta.columns:
    new['host_bin'] = new['Barcode'].astype(str).str.strip().map(meta['host_bin'])
    n_no_host = int(new['host_bin'].isna().sum())
    if n_no_host:
        say(f"  note: {n_no_host} specimen(s) have no host_bin in the crosswalk")

out = os.path.join(ODIR, "single_tick_v3_full.tsv")
new.to_csv(out, sep="\t", index=False)
say(f"wrote {out}")

# ---- 4. prevalence + co-infection ------------------------------------------
n = len(new)
present = [c for c in PATHOGEN_GENUS if c in new.columns]

say(f"\n=== single-tick prevalence (n = {n}) ===")
prev = {}
for c in present:
    pos = int((pd.to_numeric(new[c], errors='coerce').fillna(0) > 0).sum())
    prev[c] = (pos, 100 * pos / n)
    say(f"  {c:<24} {pos:>3}/{n} = {100*pos/n:5.1f}%")

say(f"\n=== endosymbionts (context, NOT pathogens) ===")
for c in ENDOSYMBIONTS:
    if c in new.columns:
        pos = int((pd.to_numeric(new[c], errors='coerce').fillna(0) > 0).sum())
        say(f"  {c:<24} {pos:>3}/{n} = {100*pos/n:5.1f}%")

say(f"\n=== co-detections (context, NOT pathogens, NOT in co-infection) ===")
for c in CO_DETECTIONS:
    if c in new.columns:
        pos = int((pd.to_numeric(new[c], errors='coerce').fillna(0) > 0).sum())
        say(f"  {c:<24} {pos:>3}/{n} = {100*pos/n:5.1f}%")
        if '16S_Bartonella' == c and 'cytb_Deer_ked' in new.columns:
            bart = pd.to_numeric(new[c], errors='coerce').fillna(0) > 0
            ked = pd.to_numeric(new['cytb_Deer_ked'], errors='coerce').fillna(0) > 0
            both = int((bart & ked).sum())
            say(f"      of which {both} also carry deer ked reads "
                f"({100*both/max(pos,1):.0f}% of Bartonella positives); "
                f"ked is {int(ked.sum())}/{n}")
say("  Reported in Results 3.7 as a co-detection, not in 3.2 or 3.3.")

genera = sorted(set(PATHOGEN_GENUS[c] for c in present))
gmat = pd.DataFrame(index=new.index)
for g in genera:
    cols = [c for c in present if PATHOGEN_GENUS[c] == g]
    gmat[g] = (new[cols].apply(pd.to_numeric, errors='coerce').fillna(0) > 0).any(axis=1)
ngen = gmat.sum(axis=1)

say(f"\n=== co-infection, distinct pathogen GENERA (n = {n}) ===")
say("  Bartonella EXCLUDED as a co-detection (see CO_DETECTIONS above).")
dist = ngen.value_counts().sort_index()
for k, v in dist.items():
    say(f"  {k} genera: {v:>3}  ({100*v/n:4.1f}%)")
say(f"  >=1 genus : {int((ngen>=1).sum()):>3}  ({100*(ngen>=1).mean():4.1f}%)")
say(f"  >=2 genera: {int((ngen>=2).sum()):>3}  ({100*(ngen>=2).mean():4.1f}%)")

# ---- 5. the figure ---------------------------------------------------------
# Top-level import on purpose. A swallowed plotting failure would leave a report
# that reads like a successful run with no figure in it, which is exactly the
# silent-wrong-answer case I want to avoid.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if True:   # keeps the plotting block indented
    BLUE, GREEN, GREY = "#0072B2", "#009E73", "#999999"
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5))

    order = sorted(present, key=lambda c: prev[c][0], reverse=True)
    labels = [c.replace('16S_', '').replace('18S_', '').replace('_', ' ') for c in order]
    vals   = [prev[c][1] for c in order]
    cnts   = [prev[c][0] for c in order]
    ypos   = range(len(order))
    axA.barh(list(ypos), vals, color=BLUE)
    axA.set_yticks(list(ypos)); axA.set_yticklabels(labels, fontstyle='italic')
    axA.invert_yaxis()
    axA.set_xlabel(f"Prevalence (%), n = {n}")
    axA.set_title("A. Pathogen prevalence in single ticks", loc='left', fontweight='bold')
    for i, (v, c) in enumerate(zip(vals, cnts)):
        axA.text(v + 1, i, f"n={c}", va='center', fontsize=9, color='#444444')
    axA.set_xlim(0, max(vals) * 1.25 if vals else 1)
    for sp in ('top', 'right'):
        axA.spines[sp].set_visible(False)

    ks = list(range(0, int(ngen.max()) + 1))
    pcts = [100 * (ngen == k).mean() for k in ks]
    axB.bar([str(k) for k in ks], pcts, color=GREEN)
    for k, p in zip(ks, pcts):
        axB.text(str(k), p + 0.8, f"{p:.1f}%", ha='center', fontsize=9)
    axB.set_xlabel("Co-infecting pathogen genera per tick")
    axB.set_ylabel("Ticks (%)")
    axB.set_title("B. Co-infection (distinct pathogen genera)", loc='left', fontweight='bold')
    axB.set_ylim(0, max(pcts) * 1.2 if pcts else 1)
    axB.text(0.97, 0.95,
             f"{100*(ngen>=1).mean():.1f}% carry \u22651\n{100*(ngen>=2).mean():.1f}% carry \u22652",
             transform=axB.transAxes, ha='right', va='top', fontsize=10,
             bbox=dict(boxstyle='round', fc='#EEEEEE', ec='none'))
    for sp in ('top', 'right'):
        axB.spines[sp].set_visible(False)

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        p = os.path.join(FDIR, f"fig2_single_tick_prevalence_coinfection.{ext}")
        fig.savefig(p, dpi=300, bbox_inches='tight')
        say(f"wrote {p}")
    plt.close(fig)

say("\nNOTE: panel B bars are PERCENTAGES. The old figure labelled the bars with")
say("COUNTS while the axis read '%', which was a real labelling bug.")

with open(os.path.join(ODIR, "single_tick_v3_report.txt"), "w") as fh:
    fh.write("\n".join(rep))
say(f"\nreport -> {ODIR}/single_tick_v3_report.txt")
