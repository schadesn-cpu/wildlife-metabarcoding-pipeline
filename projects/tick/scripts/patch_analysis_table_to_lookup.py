#!/usr/bin/env python3
"""patch_analysis_table_to_lookup.py

I put build_analysis_table.py on the shared library lookup and bring its
pathogen list into line with the rest of the project. It is the last consumer
carrying a private resolver.

THREE CHANGES
-------------
1. THE JOIN. cwmap is replaced by results/metadata/library_lookup.tsv. The
   private version was the most careful of them, and its own comment says why:
   it casefolds and strips the marker tag by regex, after the join "has been
   got wrong three times: 2026-08-05 x2, 2026-08-06". But MARKER_TAG is

       -(16s|18s|cytob|cb)([A-Za-z]*\\d*)(-\\d+)?$

   which strips '-16s-2' and leaves a bare '-2' alone, so a replicate suffix
   only matches when the crosswalk carries that exact string.

   Measured 2026-09-18 by diagnose_analysis_table_join.py, on 672 qualifying
   samples: cwmap resolved 623, the lookup resolves 667. 44 samples across 10
   specimens were being dropped, ZERO of them single ticks, all of them the
   replicate-suffix shape (EQ02441291-2 through -10 and the like).

   analysis_table_16S.tsv therefore goes from 558 samples to 602, and every
   pooled analysis downstream of it was short by the same 44: PERMANOVA,
   PERMDISP, pathogen_cooccurrence.py, fig1a, fig1b.

2. BARTONELLA leaves PATHOGENS for a CO_DETECTIONS list. Section 7 of the
   manuscript classifies these as ruminant-associated Bartonella, not
   tick-borne human pathogens, and build_single_tick_set.py was changed to
   match on 2026-09-18 (ground truth 1.16). This script still called it "real,
   reportable". Reclassified, not deleted: the reads and the call are still
   written, under their own heading.

3. THE CLOSING NOTE says the locked single-tick figure is n=91. It is 93, and
   has been since the 2026-08-07 exclusion rule. I also point it at the
   current canonical file rather than describing it.

WHAT MUST NOT CHANGE
--------------------
No single tick was among the 44, so nothing in the single-tick cohort moves.
The per-pathogen positivity percentages will shift slightly because the
denominator grows from 558 to 602.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_analysis_table_to_lookup.py
    python3 scripts/build_analysis_table.py
"""
import os
import sys

TARGET = "scripts/build_analysis_table.py"
NEEDS = "scripts/tick_lookup.py"

OLD_DOC = """PATHOGENS (real, reportable) — Anaplasma, Borreliella, Borrelia, Bartonella, Ehrlichia (16S)
  NOT counted as pathogens (verified ENDOSYMBIONTS): Rickettsia, Francisella. Reported separately."""

NEW_DOC = """PATHOGENS (real, reportable) — Anaplasma, Borreliella, Borrelia, Ehrlichia (16S)
  NOT counted as pathogens (verified ENDOSYMBIONTS): Rickettsia, Francisella. Reported separately.
  NOT counted as a pathogen (CO-DETECTION, 2026-09-18): Bartonella. Section 7 reports these as
  ruminant-associated and ked-associated organisms, outside the pathogen results; 6 of the 10
  Bartonella-positive single ticks also carry deer ked reads (OR 23.4, p = 0.0001). Still counted
  and still written out, under its own heading. See ground truth 1.16."""

OLD_LISTS = """PATHOGENS   = ['Anaplasma','Borreliella','Borrelia','Bartonella','Ehrlichia']
ENDOSYMBIONT= ['Rickettsia','Francisella']   # verified endosymbionts — NOT pathogens"""

NEW_LISTS = """PATHOGENS   = ['Anaplasma','Borreliella','Borrelia','Ehrlichia']
ENDOSYMBIONT= ['Rickettsia','Francisella']   # verified endosymbionts — NOT pathogens
CO_DETECTION= ['Bartonella']                 # ked/blood-meal associated — NOT a pathogen call"""

OLD_JOIN = """    # crosswalk: sample -> host/pool/barcode
    # crosswalk_4marker is one row per SPECIMEN; name_16S holds every 16S library
    # for it, semicolon-delimited, sometimes with a '-16s' suffix, never with the
    # '.{run}' suffix the feature table uses. Explode to one key per library and
    # register both the bare and '-16s' forms (this join has been got wrong three
    # times: 2026-08-05 x2, 2026-08-06).
    cw = pd.read_csv(CW, sep='\\t')
    cwmap = {}
    for _, r in cw.iterrows():
        rec = {'Barcode': r['sample_id'], 'Host': r['host_spp'],
               'Host_clean': r['host_spp'], 'Host_group': r['host_bin'],
               'NumTicks': r['num_ticks'],
               'is_single_tick': bool(r.get('is_single_tick', False))}
        for nm in str(r['name_16S']).split(';'):
            nm = nm.strip()
            if not nm or nm.lower() == 'nan':
                continue
            # Casefold and strip the marker tag by regex. The old
            # nm.endswith('-16s') test was case-sensitive and dropped 15 of the
            # 16 specimens whose crosswalk names are lowercase, 11 of them Deer.
            stem = MARKER_TAG.sub('', nm)
            for form in (nm, stem, stem + '-16s'):
                cwmap[form.lower()] = rec
    print(f"crosswalk: {len(cw)} specimens -> {len(cwmap)} library keys")"""

NEW_JOIN = """    # Name resolution is NOT done here. build_4marker_crosswalk.py resolves every
    # name in the workbook and writes results/metadata/library_lookup.tsv; this
    # script reads it through scripts/tick_lookup.py.
    #
    # The private cwmap that used to sit here was the most careful of the local
    # resolvers, and it still lost 44 samples across 10 specimens. MARKER_TAG
    # strips '-16s-2' but leaves a bare '-2' alone, so the replicate-suffix shape
    # (EQ02441291-2 through -10) matched only when the crosswalk carried that
    # exact string. Measured 2026-09-18: 623 of 672 qualifying samples against
    # the lookup's 667. None of the 44 were single ticks.
    #
    # Do not reintroduce a local resolver. An unresolved name is a crosswalk
    # problem and is reported below by name.
    sys.path.insert(0, 'scripts')
    import tick_lookup
    lookup = tick_lookup.load()

    cw = pd.read_csv(CW, sep='\\t')
    specimen_meta = {}
    for _, r in cw.iterrows():
        specimen_meta[str(r['sample_id'])] = {
            'Barcode': r['sample_id'], 'Host': r['host_spp'],
            'Host_clean': r['host_spp'], 'Host_group': r['host_bin'],
            'NumTicks': r['num_ticks'],
            'is_single_tick': bool(r.get('is_single_tick', False))}
    print(f"crosswalk: {len(cw)} specimens; resolving through library_lookup.tsv")

    def resolve(sample_id, base_id):
        \"\"\"I turn a feature-table sample id into the crosswalk record, or None.\"\"\"
        for marker, value in (('16S', sample_id), ('16S', base_id),
                              ('specimen', base_id)):
            hit = tick_lookup.resolve(lookup, marker, value)
            if hit is not None:
                return specimen_meta.get(str(hit))
        return None"""

OLD_LOOKUP = """        # 2026-08-07: crosswalk_4marker keys on BASE library names (no '.{run}'
        # suffix), the old crosswalk_VERIFIED keyed on the full suffixed name.
        # Try both so this works with either.
        meta=cwmap.get(str(s).lower()) or cwmap.get(str(base).lower())
        if meta is None or pd.isna(meta['Host']) or pd.isna(meta['NumTicks']):
            excluded['no_host_or_pool']+=1; continue"""

NEW_LOOKUP = """        meta = resolve(s, base)
        # These two exclusions used to share one counter, which is how a join
        # bug hid inside a count that also held legitimate exclusions. A name
        # that resolves to nothing is a BUG; a specimen with no recorded host or
        # tick count is a DATA GAP and a correct exclusion.
        if meta is None:
            excluded['unresolved_name']+=1
            unresolved.append(base)
            continue
        if pd.isna(meta['Host']) or pd.isna(meta['NumTicks']):
            excluded['no_host_or_pool']+=1; continue"""

OLD_COUNTERS = """    rows=[]
    excluded={'control':0,'fungal_trial':0,'below_thresh':0,'no_host_or_pool':0}"""

NEW_COUNTERS = """    rows=[]
    excluded={'control':0,'fungal_trial':0,'below_thresh':0,
              'unresolved_name':0,'no_host_or_pool':0}
    unresolved=[]"""

OLD_ENDO_LOOP = """        for e,ids in endo_ids.items():
            reads=sum(col[idx[f]] for f in ids if f in idx)
            rec[e+'_reads']=int(reads)   # endosymbiont reads (context, not a pathogen call)
        rows.append(rec)"""

NEW_ENDO_LOOP = """        for e,ids in endo_ids.items():
            reads=sum(col[idx[f]] for f in ids if f in idx)
            rec[e+'_reads']=int(reads)   # endosymbiont reads (context, not a pathogen call)
        for c,ids in codet_ids.items():
            reads=sum(col[idx[f]] for f in ids if f in idx)
            rec[c]=int(reads); rec[c+'_pos']=int(reads>=MIN_READS_CALL)
        rows.append(rec)"""

OLD_IDS = """    path_ids={p: ids_for(p) for p in PATHOGENS}
    endo_ids={e: ids_for(e) for e in ENDOSYMBIONT}"""

NEW_IDS = """    path_ids={p: ids_for(p) for p in PATHOGENS}
    endo_ids={e: ids_for(e) for e in ENDOSYMBIONT}
    codet_ids={c: ids_for(c) for c in CO_DETECTION}"""

OLD_PREFILTER = """    path_ids={p:(v & features) for p,v in path_ids.items()}
    endo_ids={e:(v & features) for e,v in endo_ids.items()}"""

NEW_PREFILTER = """    path_ids={p:(v & features) for p,v in path_ids.items()}
    endo_ids={e:(v & features) for e,v in endo_ids.items()}
    codet_ids={c:(v & features) for c,v in codet_ids.items()}"""

OLD_REPORT = """    print("\\n(Rickettsia/Francisella are endosymbionts — reads recorded as context, NOT pathogen calls)")"""

NEW_REPORT = """    print("\\n(Rickettsia/Francisella are endosymbionts — reads recorded as context, NOT pathogen calls)")
    if unresolved:
        print(f"\\n!! {len(unresolved)} sample(s) did not resolve to any crosswalk specimen.")
        print("   These are NOT controls and NOT data gaps; they are unresolved names.")
        for u in sorted(set(unresolved))[:20]:
            print("      " + str(u))
        if len(set(unresolved)) > 20:
            print(f"      ... and {len(set(unresolved))-20} more")
        print("   Fix the crosswalk; do not add a fallback here.")
    print("\\n=== co-detections (NOT pathogens; see ground truth 1.16) ===")
    for c in CO_DETECTION:
        pos=df[c+'_pos'].sum()
        print(f"  {c:14s}: {pos}/{len(df)} = {100*pos/len(df):.1f}%")
    print("  Bartonella is ked-associated and blood-meal associated. Reported in")
    print("  Results 3.7 as a co-detection, not in 3.2 and not in any co-infection")
    print("  count.")"""

OLD_NOTE = '''    print("\\nNOTE: single-tick prevalences here use the >=1000 filter; the LOCKED single-tick figures")
    print("use n=91 (single_tick_v3, low floor). Expect small differences — single_tick_v3_full.tsv is canonical for")
    print("single-tick descriptives. This table is for the POOLED analyses.")'''

NEW_NOTE = '''    print("\\nNOTE: single-tick prevalences here use the >=1000 filter; the LOCKED single-tick")
    print("figures use n=93 (results/single_tick_v3/single_tick_v3_full.tsv, low floor), not 91.")
    print("The cohort has been 93 since the 2026-08-07 exclusion rule; see ground truth 1.6.")
    print("Expect small differences. single_tick_v3_full.tsv is canonical for single-tick")
    print("descriptives. This table is for the POOLED analyses.")'''

REPLACEMENTS = [
    (OLD_DOC, NEW_DOC, "the docstring PATHOGENS block"),
    (OLD_LISTS, NEW_LISTS, "the PATHOGENS / ENDOSYMBIONT lists"),
    (OLD_IDS, NEW_IDS, "the id-set construction"),
    (OLD_PREFILTER, NEW_PREFILTER, "the id pre-filter"),
    (OLD_JOIN, NEW_JOIN, "the cwmap construction"),
    (OLD_COUNTERS, NEW_COUNTERS, "the exclusion counters"),
    (OLD_LOOKUP, NEW_LOOKUP, "the per-sample lookup"),
    (OLD_ENDO_LOOP, NEW_ENDO_LOOP, "the endosymbiont read loop"),
    (OLD_REPORT, NEW_REPORT, "the endosymbiont report line"),
    (OLD_NOTE, NEW_NOTE, "the closing n=91 note"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    if not os.path.isfile(NEEDS):
        sys.exit("not found: " + NEEDS + "\nWrite tick_lookup.py first.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "import tick_lookup" in source:
        sys.exit("already patched: " + TARGET + " imports tick_lookup. "
                 "Nothing changed.")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("expected exactly one occurrence of " + label
                     + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    if "import sys" not in source.split("def main")[0]:
        source = source.replace("import pandas as pd, numpy as np",
                                "import sys\nimport pandas as pd, numpy as np", 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print()
    print("Before this patch: 558 samples in analysis_table_16S.tsv.")
    print("Expect 602. None of the 44 recovered samples are single ticks, so")
    print("nothing in the single-tick cohort moves; only the pooled denominator.")
    print()
    print("run: python3 " + TARGET)


if __name__ == "__main__":
    main()
