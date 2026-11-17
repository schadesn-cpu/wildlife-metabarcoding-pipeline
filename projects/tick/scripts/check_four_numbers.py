#!/usr/bin/env python3
"""check_four_numbers.py  --  READ ONLY. I write nothing, I change nothing.

WHAT THIS IS FOR
----------------
Four numbers in the manuscript were computed before the two duplicate bear
accessions were merged, and they cannot be fixed by subtracting two, because
each has its own denominator and its own numerator. This recomputes all four
from the package you already have on your laptop, so no ron and no password.

    Section 2.2   "tick size or engorgement for 310 of 561 specimens and are
                   silent for the remaining 253"          310 + 253 = 563
    Section 3.4   "214 of 559 animals (38%)"
    Section 3.4   "207 of 552 other animals (38%)"
    Section 3.4   "the detection appeared in 13 of their 18 tested
                   subsamples against 263 of 807 elsewhere (33%)"
    Section 3.4   "All 7 black bear with a qPCR result"   against Table 1's 10

HOW IT BEHAVES
--------------
It prints the counts and names which published number each one reproduces.
Where a published number matches nothing it says so plainly instead of
picking the nearest. Six times this month a document's account of the work
has turned out not to match the work, so this reports and does not assert.

THE MERGE RULE, STATED SO IT CAN BE ARGUED WITH
-----------------------------------------------
USGSTICKS-VTBearA and USGSTICKS-VT-12-BearA are one animal, and so are the
two BearC accessions. Any animal id containing BearA..BearE collapses to that
letter; everything else is left exactly as it is. The groups are printed
before anything is counted, so if that rule catches something it should not,
you will see it in the first section rather than in a wrong total.

    cd ~/Downloads/usgs_G23AC00567 && python3 check_four_numbers.py
"""
import csv
import os
import re
import sys
from collections import defaultdict

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
RECOUNT = "analysis/qpcr_recount_by_rule.tsv"
SCREEN = "analysis/qpcr_screening_record_corrected.tsv"
COVAR = "metadata/accessioning_covariates.tsv"

PUBLISHED = [
    ("Section 2.2  size or engorgement recorded", 310),
    ("Section 2.2  silent", 253),
    ("Section 3.4  B. microti positive animals", 214),
    ("Section 3.4  animals with a qPCR result", 559),
    ("Section 3.4  B. microti positive, non-bear", 207),
    ("Section 3.4  non-bear animals", 552),
    ("Section 3.4  bear subsamples with the detection", 13),
    ("Section 3.4  bear subsamples tested", 18),
    ("Section 3.4  other subsamples with the detection", 263),
    ("Section 3.4  other subsamples tested", 807),
    ("Section 3.4  black bear with a qPCR result", 7),
    ("Table 1      black bear accessions", 10),
]

bar = "=" * 78


def rows(path):
    full = os.path.join(ROOT, path)
    if not os.path.isfile(full):
        sys.exit("missing %s\nRun me from inside the extracted package." % full)
    with open(full, encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def key(animal):
    """Collapse the duplicate bear accessions and leave everything else alone."""
    a = (animal or "").strip()
    m = re.search(r'bear\s*([a-e])(?![a-z])', a, re.I)
    return "BEAR_" + m.group(1).upper() if m else a


def truthy(v):
    v = (v or "").strip()
    if v == "" or v.lower() in ("na", "nan", "none", "null"):
        return False
    return True


def positive(v):
    v = (v or "").strip()
    try:
        return float(v) > 0
    except ValueError:
        return v.lower() in ("true", "yes", "y", "pos", "positive")


found = {}
# animal key -> is this animal a black bear, decided from the host columns of
# the recount table rather than from its accession number. USGSTICKS-LHT113 is
# a black bear whose id does not say so, and counting bears by name misses it.
BEAR_BY_HOST = {}


def note(label, value):
    found.setdefault(value, label)


# ---------------------------------------------------------------- 1. bears
print(bar)
print("1. WHICH ANIMALS THE MERGE RULE TOUCHES")
print(bar)
rec = rows(RECOUNT)
animals = sorted({r["animal"].strip() for r in rec if r.get("animal")})
groups = defaultdict(list)
for a in animals:
    groups[key(a)].append(a)
merged = {k: v for k, v in groups.items() if len(v) > 1}
print("  %d accessions in %s" % (len(animals), RECOUNT))
print("  %d of them group into %d animals" %
      (sum(len(v) for v in merged.values()), len(merged)))
for k, v in sorted(merged.items()):
    print("     %-10s  <-  %s" % (k, "  +  ".join(v)))
bearish = [a for a in animals if re.search(r'bear', a, re.I)]
print("\n  every accession whose id mentions a bear:")
for a in bearish:
    print("     %-34s -> %s" % (a, key(a)))
if not merged:
    print("\n  *** the merge rule matched nothing. The animal column may not")
    print("      hold barcodes. Check the ids above before trusting part 3. ***")
print()

# --------------------------------------------------- 2. the 559 denominator
print(bar)
print("2. THE qPCR DENOMINATOR, BEFORE AND AFTER THE MERGE")
print(bar)
targets = sorted({r["target"].strip() for r in rec if r.get("target")})
before = len(animals)
after = len(groups)
print("  targets in the table: %d   %s" % (len(targets), ", ".join(targets)))
print("  rows: %d   =  %d accessions x %d targets" % (len(rec), before, len(targets)))
print("  animals with a qPCR result, as counted now : %d" % before)
print("  animals with a qPCR result, after the merge: %d" % after)
note("animals with a qPCR result, pre-merge", before)
note("animals with a qPCR result, post-merge", after)
print()

# ----------------------------------------------- 3. B. microti, three rules
print(bar)
print("3. BABESIA MICROTI POSITIVE ANIMALS UNDER EACH COUNTING RULE")
print(bar)
RULES = [c for c in ("detected_recorded", "detected_stated", "detected_narrow")
         if c in (rec[0].keys() if rec else ())]
if not RULES:
    print("  none of the detected_* columns are present; columns are:")
    print("    %s" % ", ".join(rec[0].keys()))
else:
    micro = [r for r in rec if re.search(r'microti', r.get("target", ""), re.I)]
    print("  %d rows carry a B. microti target" % len(micro))
    host_of = {}
    for r in micro:
        host_of[key(r["animal"])] = "%s %s" % (r.get("host_bin", ""),
                                               r.get("host_spp", ""))
    is_bear = {k: bool(re.search(r'bear', v, re.I)) or k.startswith("BEAR_")
               for k, v in host_of.items()}
    BEAR_BY_HOST.update(is_bear)
    nbear = sum(1 for k in host_of if is_bear[k])
    print("  %d of the %d animals are black bear after the merge" % (nbear, len(host_of)))
    note("black bear with a qPCR result, post-merge", nbear)
    print()
    print("  %-20s %12s %12s %12s %12s" %
          ("rule", "pos pre", "pos post", "non-bear", "non-bear n"))
    for rule in RULES:
        pre = sum(1 for r in micro if positive(r.get(rule)))
        post_set = {key(r["animal"]) for r in micro if positive(r.get(rule))}
        nonbear_pos = sum(1 for k in post_set if not is_bear.get(k))
        nonbear_n = sum(1 for k in host_of if not is_bear[k])
        print("  %-20s %12d %12d %12d %12d"
              % (rule.replace("detected_", ""), pre, len(post_set),
                 nonbear_pos, nonbear_n))
        note("B. microti positive (%s), pre-merge" % rule, pre)
        note("B. microti positive (%s), post-merge" % rule, len(post_set))
        note("B. microti positive non-bear (%s)" % rule, nonbear_pos)
        note("non-bear animals with a qPCR result", nonbear_n)
print()

# ------------------------------------------------------- 4. the subsamples
print(bar)
print("4. TESTED SUBSAMPLES, BEAR AGAINST EVERYTHING ELSE")
print(bar)
scr = rows(SCREEN)
print("  %d rows in %s   (Section 2.5 says the screening record holds 568)"
      % (len(scr), SCREEN))
note("rows in the corrected screening record", len(scr))
sub_by = defaultdict(int)
mic_by = defaultdict(int)
for r in scr:
    k = key(r.get("Barcode", ""))
    try:
        n = int(float((r.get("subsample_count") or "0").strip() or 0))
    except ValueError:
        n = 0
    sub_by[k] += n
    if positive(r.get("b-microti")):
        mic_by[k] += 1
if BEAR_BY_HOST:
    bear_keys = [k for k in sub_by if BEAR_BY_HOST.get(k)]
    how = "host columns of the recount table"
else:
    bear_keys = [k for k in sub_by if k.startswith("BEAR_")]
    how = "accession name only, because the host lookup came back empty"
bsub = sum(sub_by[k] for k in bear_keys)
osub = sum(v for k, v in sub_by.items() if k not in bear_keys)
bmic = sum(mic_by[k] for k in bear_keys)
omic = sum(v for k, v in mic_by.items() if k not in bear_keys)
print("  bears identified by: %s" % how)
print("  bear accessions here: %s" % (", ".join(sorted(bear_keys)) or "(none)"))
print("  %-44s %6d" % ("subsamples tested, bear", bsub))
print("  %-44s %6d" % ("subsamples tested, other", osub))
print("  %-44s %6d" % ("subsamples tested, total", bsub + osub))
print("  %-44s %6d" % ("B. microti positive subsamples, bear", bmic))
print("  %-44s %6d" % ("B. microti positive subsamples, other", omic))
note("bear subsamples tested", bsub)
note("other subsamples tested", osub)
note("bear subsamples with the B. microti detection", bmic)
note("other subsamples with the B. microti detection", omic)
print()

# ----------------------------------------------- 5. size and engorgement
print(bar)
print("5. SIZE AND ENGORGEMENT, WHICH DEFINITION GIVES 310")
print(bar)
cov = rows(COVAR)
print("  %d rows in %s" % (len(cov), COVAR))
note("rows in accessioning_covariates", len(cov))
ckeys = {key(r.get("animal", "")) for r in cov if r.get("animal")}
print("  %d distinct animals after the merge" % len(ckeys))
note("animals in accessioning_covariates, post-merge", len(ckeys))
print()
CANDIDATES = [
    ("engorged_count alone", lambda r: truthy(r.get("engorged_count"))),
    ("male_count alone", lambda r: truthy(r.get("male_count"))),
    ("engorged_count or male_count", lambda r: truthy(r.get("engorged_count"))
     or truthy(r.get("male_count"))),
    ("sex_census_complete set", lambda r: truthy(r.get("sex_census_complete"))),
    ("any of the three", lambda r: truthy(r.get("engorged_count"))
     or truthy(r.get("male_count")) or truthy(r.get("sex_census_complete"))),
]
print("  %-32s %10s %10s %10s" % ("definition", "rows", "animals", "silent"))
for label, test in CANDIDATES:
    hit = {key(r.get("animal", "")) for r in cov if test(r)}
    nrows = sum(1 for r in cov if test(r))
    print("  %-32s %10d %10d %10d"
          % (label, nrows, len(hit), len(ckeys) - len(hit)))
    note("size/engorgement, %s, rows" % label, nrows)
    note("size/engorgement, %s, animals" % label, len(hit))
    note("silent, %s" % label, len(ckeys) - len(hit))
print()

# ------------------------------------------------------------ 6. the verdict
print(bar)
print("6. WHICH PUBLISHED NUMBER IS WHICH")
print(bar)
for label, value in PUBLISHED:
    match = found.get(value)
    print("  %-46s %5d   %s"
          % (label, value, match or "*** reproduces nothing above ***"))
print()
print("  A number that reproduces nothing is the one to chase. A number that")
print("  matches a pre-merge cell needs its post-merge twin from the same row.")
print()
print("Send this whole output back.")
