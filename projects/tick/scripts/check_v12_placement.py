#!/usr/bin/env python3
"""check_vt12_placement.py  --  READ ONLY. I write nothing.

Replaces check_vt12_vs_questing.py, which produced Rickettsia 0.0% for every
VT-12 tick while vt12_report.txt gives 51LC 76.1% and 54WC 94.0%.

WHAT WENT WRONG, SO IT CANNOT RECUR SILENTLY
--------------------------------------------
515F/926R is a universal small-subunit pair. It amplifies the tick's own 18S
along with bacterial 16S, and SILVA classifies that host rDNA as
d__Eukaryota ... Acari. The old script built its denominator by subtracting
only Mitochondria and Chloroplast, so bacterial reads were being divided by a
mostly-tick total and every percentage collapsed. 53TL came out as 99.1%
"Acari" -- which is the tick, not its microbiome.

So this version:
  * classifies on DOMAIN parsed from the taxon string, not on genus name
  * reports the domain split per library, so host co-amplification is visible
    rather than hidden in a denominator
  * refuses to report a composition for any library whose bacterial fraction
    is below the floor
  * CHECKS ITSELF against vt12_report.txt's published figures before printing
    anything else. If 51LC does not come back near 76.1% and 54WC near 94.0%,
    it says so at the top and you should not trust the rest.

That last point is the lesson from the first attempt: a number that is wrong
for a structural reason looks exactly like a number that is right, unless
something independent is checking it.

THE QUESTION
------------
Where do Cheryl's nine coded ticks sit in the collection as a whole?

Not "did the treatment do anything" -- unanswerable, because all nine sit in
one sequencing batch, so batch and group are the same variable and nothing
separates a variable from itself.

The first version compared her ticks only against questing ticks, reasoning
that hers are unfed and the deer cohort is engorged. That picked a comparator
from an assumption and discarded most of the collection to do it. This version
commits to no comparator: it summarises every stratum the metadata offers,
places each of her ticks as a percentile, and MEASURES the fed/unfed axis
through mitochondrial share instead of arguing about it.
"""
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                          else "~/meedlab/tick_project")
TAX = "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv"
TABLE = "results/decontam/16S/table_decontam_16S.tsv"
CROSSWALK = "results/metadata/crosswalk_VERIFIED.tsv"
SINGLE = "results/single_tick_v3/single_tick_host_pathogen.tsv"

VT12 = ["51LC", "51TL", "53TL", "54WC", "57TW", "59TL1", "59TL2", "59TW2", "VTW"]
PANEL = ["Anaplasma", "Borreliella", "Borrelia", "Bartonella", "Ehrlichia",
         "Francisella", "Rickettsia", "Wolbachia"]
# vt12_report.txt, regenerated 2026-10-20. The script must reproduce these.
EXPECT = {"51LC": 76.14, "53TL": 1.61, "54WC": 94.02, "57TW": 78.99,
          "59TL2": 85.70, "VTW": 83.34}
MIN_BACT = 1000
NONBACT = ("Mitochondria", "Chloroplast")


def parse(taxon):
    """I return (domain, lowest named rank). Domain decides the denominator."""
    d = re.search(r'd__([^;]+)', taxon)
    domain = d.group(1).strip() if d and d.group(1).strip() else "Unassigned"
    for rank in ("g__", "f__", "o__", "c__", "p__", "d__"):
        m = re.search(re.escape(rank) + r'([^;]+)', taxon)
        if m and m.group(1).strip():
            return domain, m.group(1).strip()
    return domain, "Unassigned"


def header_of(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("# Constructed from biom"):
                continue
            return line.rstrip("\n").lstrip("#").split("\t")
    return []


def main():
    os.chdir(ROOT)
    for need in (TAX, TABLE):
        if not os.path.isfile(need):
            sys.exit("missing %s" % need)
    bar = "=" * 78

    tax = {}
    with open(TAX, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        tcol = next(i for i, c in enumerate(head) if c.lower() == "taxon")
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) > tcol:
                tax[r[0]] = parse(r[tcol])

    with open(TABLE, encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
        if first.startswith("# Constructed from biom"):
            first = fh.readline()
        samples = first.rstrip("\n").lstrip("#").split("\t")[1:]
        bactg = [defaultdict(int) for _ in samples]
        dom = [defaultdict(int) for _ in samples]
        seen = unmapped = 0
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) < 2:
                continue
            seen += 1
            hit = tax.get(r[0])
            if hit is None:
                unmapped += 1
                hit = ("Unassigned", "Unassigned")
            domain, name = hit
            for k, v in enumerate(r[1:]):
                if v and v not in ("0", "0.0"):
                    try:
                        n = int(float(v))
                    except ValueError:
                        continue
                    if not n:
                        continue
                    dom[k][domain] += n
                    if domain == "Bacteria" and name not in NONBACT:
                        bactg[k][name] += n

    print(bar); print("DOES THE TABLE MATCH THE TAXONOMY?"); print(bar)
    print("  %d features in the table, %d not found in the taxonomy (%.1f%%)"
          % (seen, unmapped, 100.0 * unmapped / max(seen, 1)))
    if unmapped > seen * 0.05:
        print("  *** More than 5% unmapped. The table and the taxonomy are")
        print("      probably not the same feature space. Stop here. ***")
    print()

    def measures(i):
        b = bactg[i]
        bact = sum(b.values())
        tot = sum(dom[i].values())
        euk = dom[i].get("Eukaryota", 0)
        return dict(bact=bact, total=tot,
                    euk_pct=100.0 * euk / tot if tot else 0.0,
                    rick=100.0 * b.get("Rickettsia", 0) / bact if bact else 0.0,
                    rich=sum(1 for g, n in b.items() if n > 0),
                    panel=[g for g in PANEL if b.get(g, 0) > 0])

    idx = {c: next((i for i, s in enumerate(samples)
                    if "VT-12" in s and re.search(r'[-.]%s[-.]' % re.escape(c), s)),
                   None) for c in VT12}

    print(bar); print("SELF-CHECK AGAINST vt12_report.txt"); print(bar)
    bad = 0
    for c, want in sorted(EXPECT.items()):
        i = idx.get(c)
        if i is None:
            print("  %-6s no library found" % c); bad += 1; continue
        got = measures(i)["rick"]
        ok = abs(got - want) < 1.5
        bad += 0 if ok else 1
        print("  %-6s expect %6.2f%%   got %6.2f%%   %s"
              % (c, want, got, "ok" if ok else "MISMATCH"))
    print()
    if bad:
        print("  *** %d mismatches. Something upstream differs from the report;" % bad)
        print("      do not use the numbers below. ***")
        print()
    else:
        print("  All match. The denominator is right and the rest can be trusted.")
        print()

    print(bar); print("CHERYL'S NINE"); print(bar)
    print("  %-7s %10s %9s %9s %7s  %s"
          % ("code", "bact reads", "euk share", "Rickettsia", "genera", "panel"))
    for c in VT12:
        i = idx.get(c)
        if i is None:
            print("  %-7s %10s" % (c, "no library")); continue
        m = measures(i)
        thin = "  THIN" if m["bact"] < MIN_BACT else ""
        print("  %-7s %10d %8.1f%% %8.1f%% %7d  %s%s"
              % (c, m["bact"], m["euk_pct"], m["rick"], m["rich"],
                 ", ".join(m["panel"]) or "-", thin))
    print()
    print("  euk share is the tick's own 18S, co-amplified by 515F/926R. It is")
    print("  excluded from every percentage to its right.")
    print()

    # ---- strata --------------------------------------------------------
    strata, ticks_of = defaultdict(list), {}
    if os.path.isfile(CROSSWALK):
        head = header_of(CROSSWALK)
        ic = next((i for i, c in enumerate(head)
                   if re.search(r'table_sample_16S', c, re.I)), 0)
        gc = next((i for i, c in enumerate(head)
                   if re.search(r'host_group', c, re.I)), None)
        nc = next((i for i, c in enumerate(head)
                   if re.search(r'numticks', c, re.I)), None)
        by = {}
        with open(CROSSWALK, encoding="utf-8", errors="replace") as fh:
            next(fh)
            for line in fh:
                r = line.rstrip("\n").split("\t")
                if ic < len(r) and r[ic].strip():
                    by[r[ic].strip()] = r
        mine = {i for i in idx.values() if i is not None}
        for i, name in enumerate(samples):
            r = by.get(name)
            if not r or i in mine:
                continue
            g = (r[gc].strip() if gc is not None and gc < len(r) else "") or "(blank)"
            strata[g].append(i)
            if nc is not None and nc < len(r):
                try:
                    ticks_of[i] = float(r[nc])
                except ValueError:
                    pass

    def ok(ii):
        return [i for i in ii if measures(i)["bact"] >= MIN_BACT]

    def mid(v):
        return sorted(v)[len(v) // 2]

    def row(label, ii):
        ms = [measures(i) for i in ok(ii)]
        if not ms:
            print("  %-26s  nothing above the floor" % label[:26]); return
        print("  %-26s  n=%-4d Rick %5.1f%% (%.1f-%.1f)  genera %3d  euk %5.1f%%"
              % (label[:26], len(ms), mid([m["rick"] for m in ms]),
                 min(m["rick"] for m in ms), max(m["rick"] for m in ms),
                 mid([m["rich"] for m in ms]), mid([m["euk_pct"] for m in ms])))

    print(bar); print("THE COLLECTION, BY STRATUM"); print(bar)
    for g in sorted(strata, key=lambda k: -len(strata[k])):
        row(g, strata[g])
    singles = [i for i, n in ticks_of.items() if n == 1]
    pooled = [i for i, n in ticks_of.items() if n > 1]
    print()
    row("ALL single ticks", singles)
    row("ALL pooled", pooled)
    print()
    print("  Read the euk column first. If single ticks and pools differ there,")
    print("  collection mode is doing real work and her ticks belong against the")
    print("  comparable stratum. If they do not, the whole collection is a fair")
    print("  reference and my earlier caution was unnecessary.")
    print()

    ref = ok(singles) or ok([i for i in range(len(samples))
                             if i not in set(idx.values())])
    rr = sorted(measures(i)["rick"] for i in ref)
    rh = sorted(measures(i)["rich"] for i in ref)

    def pct(v, x):
        return 100.0 * sum(1 for y in v if y < x) / len(v) if v else float("nan")

    print(bar)
    print("WHERE EACH OF HERS FALLS, against %d reference libraries" % len(ref))
    print(bar)
    for c in VT12:
        i = idx.get(c)
        if i is None:
            continue
        m = measures(i)
        if m["bact"] < MIN_BACT:
            print("  %-7s too thin to place (%d bacterial reads)" % (c, m["bact"]))
            continue
        print("  %-7s Rickettsia %5.1f%% -> %5.1f pct   richness %3d -> %5.1f pct"
              % (c, m["rick"], pct(rr, m["rick"]), m["rich"], pct(rh, m["rich"])))
    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
