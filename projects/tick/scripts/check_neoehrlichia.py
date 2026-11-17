#!/usr/bin/env python3
"""check_neoehrlichia.py  --  READ ONLY. I write nothing.

WHY
---
droplist_report.txt reassigned one ASV, prefix 865afed6, from g__uncultured to
g__Candidatus_Neoehrlichia at 95.5% identity, "off-panel genus, genus level
only". Neoehrlichia mikurensis is an emerging human pathogen (neoehrlichiosis),
so this is worth resolving rather than leaving in a droplist file.

But the species cannot be read off a genus call, and the species is the whole
question:

  N. mikurensis          Eurasian; Ixodes ricinus / I. persulcatus; human
                         pathogen, especially in immunocompromised patients.
                         Not established in North America.
  Ca. N. lotoris         North American; described from raccoon (Procyon
                         lotor); found in US Ixodes and Dermacentor; no known
                         human disease association.

95.5% clears this paper's own genus floor (94.5%, Yarza et al. 2014) by one
point and misses its species floor (98.7%) by more than three. So the honest
ceiling on a sequence alone is "a Neoehrlichia-genus detection", and the prior
for New England deer ticks favours lotoris.

WHAT I DO ABOUT THAT
--------------------
This dataset can bring independent evidence, because it sequences blood meals.
If the Neoehrlichia-positive animals carry raccoon cytb, that points at
lotoris. If they are clean deer with no raccoon signal, it stays open. Same
move Section 3.4 already makes for the Martha's Vineyard Coxiella.

So I report, for every Neoehrlichia-labelled ASV:
  1. its full feature id, taxon string and confidence
  2. reads per sample from the decontam table, and the total
  3. the host on record for each positive animal, from the crosswalk
  4. the cytb blood-meal call for each positive animal, raccoon flagged
  5. the ASV sequence, printed for you to BLAST against nt yourself

I print the sequence rather than writing a fasta so that this stays read-only.

WHAT WOULD ACTUALLY SETTLE THE SPECIES
--------------------------------------
Not this script. A targeted assay: groEL or gltA on the positive extractions,
which is how the European literature separates Neoehrlichia species, or
long-read 16S. Note that in the manuscript's own terms this would be a
genus-level detection reported as such, exactly like Borreliella -- which is
the paper's whole argument about where each marker's resolution stops.
"""
import os
import re
import sys

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                          else "~/meedlab/tick_project")
TAX = "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv"
DECONTAM = "results/decontam/16S/table_decontam_16S.tsv"
SEQS = "results/tables/asv_sequences_16S.fasta"
CROSSWALK = "results/metadata/crosswalk_VERIFIED.tsv"
BLOODMEAL = ["results/cytb/bloodmeal_per_sample_flagged.tsv",
             "results/cytb/bloodmeal_per_sample.tsv"]
NEEDLE = re.compile(r'neoehrlichia', re.I)
RACCOON = re.compile(r'procyon|raccoon', re.I)


def rows(path, comment_prefixes=("#",)):
    """I yield (header, list-of-rows) from a TSV, skipping biom comment lines
    but keeping the one that carries the column names."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        header, out = None, []
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if header is None:
                if line.startswith("# Constructed from biom"):
                    continue
                header = line.lstrip("#").split("\t")
                header[0] = header[0].strip()
                continue
            out.append(line.split("\t"))
        return header, out


def norm(name):
    """The project's join key: lowercase, hyphens out, marker suffix off."""
    key = name.lower().split(".")[0]
    key = re.sub(r'-?16s$|-?18s$|-?cytb$|-?its$', '', key)
    return key.replace("-", "").replace("_", "")


def main():
    os.chdir(ROOT)
    for need in (TAX, DECONTAM):
        if not os.path.isfile(need):
            sys.exit("missing %s" % need)

    # ---- 1. which ASVs carry the label -------------------------------------
    head, tax_rows = rows(TAX)
    tcol = next(i for i, c in enumerate(head) if c.lower() == "taxon")
    hits = [r for r in tax_rows if len(r) > tcol and NEEDLE.search(r[tcol])]
    print("=" * 78)
    print("NEOEHRLICHIA-LABELLED ASVs IN %s" % TAX)
    print("=" * 78)
    if not hits:
        print("  none. The label is not in this taxonomy table.")
        return
    for r in hits:
        print("  %s" % r[0])
        print("      %s" % r[tcol])
        if len(r) > tcol + 1:
            print("      confidence %s" % r[tcol + 1])
    ids = {r[0] for r in hits}
    print()

    # ---- 2. reads per sample ------------------------------------------------
    dhead, drows = rows(DECONTAM)
    samples = dhead[1:]
    per_sample, totals = {}, {}
    for r in drows:
        if r[0] not in ids:
            continue
        counts = {}
        for name, value in zip(samples, r[1:]):
            try:
                v = float(value)
            except ValueError:
                continue
            if v > 0:
                counts[name] = v
        per_sample[r[0]] = counts
        totals[r[0]] = sum(counts.values())

    print("=" * 78)
    print("READS, FROM THE DECONTAM TABLE")
    print("=" * 78)
    if not per_sample:
        print("  These ASVs are not in the decontam table. They were removed")
        print("  as contaminants or fell out of the feature table, in which")
        print("  case there is nothing to report and that is the answer.")
        return
    for fid, counts in per_sample.items():
        print("  %s" % fid)
        print("      %d reads across %d libraries"
              % (totals[fid], len(counts)))
        for name, v in sorted(counts.items(), key=lambda kv: -kv[1]):
            print("        %-52s %10d" % (name, v))
    print()

    positives = sorted({s for c in per_sample.values() for s in c})
    keys = {norm(s): s for s in positives}

    # ---- 3. host on record --------------------------------------------------
    print("=" * 78)
    print("HOST ON RECORD FOR EACH POSITIVE LIBRARY")
    print("=" * 78)
    if os.path.isfile(CROSSWALK):
        chead, crows = rows(CROSSWALK)
        print("  crosswalk columns: %s" % ", ".join(chead[:14]))
        hostcols = [i for i, c in enumerate(chead)
                    if re.search(r'host', c, re.I)]
        idcols = [i for i, c in enumerate(chead)
                  if re.search(r'sample|barcode|librar|^num$', c, re.I)]
        shown = 0
        for r in crows:
            cell = {norm(r[i]) for i in idcols if i < len(r) and r[i]}
            if not (cell & set(keys)):
                continue
            hosts = ", ".join("%s=%s" % (chead[i], r[i])
                              for i in hostcols if i < len(r) and r[i])
            print("    %-46s %s" % (r[idcols[0]] if idcols else "?", hosts))
            shown += 1
        if not shown:
            print("    no crosswalk row matched. Join by hand from the "
                  "library names above.")
    else:
        print("  %s not found" % CROSSWALK)
    print()

    # ---- 4. cytb blood meal, raccoon flagged -------------------------------
    print("=" * 78)
    print("CYTB BLOOD MEAL FOR EACH POSITIVE ANIMAL")
    print("=" * 78)
    print("  Raccoon here would point at Ca. N. lotoris rather than")
    print("  N. mikurensis. Absence of raccoon does not settle it.")
    print()
    bm = next((p for p in BLOODMEAL if os.path.isfile(p)), None)
    if bm:
        bhead, brows = rows(bm)
        print("  from %s" % bm)
        print("  columns: %s" % ", ".join(bhead[:12]))
        found = False
        for r in brows:
            if norm(r[0]) not in keys:
                continue
            found = True
            line = "\t".join(r[1:8])
            flag = "   <-- RACCOON" if RACCOON.search("\t".join(r)) else ""
            print("    %-40s %s%s" % (r[0], line, flag))
        if not found:
            print("    no blood-meal row for any positive library. Those")
            print("    animals may have no cytb library, or none above the")
            print("    2,500-read floor.")
    else:
        print("  no blood-meal table found at either expected path")
    print()

    # ---- 5. the sequence ----------------------------------------------------
    print("=" * 78)
    print("SEQUENCE, FOR YOU TO BLAST AGAINST nt")
    print("=" * 78)
    if os.path.isfile(SEQS):
        keep, name = False, None
        with open(SEQS, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    keep = name in ids
                    if keep:
                        print("  >%s" % name)
                elif keep:
                    print("  %s" % line)
        print()
        print("  Against nt, megablast. What to look for: whether the top hits")
        print("  are Ca. N. lotoris from North American hosts, N. mikurensis")
        print("  from Eurasian Ixodes, or neither at a convincing identity.")
        print("  At 95.5% to its current best match, do not expect a species.")
    else:
        print("  %s not found" % SEQS)
    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
