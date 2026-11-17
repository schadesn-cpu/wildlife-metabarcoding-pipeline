#!/usr/bin/env python3
"""build_s2_table.py -- S2 Table, the two SILVA classifiers compared.

    cd ~/meedlab/tick_project
    python3 -I scripts/build_s2_table.py

THIS ONE WRITES. Everything else I have sent you today was read-only. It
creates results/tables/S2_classifier_comparison.tsv and refuses to run if that
file already exists, so it cannot overwrite anything.

It reads the two exported relative-abundance tables, pulls the genus out of
each lineage string, and reports per genus under each classifier. These are
relative abundances, so there is no read count here and no "reads in dispute"
figure. Section 3.3 currently claims 0.001% of reads in dispute, which these
files cannot support either way.
"""
import os
import sys

FULL = "results/16S/all/taxonomy/relabund_taxonomy.tsv"
REGION = "results/16S/all/taxonomy/relabund_taxonomy_v4v5_regionmatched.tsv"
OUT = "results/tables/S2_classifier_comparison.tsv"

# every genus named in Sections 3.3 and 3.4
NAMED = ["Rickettsia", "Francisella", "Anaplasma", "Borreliella", "Borrelia",
         "Bartonella", "Arsenophonus", "Escherichia-Shigella", "Ehrlichia",
         "Coxiella", "Nocardia", "Diplorickettsia", "Rickettsiella",
         "Neoehrlichia", "Enhydrobacter", "Alkanindiges", "Chryseobacterium",
         "Deinococcus", "Sphingobacterium"]

if os.path.exists(OUT):
    sys.exit("%s already exists. Move it aside first; I will not overwrite." % OUT)


def load(path):
    header, rows = None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("# ") and "OTU" not in ln:
                continue
            parts = ln.rstrip("\n").split("\t")
            if header is None:
                header = parts[1:]
                continue
            vals = []
            for v in parts[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(0.0)
            rows.append((parts[0], vals))
    return header, rows


def genus(lineage):
    for field in lineage.split(";"):
        f = field.strip()
        if f.startswith("g__"):
            name = f[3:].strip()
            return name or None
    return None


def summarize(path):
    header, rows = load(path)
    print("   %s" % path)
    print("      %d taxa rows, %d sample columns" % (len(rows), len(header)))
    agg = {}
    for lineage, vals in rows:
        g = genus(lineage)
        if not g:
            continue
        n_present = sum(1 for v in vals if v > 0)
        total = sum(vals)
        cur = agg.setdefault(g, [0, 0, 0.0])
        cur[0] += 1                      # taxa rows folded into this genus
        cur[1] = max(cur[1], n_present)  # samples carrying the genus
        cur[2] += total                  # summed relative abundance
    print("      %d named genera" % len(agg))
    return agg, len(header)


print("S2 Table, the two SILVA classifiers compared")
print()
full, n_samples = summarize(FULL)
region, _ = summarize(REGION)
print()

only_full = sorted(set(full) - set(region))
only_region = sorted(set(region) - set(full))
both = sorted(set(full) & set(region))
print("   genera in both classifiers        %d" % len(both))
print("   genera only under full-length     %d" % len(only_full))
print("   genera only under region-matched  %d" % len(only_region))
print()

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as out:
    out.write("genus\tnamed_in_manuscript\tfull_taxa\tfull_samples\t"
              "full_mean_relabund\tregion_taxa\tregion_samples\t"
              "region_mean_relabund\tstatus\n")
    for g in sorted(set(full) | set(region)):
        f = full.get(g, [0, 0, 0.0])
        r = region.get(g, [0, 0, 0.0])
        status = ("both" if g in full and g in region
                  else "full-length only" if g in full else "region-matched only")
        named = "yes" if any(n.lower() in g.lower() for n in NAMED) else ""
        out.write("%s\t%s\t%d\t%d\t%.6f\t%d\t%d\t%.6f\t%s\n" % (
            g, named, f[0], f[1], f[2] / n_samples,
            r[0], r[1], r[2] / n_samples, status))

print("wrote %s" % OUT)
print()
print("The genera named in Sections 3.3 and 3.4")
print("%-24s %22s %22s" % ("", "full-length", "region-matched"))
print("%-24s %10s %11s %10s %11s" % ("genus", "samples", "mean%", "samples", "mean%"))
for g in sorted(set(full) | set(region)):
    if not any(n.lower() in g.lower() for n in NAMED):
        continue
    f = full.get(g, [0, 0, 0.0])
    r = region.get(g, [0, 0, 0.0])
    print("%-24s %10d %10.4f %10d %10.4f" % (
        g[:24], f[1], 100 * f[2] / n_samples, r[1], 100 * r[2] / n_samples))
print()
print("A genus showing samples under one classifier and 0 under the other is")
print("the case Section 3.4's clause is about. If every named genus shows the")
print("same or more under region-matched, the superset claim in 3.3 holds.")
