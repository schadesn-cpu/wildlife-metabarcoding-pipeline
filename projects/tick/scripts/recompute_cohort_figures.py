#!/usr/bin/env python3
"""I recompute the figures in Section 3.3 that the cohort change moved.

WHY
---
Excluding five libraries that held keds and no ticks took the deer cohort
from 335 animals to 330. Most of Section 3.3 reads its numbers from
core_microbiome_ixodes_deer.tsv, which has been rebuilt. Three do not:

  the count of animals contributing more than one qualifying library,
  Rickettsia as a share of non-host reads rather than of classified
  bacterial reads, and the two-ASV sentence with its read totals and
  per-animal occurrence.

I compute those three on the current cohort so the paragraph can be closed
rather than marked. I reuse build_core_microbiome for the counts, the host
filtering and the cohort, so these are the same 330 animals and the same
reads as the rebuilt table rather than a parallel calculation.

A NUMBER I CANNOT REPRODUCE
---------------------------
Paragraphs 41 and 70 refer to "343 deer whose Ixodes library cleared the
read floor, the set from which these 335 are drawn". The pipeline as it
stands has one floor, applied to non-host reads, which takes the deer
carrying Ixodes straight to the cohort. I report both numbers it does
produce and name the gap rather than inventing a replacement for 343.

I write nothing. Everything goes to stdout.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import build_core_microbiome as core
import check_batch_arms as arms

DROP = ("Mitochondria", "Chloroplast", "Unassigned")


def main():
    if not os.path.isdir("results/analysis"):
        sys.exit("run me from the project root.")
    table, module = core.lookup()
    counts = core.read_counts(table, module)
    animals = [core.animal_of(table, module, s) for s in counts.index]

    # the same block build_core_microbiome.main() runs
    _, _, mapping = arms.read_taxonomy(list(counts.columns))
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in ("mitochondri", "chloroplast", "eukaryota"))]
    counts = counts.drop(columns=host_like)
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})
    deer, moose, species = core.group_members(table, module)
    chosen = core.one_per_animal(counts, animals)
    deer_chosen = chosen[chosen["animal"].isin(deer)]
    print()
    print("cohort: %d deer carrying Ixodes and not Dermacentor, of which %d "
          "have a library above the %d-read floor"
          % (len(deer), len(deer_chosen), core.MIN_READS))

    # --- 1. animals contributing more than one qualifying library ---------
    depth = counts.sum(axis=1)
    per_animal = pd.DataFrame({"animal": animals, "depth": depth.values},
                              index=counts.index)
    qualifying = per_animal[(per_animal["animal"].isin(deer))
                            & (per_animal["depth"] >= core.MIN_READS)]
    many = qualifying.groupby("animal").size()
    print()
    print("1. qualifying libraries per animal, within the cohort")
    print("   animals with a qualifying library     %d" % len(many))
    print("   animals contributing more than one    %d" % int((many > 1).sum()))
    print("   (paragraph 41 said 102 of 335)")

    # --- 2. Rickettsia as a share of non-host reads ------------------------
    part = counts.loc[deer_chosen.index]
    by_genus = part.T.groupby(genus).sum().T
    nonhost = part.sum(axis=1)
    rick = by_genus["Rickettsia"] if "Rickettsia" in by_genus.columns else None
    print()
    print("2. Rickettsia against the two denominators")
    if rick is None:
        print("   no Rickettsia column; nothing to report")
    else:
        share_nonhost = 100 * (rick / nonhost)
        kept = by_genus.drop(columns=[c for c in by_genus.columns
                                      if c in DROP or c.endswith("(domain)")])
        share_class = 100 * (kept["Rickettsia"] / kept.sum(axis=1))
        print("   median share of non-host reads            %.1f%%"
              % share_nonhost.median())
        print("   median share of classified bacterial reads %.1f%%"
              % share_class.median())
        print("   mean share of classified bacterial reads   %.1f%%"
              % share_class.mean())
        print("   (paragraph 70 said 91.6%, 98.2% and 84.0%)")

    # --- 3. the two-ASV sentence -------------------------------------------
    rick_features = [f for f in counts.columns if genus.get(f) == "Rickettsia"]
    print()
    print("3. Rickettsia ASVs across the cohort")
    if not rick_features:
        print("   none found")
        return
    block = part[rick_features]
    totals = block.sum(axis=0).sort_values(ascending=False)
    grand = float(totals.sum())
    top = totals.head(2)
    print("   Rickettsia ASVs                     %d" % len(totals))
    print("   Rickettsia reads in the cohort      %.1f million" % (grand / 1e6))
    print("   the two largest                     %.1f and %.1f million"
          % (top.iloc[0] / 1e6, top.iloc[1] / 1e6))
    print("   their share of Rickettsia reads     %.1f%%"
          % (100 * top.sum() / grand if grand else 0))
    for name in top.index:
        present = int((block[name] > 0).sum())
        print("   %-34s occurs in %d of the %d cohort animals"
              % (str(name)[:34], present, len(deer_chosen)))
    print("   (paragraph 70 said 454 ASVs, 12.8 million reads, 6.3 and 6.1")
    print("    million, 97.2%, occurring in 330 and 332 of 343)")


if __name__ == "__main__":
    main()
