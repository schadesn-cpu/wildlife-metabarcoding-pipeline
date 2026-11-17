#!/usr/bin/env python3
"""I collapse asv_table_16S.tsv to one row per taxon for sharing.

This is a summary for people who want the species list without the ASV-level
detail. It is derived, not primary: asv_table_16S.tsv stays the source of truth.
"""
import os
import sys

import pandas as pd

SRC = "results/tables/asv_table_16S.tsv"
OUT = "results/tables/taxa_16S_summary.tsv"

RANKS = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
         "f": "family", "g": "genus", "s": "species"}


def deepest_named(lineage):
    """I return the deepest rank that actually carries a name."""
    if not isinstance(lineage, str):
        return "(unassigned)"
    named = [part for part in lineage.split(";") if part and not part.endswith("__")]
    return named[-1] if named else "(unassigned)"


def main():
    if not os.path.isfile(SRC):
        sys.exit("input not found: " + SRC)
    table = pd.read_csv(SRC, sep="\t", dtype=str)
    for column in ("taxon_lineage", "category", "total_reads", "n_libraries_present"):
        if column not in table.columns:
            sys.exit("column '" + column + "' missing from " + SRC)

    table["total_reads"] = pd.to_numeric(table["total_reads"], errors="coerce").fillna(0)
    table["n_libraries_present"] = pd.to_numeric(
        table["n_libraries_present"], errors="coerce").fillna(0)
    table["taxon"] = table["taxon_lineage"].map(deepest_named)

    grouped = table.groupby("taxon").agg(
        n_asvs=("feature_id", "size"),
        total_reads=("total_reads", "sum"),
        max_libraries_present=("n_libraries_present", "max"),
        category=("category", lambda values: ";".join(
            sorted(set(v for v in values if isinstance(v, str))))),
        lineage=("taxon_lineage", "first"),
    ).reset_index()

    grouped["rank"] = grouped["taxon"].map(
        lambda t: RANKS.get(t.split("__")[0], "unresolved") if "__" in t else "unresolved")
    total = grouped["total_reads"].sum()
    grouped["pct_of_reads"] = (100.0 * grouped["total_reads"] / total).round(4)
    grouped["taxon"] = grouped["taxon"].str.replace(r"^[a-z]__", "", regex=True)

    grouped = grouped[["taxon", "rank", "category", "n_asvs", "total_reads",
                       "pct_of_reads", "max_libraries_present", "lineage"]]
    grouped = grouped.sort_values("total_reads", ascending=False)
    grouped.to_csv(OUT, sep="\t", index=False)

    print("taxa: " + str(len(grouped)))
    print("total reads: " + str(int(total)))
    for name in ("Rickettsia", "Anaplasma", "Borreliella", "Bartonella", "Wolbachia", "Ehrlichia"):
        row = grouped[grouped["taxon"] == name]
        if len(row):
            row = row.iloc[0]
            print("  " + name.ljust(13) + str(int(row["total_reads"])).rjust(10)
                  + " reads  " + str(row["n_asvs"]).rjust(5) + " ASVs  " + str(row["category"]))
    print("wrote " + OUT)


if __name__ == "__main__":
    main()
