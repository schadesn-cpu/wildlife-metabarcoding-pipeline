#!/usr/bin/env python3
"""Does the tick microbiome differ by pathogen-carriage confidence?

Runs the comparison specified in ANALYSIS_PLAN_pathogen_carriage_vs_microbiome.md, on the
69 deer-collected ticks screened by both qPCR and sequencing:

    confirmed_carrier   34   positive by BOTH methods
    single_method       21   positive by one only
    confirmed_negative  14   negative by both

Three decisions are baked in here and each of them was made BEFORE seeing a result. They are
in the code so nobody has to take my word for it later.

1. DEER ONLY. All 34 confirmed carriers are deer-collected and all 8 other-wildlife ticks are
   confirmed negative, so across hosts "carriage" and "deer" are nearly the same variable.
   Deer-collected ticks are treated as Ixodes: 0 of 273 deer specimens carried Dermacentor.

2. PATHOGEN GENERA ARE STRIPPED from the community before testing. Anaplasma, Borreliella,
   Borrelia and Bartonella are members of the community being compared, so leaving them in
   would partly reduce the question to "does the community contain the thing we grouped on".

3. NO BATCH TERM. Only 27 of the 69 carry a batch label, so adjusting drops the analysis to
   11/7/9. Batch is reported as a limitation with its effect size instead. See the plan.

I report PERMDISP whatever PERMANOVA says. With groups of 34/21/14, a significant PERMANOVA
can reflect unequal spread rather than a difference in community composition, and only
PERMDISP distinguishes those.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/carriage_vs_microbiome.py \
        --design results/analysis/carriage_design_check.tsv \
        --outdir results/analysis/carriage_microbiome

Writes:
    <outdir>/carriage_microbiome_report.txt
    <outdir>/alpha_by_group.tsv
    <outdir>/permanova_permdisp.tsv
    <outdir>/genus_by_group.tsv
"""

import argparse
import csv
import glob
import os
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd

TABLE_QZA = "results/decontam/16S/table_decontam_16S.qza"
TAXONOMY = "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv"
DECISIONS = "results/tables/taxon_decision_table.tsv"

# Batch effect from the locked diversity rebuild, reported alongside the carriage result so
# the reader can see the relative size rather than being told it is small.
BATCH_R2_PCT = 10.7
BATCH_P = 0.001

GROUP_ORDER = ["confirmed_negative", "single_method", "confirmed_carrier"]

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def back_up_existing(out_path):
    """Moves an existing output aside to .bak. Refuses to clobber a backup."""
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def load_qza_table(path):
    """Feature table from a .qza, as a DataFrame with ASVs on rows."""
    if not os.path.exists(path):
        sys.exit(f"no such file: {path}")
    import biom
    with tempfile.TemporaryDirectory() as temp_dir:
        with zipfile.ZipFile(path) as archive:
            archive.extractall(temp_dir)
        matches = glob.glob(f"{temp_dir}/*/data/feature-table.biom")
        if not matches:
            sys.exit(f"no feature-table.biom inside {path}")
        return biom.load_table(matches[0]).to_dataframe(dense=True)


def asvs_to_strip():
    """ASVs whose taxon is a pathogen, an endosymbiont, or not bacterial.

    Pathogens go because the grouping is partly derived from them. Endosymbionts go because
    Rickettsia alone is ~65% of these communities and would dominate any compositional
    method; its variation is not what the question is about. Non-bacterial ASVs go because
    the 515F/926R primers co-amplify 18S and eukaryotic sequence is not part of a bacterial
    community comparison.
    """
    if not os.path.exists(DECISIONS):
        sys.exit(f"{DECISIONS} not found. Run make_taxon_decision_table.py first.")

    strip_categories = {"pathogen_bacterial", "pathogen_protozoan", "endosymbiont",
                        "ked_associated_bacterium", "offpanel_pathogen_genus",
                        "piroplasm_nonhuman", "crossamp", "mitochondrial",
                        "tick_host", "vertebrate_host", "human_dna",
                        "eukaryotic_in_16S"}
    strip_lineages = set()
    with open(DECISIONS, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("marker") == "16S" and row.get("category") in strip_categories:
                strip_lineages.add(row["full_lineage"])

    # Belt and braces. The category strip above only works if the decision table has been
    # regenerated since eukaryotic_in_16S was added. These lineage patterns catch the same
    # sequence directly, so the script is correct against a stale table too.
    NON_BACTERIAL = ("f__Mitochondria", "f__Chloroplast", ";Mitochondria", ";Chloroplast",
                     "d__Eukaryota", "d__Archaea")

    strip_ids = set()
    kept_lineage_by_id = {}
    with open(TAXONOMY, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            lineage = row[1].strip()
            if lineage in strip_lineages or any(m in lineage for m in NON_BACTERIAL):
                strip_ids.add(row[0].strip())
            else:
                kept_lineage_by_id[row[0].strip()] = row[1].strip()
    return strip_ids, kept_lineage_by_id


def genus_of(lineage):
    for part in reversed([p.strip() for p in lineage.split(";")]):
        name = part.split("__")[-1].strip()
        if name and name.lower() not in ("", "uncultured", "unassigned", "unclassified"):
            return name
    return "Unassigned"


def bray_curtis(matrix):
    """Bray-Curtis distances from a samples x features array of relative abundances."""
    n = matrix.shape[0]
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            numerator = np.abs(matrix[i] - matrix[j]).sum()
            denominator = (matrix[i] + matrix[j]).sum()
            value = numerator / denominator if denominator else 0.0
            distances[i, j] = distances[j, i] = value
    return distances


def permanova(distances, labels, permutations=9999, seed=1):
    """One-way PERMANOVA. Returns (pseudo_F, R2, p).

    Written out rather than imported so the calculation is visible: this is the number the
    whole analysis turns on and it should not sit behind a library call nobody checks.
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    n = len(labels)
    squared = distances ** 2
    total_ss = squared.sum() / (2 * n)

    def within_ss(assignment):
        total = 0.0
        for group in np.unique(assignment):
            index = np.where(assignment == group)[0]
            if len(index) < 2:
                continue
            block = squared[np.ix_(index, index)]
            total += block.sum() / (2 * len(index))
        return total

    n_groups = len(np.unique(labels))
    observed_within = within_ss(labels)
    observed_between = total_ss - observed_within
    f_observed = ((observed_between / (n_groups - 1)) /
                  (observed_within / (n - n_groups)))
    r2 = observed_between / total_ss

    count = 0
    for _ in range(permutations):
        shuffled = rng.permutation(labels)
        within = within_ss(shuffled)
        between = total_ss - within
        f = (between / (n_groups - 1)) / (within / (n - n_groups))
        if f >= f_observed:
            count += 1
    return f_observed, r2, (count + 1) / (permutations + 1)


def permdisp(distances, labels, permutations=9999, seed=1):
    """Distance-to-group-centroid homogeneity test. Returns (F, p, mean distance per group).

    This is the check that a significant PERMANOVA is a location difference and not unequal
    spread. With groups of 34/21/14 it is not optional.
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    squared = distances ** 2
    n = len(labels)

    # classical MDS to get coordinates, then Euclidean distance to each group centroid
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ squared @ centering
    values, vectors = np.linalg.eigh(gram)
    keep = values > 1e-9
    coordinates = vectors[:, keep] * np.sqrt(values[keep])

    def spreads(assignment):
        out = np.zeros(n)
        for group in np.unique(assignment):
            index = np.where(assignment == group)[0]
            centroid = coordinates[index].mean(axis=0)
            out[index] = np.sqrt(((coordinates[index] - centroid) ** 2).sum(axis=1))
        return out

    observed = spreads(labels)
    groups = np.unique(labels)
    grand_mean = observed.mean()
    between = sum(len(observed[labels == g]) *
                  (observed[labels == g].mean() - grand_mean) ** 2 for g in groups)
    within = sum(((observed[labels == g] - observed[labels == g].mean()) ** 2).sum()
                 for g in groups)
    f_observed = (between / (len(groups) - 1)) / (within / (n - len(groups)))

    count = 0
    for _ in range(permutations):
        shuffled = rng.permutation(labels)
        values_ = spreads(shuffled)
        gm = values_.mean()
        b = sum(len(values_[shuffled == g]) * (values_[shuffled == g].mean() - gm) ** 2
                for g in groups)
        w = sum(((values_[shuffled == g] - values_[shuffled == g].mean()) ** 2).sum()
                for g in groups)
        f = (b / (len(groups) - 1)) / (w / (n - len(groups)))
        if f >= f_observed:
            count += 1

    per_group = {g: float(observed[labels == g].mean()) for g in groups}
    return f_observed, (count + 1) / (permutations + 1), per_group


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--design", default="results/analysis/carriage_design_check.tsv")
    parser.add_argument("--outdir", default="results/analysis/carriage_microbiome")
    parser.add_argument("--permutations", type=int, default=9999)
    args = parser.parse_args()

    if not os.path.isfile(args.design):
        sys.exit(f"no such file: {args.design}. Run carriage_design_check.py first.")
    os.makedirs(args.outdir, exist_ok=True)

    design = pd.read_csv(args.design, sep="\t")
    deer = design[design["host"] == "Deer"].copy()
    say(f"design file: {len(design)} ticks; deer-collected: {len(deer)}")
    if len(deer) < 20:
        sys.exit("fewer than 20 deer ticks; not worth testing")

    say("\ngroups (deer only):")
    for group in GROUP_ORDER:
        say(f"  {group:<20} {int((deer['group'] == group).sum()):>3}")

    # ---- community table ----------------------------------------------------
    say(f"\nreading {TABLE_QZA}")
    table = load_qza_table(TABLE_QZA)
    say(f"  {table.shape[0]:,} ASVs x {table.shape[1]:,} libraries")

    strip_ids, kept_lineage = asvs_to_strip()
    say(f"  stripping {len(strip_ids):,} ASVs (pathogens, endosymbionts, non-bacterial)")
    table = table.loc[[i for i in table.index if i not in strip_ids]]
    say(f"  {table.shape[0]:,} ASVs remain")

    # libraries -> specimen barcode
    library_to_barcode = {}
    for library in table.columns:
        library_to_barcode[library] = str(library).split(".")[0]

    wanted = {str(b): b for b in deer["barcode"]}
    columns = [c for c in table.columns if library_to_barcode[c] in wanted]
    if not columns:
        sys.exit("no library in the feature table matched a deer barcode; "
                 f"table example: {list(table.columns)[:2]}; "
                 f"design example: {list(wanted)[:2]}")

    # sum a specimen's libraries, then convert to relative abundance
    summed = {}
    for column in columns:
        barcode = library_to_barcode[column]
        summed[barcode] = summed.get(barcode, 0) + table[column]
    community = pd.DataFrame(summed).T
    say(f"\nspecimens with community data: {len(community)} of {len(deer)}")

    depth = community.sum(axis=1)
    say(f"  read depth: min {int(depth.min()):,}  median {int(depth.median()):,}  "
        f"max {int(depth.max()):,}")
    low = depth[depth < 1000]
    if len(low):
        say(f"  NOTE: {len(low)} specimen(s) below 1,000 reads; retained but they will "
            f"sit apart in ordination space")

    labels = deer.set_index(deer["barcode"].astype(str))["group"].reindex(community.index)
    keep = labels.notna() & (depth > 0)
    community, labels = community[keep], labels[keep]
    say(f"  testable: {len(community)}")
    say("  " + "  ".join(f"{g}={int((labels == g).sum())}" for g in GROUP_ORDER))

    if len(community) < 20 or labels.nunique() < 2:
        sys.exit("too few specimens or groups to test")

    relative = community.div(community.sum(axis=1), axis=0).fillna(0).values

    # ---- alpha --------------------------------------------------------------
    say("\n=== alpha diversity ===")
    richness = (community > 0).sum(axis=1)
    proportions = np.where(relative > 0, relative, 1)
    shannon = -(relative * np.log(proportions)).sum(axis=1)
    alpha = pd.DataFrame({"group": labels.values, "observed_asvs": richness.values,
                          "shannon": shannon}, index=community.index)
    say(alpha.groupby("group")[["observed_asvs", "shannon"]]
             .agg(["median", "mean", "count"]).round(2).to_string())

    try:
        from scipy.stats import kruskal
        for metric in ("observed_asvs", "shannon"):
            samples = [alpha[alpha["group"] == g][metric].values
                       for g in GROUP_ORDER if (alpha["group"] == g).any()]
            statistic, p = kruskal(*samples)
            say(f"  Kruskal-Wallis {metric}: H = {statistic:.2f}, p = {p:.4f}")
    except ImportError:
        say("  scipy unavailable; alpha tests skipped")
    alpha.to_csv(f"{args.outdir}/alpha_by_group.tsv", sep="\t")

    # ---- beta ---------------------------------------------------------------
    say("\n=== beta diversity (Bray-Curtis) ===")
    distances = bray_curtis(relative)
    f_stat, r2, p_value = permanova(distances, labels.values, args.permutations)
    say(f"  PERMANOVA: pseudo-F = {f_stat:.3f}, R2 = {100*r2:.1f}%, p = {p_value:.4f}")
    say(f"  for comparison, BATCH explains {BATCH_R2_PCT}% (p = {BATCH_P}) in the deer set")
    if 100 * r2 < BATCH_R2_PCT:
        say(f"  -> carriage explains LESS than batch. Say so plainly.")

    d_stat, d_p, spread = permdisp(distances, labels.values, args.permutations)
    say(f"  PERMDISP:  F = {d_stat:.3f}, p = {d_p:.4f}")
    for group in GROUP_ORDER:
        if group in spread:
            say(f"    mean distance to centroid, {group:<20} {spread[group]:.3f}")
    if d_p < 0.05:
        say("  ** PERMDISP SIGNIFICANT: group dispersions differ. A significant PERMANOVA")
        say("     here may reflect unequal spread rather than a difference in composition. **")

    pd.DataFrame([{"test": "PERMANOVA", "statistic": f_stat, "R2_pct": 100 * r2,
                   "p": p_value},
                  {"test": "PERMDISP", "statistic": d_stat, "R2_pct": "", "p": d_p}]
                 ).to_csv(f"{args.outdir}/permanova_permdisp.tsv", sep="\t", index=False)

    # ---- which genera, descriptively ---------------------------------------
    say("\n=== genus-level mean relative abundance by group ===")
    say("(descriptive; no differential-abundance test is run here)")
    genus_names = [genus_of(kept_lineage.get(a, "")) for a in community.columns]
    by_genus = pd.DataFrame(relative, index=community.index, columns=genus_names)
    by_genus = by_genus.T.groupby(level=0).sum().T
    by_genus["group"] = labels.values
    means = by_genus.groupby("group").mean().T
    means["spread"] = means.max(axis=1) - means.min(axis=1)
    top = means.sort_values("spread", ascending=False).head(15)
    say((100 * top[[g for g in GROUP_ORDER if g in top.columns]]).round(2).to_string())
    means.to_csv(f"{args.outdir}/genus_by_group.tsv", sep="\t")

    # ---- what this can and cannot say --------------------------------------
    say("\n=== READ THIS BEFORE QUOTING ANY OF THE ABOVE ===")
    say("  - Deer-collected ticks only, treated as Ixodes because 0 of 273 deer specimens")
    say("    carried Dermacentor. This is an inference from the host, not a per-specimen")
    say("    18S call: only 33 of 83 ticks have one and 42 have no 18S library at all.")
    say("  - Pathogen and endosymbiont genera were REMOVED before testing. Without that,")
    say("    the comparison is partly circular, since the grouping derives from pathogens.")
    say(f"  - Batch is NOT in the model. Only 27 of 69 deer ticks carry a batch label, and")
    say(f"    batch explains {BATCH_R2_PCT}% of community variation with PERMDISP also")
    say("    significant. A carriage effect of comparable or smaller size cannot be")
    say("    distinguished from a processing-run difference in this design.")
    say("  - Groups are unbalanced (34/21/14). PERMDISP is reported for that reason.")
    say("  - A null result is a legitimate answer and should be reported as one.")

    report_path = f"{args.outdir}/carriage_microbiome_report.txt"
    with open(report_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"\nreport -> {report_path}")


if __name__ == "__main__":
    main()
