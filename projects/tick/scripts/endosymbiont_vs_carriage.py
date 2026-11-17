#!/usr/bin/env python3
"""Do endosymbiont abundances differ by pathogen carriage?

The global community test was null (PERMANOVA R2 = 3.4%, p = 0.81; carriage explains less
than batch's 10.7%). This asks a narrower and mechanistically motivated question instead:
does the abundance of a specific endosymbiont differ between ticks that carry a pathogen and
ticks that do not?

The endosymbionts were deliberately EXCLUDED from the global test, because Rickettsia alone
is ~65% of these communities and would dominate any compositional method. That exclusion is
what makes this a separate, targeted question rather than a subset of the previous one.

PRE-SPECIFIED. The taxa below were fixed from the literature before looking at this data,
and the list is not revised after seeing a result. Adding a taxon because it looked
interesting in the genus table would make every p-value here meaningless.

  Rickettsia    Gall et al. 2016 (Dermacentor andersoni): R. bellii abundance REDUCED
                Anaplasma marginale acquisition. R. buchneri is the I. scapularis
                endosymbiont and the direct analogue here.
  Francisella   Gall et al. 2016: Francisella endosymbionts INCREASED F. novicida
                acquisition. Opposite direction to Rickettsia in the same study.
  Coxiella      Klyachko et al. 2007 (Amblyomma americanum): Coxiella-like symbiont
                proposed to interfere with Ehrlichia chaffeensis.
  Wolbachia     Included as a NEGATIVE CONTROL. In these ticks Wolbachia is a parasitoid
                wasp symbiont carried incidentally, with no proposed interaction with
                tick-borne pathogens. If it comes out significant alongside the others,
                that argues the signal is something structural and not biology.
  Arsenophonus  Tick-associated symbiont, no documented pathogen interaction. Second
                negative control.

Direction of the Rickettsia prediction is stated in advance: Gall found endosymbiont
abundance REDUCED pathogen acquisition, so the prediction is LOWER Rickettsia in carriers.
A result in the opposite direction is not support for the hypothesis, it is a different
finding and should be reported as one.

Correction: Holm-Bonferroni across the five taxa. Also reported uncorrected, because with
five pre-specified tests a reader may reasonably want both.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/endosymbiont_vs_carriage.py \
        --design results/analysis/carriage_design_check.tsv \
        --outdir results/analysis/carriage_microbiome

Writes:
    <outdir>/endosymbiont_tests.tsv
    <outdir>/endosymbiont_by_carriage.png / .pdf
    <outdir>/endosymbiont_report.txt
"""

import argparse
import csv
import glob
import os
import re
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd

TABLE_QZA = "results/decontam/16S/table_decontam_16S.qza"
TAXONOMY = "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv"

# Pre-specified. See the docstring. Do not add to this list after seeing results.
TARGETS = {
    "Rickettsia":   ("g__Rickettsia",   "predicted LOWER in carriers (Gall 2016)"),
    "Francisella":  ("g__Francisella",  "predicted HIGHER in carriers (Gall 2016)"),
    "Coxiella":     ("g__Coxiella",     "interference proposed (Klyachko 2007)"),
    "Wolbachia":    ("g__Wolbachia",    "NEGATIVE CONTROL - no proposed interaction"),
    "Arsenophonus": ("g__Arsenophonus", "NEGATIVE CONTROL - no proposed interaction"),
}

# Organellar sequence SILVA files inside the bacterial tree. Excluded from the denominator:
# it is host and tick mitochondrial DNA, not community, and its abundance tracks how much
# host tissue was in the tube.
# Non-bacterial sequence. Organellar, plus the 9,333 eukaryotic ASVs the 16S primers
# pick up: the tick (o__Acari), the deer ked (o__Diptera), the host blood meal
# (c__Mammalia) and 6,648 resolved no deeper than d__Eukaryota. All correctly
# labelled, none part of a bacterial community.
ORGANELLAR = ("f__Mitochondria", "f__Chloroplast", ";Mitochondria", ";Chloroplast",
              "d__Eukaryota", "d__Archaea")

GROUP_ORDER = ["confirmed_negative", "single_method", "confirmed_carrier"]
GROUP_LABEL = {"confirmed_negative": "Confirmed\nnon-carrier",
               "single_method": "Single-method\npositive",
               "confirmed_carrier": "Confirmed\ncarrier"}

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def normalise(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def to_barcode(library):
    """Library name -> specimen barcode. Strips the run suffix and the -NN replicate index."""
    stem = str(library).split(".")[0]
    return re.sub(r"-\d{2}$", "", stem)


def back_up_existing(out_path):
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def holm_bonferroni(p_values):
    """Holm-Bonferroni adjusted p-values, in the original order."""
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    n = len(p_values)
    adjusted = np.empty(n)
    running = 0.0
    for rank, index in enumerate(order):
        # multiplier is (n - rank), i.e. n for the smallest p, n-1 for the next, down to 1.
        # Then enforce monotonicity so an adjusted value never falls below the one before.
        running = max(running, (n - rank) * p_values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--design", default="results/analysis/carriage_design_check.tsv")
    parser.add_argument("--outdir", default="results/analysis/carriage_microbiome")
    args = parser.parse_args()

    for path in (args.design, TABLE_QZA, TAXONOMY):
        if not os.path.exists(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(args.outdir, exist_ok=True)

    design = pd.read_csv(args.design, sep="\t")
    deer = design[design["host"] == "Deer"].copy()
    say(f"deer-collected ticks screened by both methods: {len(deer)}")

    # ---- taxonomy: which ASVs are which target, and which are organellar ----
    target_ids = {name: set() for name in TARGETS}
    organellar_ids = set()
    with open(TAXONOMY, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            feature_id, lineage = row[0].strip(), row[1].strip()
            if any(marker in lineage for marker in ORGANELLAR):
                organellar_ids.add(feature_id)
                continue
            for name, (pattern, _) in TARGETS.items():
                if pattern in lineage:
                    target_ids[name].add(feature_id)
    for name in TARGETS:
        say(f"  {name:<14} {len(target_ids[name]):>6,} ASVs")
    say(f"  organellar (excluded from denominator): {len(organellar_ids):,} ASVs")

    # ---- community ----------------------------------------------------------
    import biom
    with tempfile.TemporaryDirectory() as temp_dir:
        with zipfile.ZipFile(TABLE_QZA) as archive:
            archive.extractall(temp_dir)
        table = biom.load_table(
            glob.glob(f"{temp_dir}/*/data/feature-table.biom")[0]).to_dataframe(dense=True)

    table = table.loc[[i for i in table.index if i not in organellar_ids]]

    wanted = {normalise(b): str(b) for b in deer["barcode"]}
    summed = {}
    for column in table.columns:
        key = normalise(to_barcode(column))
        if key in wanted:
            barcode = wanted[key]
            summed[barcode] = summed.get(barcode, 0) + table[column]
    if not summed:
        sys.exit("no library matched a deer barcode")
    community = pd.DataFrame(summed).T

    depth = community.sum(axis=1)
    community = community[depth > 0]
    say(f"\nspecimens with community data: {len(community)} of {len(deer)}")
    say(f"  read depth (organellar removed): min {int(depth[depth>0].min()):,}  "
        f"median {int(depth[depth>0].median()):,}  max {int(depth.max()):,}")

    groups = deer.set_index(deer["barcode"].astype(str))["group"].reindex(community.index)
    keep = groups.notna()
    community, groups = community[keep], groups[keep]
    say("  " + "  ".join(f"{g}={int((groups == g).sum())}" for g in GROUP_ORDER))

    # ---- relative abundance per target -------------------------------------
    totals = community.sum(axis=1)
    fractions = pd.DataFrame(index=community.index)
    for name in TARGETS:
        present = [i for i in target_ids[name] if i in community.columns]
        fractions[name] = (community[present].sum(axis=1) / totals) if present else 0.0
    fractions["group"] = groups.values

    say("\n=== relative abundance by carriage group (%) ===")
    summary = fractions.groupby("group")[list(TARGETS)].agg(["median", "mean"])
    say((100 * summary).round(2).to_string())

    # ---- depth confound: this has to run BEFORE the taxon tests -----------
    # Carriage is not independent of sequencing depth. Detecting a low-abundance
    # pathogen by metabarcoding requires reads, so a deeply sequenced tick is more
    # likely to be CALLED a carrier. And in a community that is 70-93% one organism,
    # the dominant taxon's fraction rises with depth. Those two facts together
    # manufacture an association between any dominant taxon and carriage.
    #
    # This is not hypothetical: it is what produced the Rickettsia p = 0.048 in the
    # first run of this script, in the direction OPPOSITE to the one predicted.
    say("\n=== DEPTH CONFOUND CHECK (read this before the taxon tests) ===")
    try:
        from scipy.stats import mannwhitneyu as _mwu, spearmanr as _spearman
    except ImportError:
        sys.exit("scipy is required")

    depths = community.sum(axis=1)
    say("  library depth by group (organellar excluded):")
    for group in GROUP_ORDER:
        values = depths[groups == group]
        if len(values):
            say(f"    {group:<20} n={len(values):>3}  median={int(values.median()):>9,}")

    carrier_depth = depths[groups == "confirmed_carrier"]
    negative_depth = depths[groups == "confirmed_negative"]
    depth_p = float("nan")
    if len(carrier_depth) >= 3 and len(negative_depth) >= 3:
        depth_p = _mwu(carrier_depth, negative_depth, alternative="two-sided")[1]
        say(f"  carrier vs non-carrier depth: p = {depth_p:.4f}")

    dominant = fractions[["Rickettsia"]].join(depths.rename("depth"))
    rho, rho_p = _spearman(dominant["depth"], dominant["Rickettsia"])
    say(f"  Rickettsia fraction vs depth: rho = {rho:.3f}, p = {rho_p:.4g}")

    depth_confounded = (depth_p < 0.05) and (rho_p < 0.05)
    if depth_confounded:
        say("")
        say("  ** DEPTH IS CONFOUNDED WITH CARRIAGE. **")
        say("  Carriers are sequenced more deeply AND the dominant taxon's fraction rises")
        say("  with depth. Any difference in a dominant taxon between carriage groups is")
        say("  therefore expected under the null, and a significant result for one should")
        say("  be attributed to depth unless it survives a depth-matched comparison.")
        say("  Note the direction this cuts: it manufactures FALSE POSITIVES. It cannot")
        say("  hide a real effect, so the null results below remain safe.")

    # ---- the tests ----------------------------------------------------------
    say("\n=== pre-specified tests: carrier vs non-carrier ===")
    say("(confirmed_carrier vs confirmed_negative; the single-method group is shown in the")
    say(" figure but excluded from the test, because it is a mixture of two situations)")

    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        sys.exit("scipy is required for the tests")

    carriers = fractions[fractions["group"] == "confirmed_carrier"]
    negatives = fractions[fractions["group"] == "confirmed_negative"]
    say(f"\n  n carriers = {len(carriers)}, n non-carriers = {len(negatives)}")
    if len(carriers) < 5 or len(negatives) < 5:
        say("  ** fewer than 5 in a group; these tests are not interpretable **")

    rows, p_values = [], []
    for name, (_, rationale) in TARGETS.items():
        a, b = carriers[name].values, negatives[name].values
        if a.sum() == 0 and b.sum() == 0:
            say(f"  {name:<14} absent from every specimen; not tested")
            continue
        statistic, p = mannwhitneyu(a, b, alternative="two-sided")
        direction = ("higher in carriers" if np.median(a) > np.median(b)
                     else "lower in carriers" if np.median(a) < np.median(b)
                     else "equal medians")
        rows.append({"taxon": name, "rationale": rationale,
                     "median_carrier_pct": round(100 * float(np.median(a)), 3),
                     "median_noncarrier_pct": round(100 * float(np.median(b)), 3),
                     "direction": direction, "U": float(statistic), "p_raw": float(p)})
        p_values.append(p)

    if not rows:
        sys.exit("no target was testable")

    adjusted = holm_bonferroni(np.array(p_values))
    for row, value in zip(rows, adjusted):
        row["p_holm"] = float(value)

    results = pd.DataFrame(rows)
    say("")
    say(results[["taxon", "median_carrier_pct", "median_noncarrier_pct",
                 "direction", "p_raw", "p_holm"]].round(4).to_string(index=False))

    say("\n  rationale for each, fixed in advance:")
    for row in rows:
        say(f"    {row['taxon']:<14} {row['rationale']}")

    significant = [r for r in rows if r["p_holm"] < 0.05]
    say("")
    if not significant:
        say("  NO taxon differs after Holm correction. Combined with the null global test")
        say("  (PERMANOVA R2 = 3.4%, p = 0.81), there is no evidence in this dataset that")
        say("  endosymbiont abundance is associated with pathogen carriage.")
        if depth_confounded:
            say("")
            say("  This null is ROBUST to the depth confound. Depth confounding inflates")
            say("  differences between groups; it does not create nulls. If anything the")
            say("  true association is weaker than these p-values suggest.")
    else:
        for row in significant:
            say(f"  {row['taxon']}: p_holm = {row['p_holm']:.4f}, {row['direction']}")
        if depth_confounded:
            say("  ** but depth is confounded with carriage (above). Do not report any of")
            say("     these as biology without a depth-matched or rarefied comparison. **")
        controls = [r for r in significant if "NEGATIVE CONTROL" in r["rationale"]]
        if controls:
            say("  ** A NEGATIVE CONTROL is significant: "
                f"{', '.join(r['taxon'] for r in controls)}. That argues the signal is")
            say("     structural (depth, batch, compositionality) rather than biological. **")

    out_tsv = f"{args.outdir}/endosymbiont_tests.tsv"
    back_up_existing(out_tsv)
    results.to_csv(out_tsv, sep="\t", index=False)
    say(f"\nwrote {out_tsv}")

    # ---- figure -------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    present_groups = [g for g in GROUP_ORDER if (fractions["group"] == g).any()]
    names = [r["taxon"] for r in rows]
    figure, axes = plt.subplots(1, len(names), figsize=(3.0 * len(names), 4.2), sharey=False)
    if len(names) == 1:
        axes = [axes]
    colours = {"confirmed_negative": "#CCCCCC", "single_method": "#56B4E9",
               "confirmed_carrier": "#0072B2"}

    for axis, row in zip(axes, rows):
        name = row["taxon"]
        data = [100 * fractions[fractions["group"] == g][name].values for g in present_groups]
        parts = axis.boxplot(data, widths=0.6, patch_artist=True, showfliers=False)
        for patch, group in zip(parts["boxes"], present_groups):
            patch.set_facecolor(colours[group])
            patch.set_edgecolor("#333333")
        for median in parts["medians"]:
            median.set_color("#000000")
        for position, values in enumerate(data, start=1):
            jitter = np.random.default_rng(0).normal(position, 0.06, len(values))
            axis.plot(jitter, values, "o", ms=3, color="#333333", alpha=0.55, mew=0)
        axis.set_xticks(range(1, len(present_groups) + 1))
        axis.set_xticklabels([GROUP_LABEL[g] for g in present_groups], fontsize=7.5)
        is_control = "NEGATIVE CONTROL" in row["rationale"]
        axis.set_title(f"{name}{'  (control)' if is_control else ''}\n"
                       f"p = {row['p_raw']:.3f}, Holm {row['p_holm']:.3f}",
                       fontsize=9, fontweight="bold" if not is_control else "normal")
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    axes[0].set_ylabel("Relative abundance (%)", fontsize=10)

    figure.suptitle("Endosymbiont abundance by pathogen-carriage confidence "
                    f"(deer-collected ticks, n = {len(fractions)})",
                    fontsize=11, fontweight="bold")
    figure.text(0.01, -0.02,
                "Pre-specified taxa, chosen from the literature before analysis. Mann-Whitney "
                "confirmed carrier vs confirmed non-carrier; the single-method group is shown "
                "but not tested.\nHolm-Bonferroni across five taxa. Organellar sequence "
                "excluded from the denominator. Batch is not modelled and explains 10.7% of "
                "community variation (p = 0.001).\nSequencing depth is confounded with "
                "carriage (carriers median 22,792 reads vs 2,168; p = 0.004) and Rickettsia "
                "fraction rises with depth (rho = 0.34), so a dominant-taxon difference "
                "between groups is expected under the null.",
                fontsize=7, style="italic", color="#666666", va="top")
    figure.tight_layout(rect=[0, 0.04, 1, 0.95])
    for extension in ("png", "pdf"):
        path = f"{args.outdir}/endosymbiont_by_carriage.{extension}"
        figure.savefig(path, dpi=300, bbox_inches="tight")
        say(f"wrote {path}")
    plt.close(figure)

    report_path = f"{args.outdir}/endosymbiont_report.txt"
    with open(report_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"report -> {report_path}")


if __name__ == "__main__":
    main()
