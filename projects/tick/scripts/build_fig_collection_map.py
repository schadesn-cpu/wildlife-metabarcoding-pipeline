#!/usr/bin/env python3
"""build_fig_collection_map.py

I build the collection map: animals sampled per county across Vermont, New
Hampshire and Maine.

WHY COUNTY
----------
These are hunter-harvest records. A town-level map publishes where individual
animals were taken. County is coarse enough that no harvest location is
recoverable and fine enough to show that sampling followed the trailhead
check stations rather than the geography.

WHY A COUNT AND NOT A RATE
--------------------------
The number shaded is animals sampled, which is a property of our effort and not
of the landscape. Nothing here is a prevalence, a density or a risk, and the
legend says animals so that it cannot be read as one.

The count is also printed inside every county that has animals, because a
choropleth over counties of wildly different size invites the eye to read area
as quantity. Aroostook is larger than Vermont and holds almost nothing.

COLOUR
------
One hue, light to dark, which is what a sequential quantity takes. Counties
with no animals are a neutral grey and are drawn, not omitted, so a reader can
see where we did not sample. Classed rather than continuous, because five
labelled bins are read more accurately than a gradient bar.

WHAT I READ
-----------
    results/metadata/crosswalk_4marker.tsv        the animal set
    metadata/tick_full_metadata.xlsx, 'Metadata'  barcode, town, state
    results/metadata/town_aliases.tsv             confirmed spelling fixes
    metadata/tl_2023_{50,33,23}_cousub.zip        town to county
    metadata/tl_2023_us_county.zip                county outlines

WHAT I WRITE
------------
    results/figures/fig1_collection_map.png
    results/figures/fig1_collection_map.pdf
    results/analysis/animals_by_county.tsv

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_fig_collection_map.py
    python3 scripts/build_fig_collection_map.py --write
"""
import io
import math
import os
import sys
import zipfile

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

COUNTY = "metadata/tl_2023_us_county.zip"
ALIASES = "results/metadata/town_aliases.tsv"
OUT = "results/figures/fig1_collection_map"
TABLE = "results/analysis/animals_by_county.tsv"
BINS = [(1, 2), (3, 9), (10, 24), (25, 49), (50, None)]
NOT_RECORDED = ("not recorded", "unrecorded", "unknown", "none", "na", "n a")
RAMP = ["#d3e3f7", "#a8c8ee", "#6ea5e2", "#3a83d8", "#1e5aa8"]
EMPTY = "#eceae5"
EDGE = "#ffffff"
INK = "#0b0b0b"
MUTED = "#52514e"
SURFACE = "#fcfcfb"
PALE = "#c9c6bf"


def resolver():
    """I return check_town_to_county, whose loading I reuse whole."""
    try:
        import check_town_to_county as module
    except ImportError as problem:
        sys.exit("cannot import check_town_to_county (%s). It must sit in "
                 "scripts/ beside me." % problem)
    return module


def read_aliases():
    """I return {(state, tidied recorded name): census name}."""
    if not os.path.isfile(ALIASES):
        sys.exit("missing %s. Run propose_town_aliases.py and confirm the "
                 "pairs first." % ALIASES)
    frame = pd.read_csv(ALIASES, sep="\t", dtype=str)
    for column in ("state", "recorded", "census"):
        if column not in frame.columns:
            sys.exit("%s has no '%s' column; it has %s"
                     % (ALIASES, column, ", ".join(map(str, frame.columns))))
    return frame


def read_county_shapes(shapefile, states):
    """I yield (statefp, countyfp, name, [parts]) for the states I am given."""
    if not os.path.isfile(COUNTY):
        sys.exit("missing %s. Download the Census county file first." % COUNTY)
    with zipfile.ZipFile(COUNTY) as archive:
        names = {os.path.splitext(n)[1].lower(): n for n in archive.namelist()}
        for needed in (".shp", ".dbf"):
            if needed not in names:
                sys.exit("%s has no %s inside it." % (COUNTY, needed))
        reader = shapefile.Reader(
            shp=io.BytesIO(archive.read(names[".shp"])),
            dbf=io.BytesIO(archive.read(names[".dbf"])))
        fields = [f[0] for f in reader.fields[1:]]
        for wanted in ("STATEFP", "COUNTYFP", "NAME"):
            if wanted not in fields:
                sys.exit("%s has no %s field; it has %s"
                         % (COUNTY, wanted, ", ".join(fields)))
        for item in reader.iterShapeRecords():
            row = dict(zip(fields, item.record))
            if str(row["STATEFP"]) not in states:
                continue
            points = list(item.shape.points)
            breaks = list(item.shape.parts) + [len(points)]
            parts = [points[breaks[i]:breaks[i + 1]]
                     for i in range(len(breaks) - 1)]
            yield (str(row["STATEFP"]), str(row["COUNTYFP"]),
                   str(row["NAME"]), [p for p in parts if len(p) > 2])


def centroid(part):
    """I return the area-weighted centre of one ring, by the shoelace rule.

    The mean of the points would sit wherever the vertices happen to be dense,
    which on a coastline is out at sea.
    """
    total = across = down = 0.0
    for i in range(len(part) - 1):
        x1, y1 = part[i]
        x2, y2 = part[i + 1]
        step = x1 * y2 - x2 * y1
        total += step
        across += (x1 + x2) * step
        down += (y1 + y2) * step
    if abs(total) < 1e-12:
        xs = [p[0] for p in part]
        ys = [p[1] for p in part]
        return sum(xs) / len(xs), sum(ys) / len(ys)
    return across / (3 * total), down / (3 * total)


def band_of(count):
    """I return the index of the bin a count falls in, or None for zero."""
    if count < 1:
        return None
    for index, (low, high) in enumerate(BINS):
        if count >= low and (high is None or count <= high):
            return index
    return len(BINS) - 1


def band_label(low, high):
    return "%d or more" % low if high is None else "%d to %d" % (low, high)


def count_animals(module, shapefile):
    """I return {(statefp, countyfp): animals} and a tally of what was lost."""
    aliases = read_aliases()
    _, merged = module.read_animals()
    lookup, code_of = {}, {}
    for code, abbreviation in module.STATES.items():
        lookup[abbreviation] = module.read_towns(shapefile, code)
        code_of[abbreviation] = code

    replace = {}
    for state, recorded, census in zip(aliases["state"], aliases["recorded"],
                                       aliases["census"]):
        replace[(str(state).strip().upper(),
                 module.tidy(recorded))] = module.tidy(census)

    counts, lost, unplaced = {}, {"no town": 0, "no state": 0}, {}
    for state, town in zip(merged["state"], merged["town"]):
        if pd.isna(town) or not str(town).strip():
            lost["no town"] += 1
            continue
        if module.tidy(town) in NOT_RECORDED:
            # A town field saying 'not recorded' is a missing town, not a
            # town I failed to find, and an alias cannot fix it because an
            # alias has to point at a real Census place.
            lost["no town"] += 1
            continue
        abbreviation = str(state).strip().upper()
        if abbreviation not in lookup:
            lost["no state"] += 1
            continue
        key = module.tidy(town)
        key = replace.get((abbreviation, key), key)
        hit = lookup[abbreviation].get(key)
        if hit is None:
            unplaced[(abbreviation, str(town).strip())] = \
                unplaced.get((abbreviation, str(town).strip()), 0) + 1
            continue
        where = (code_of[abbreviation], hit[0])
        counts[where] = counts.get(where, 0) + 1
    return counts, lost, unplaced


def main():
    write = "--write" in sys.argv[1:]
    module = resolver()
    shapefile = module.shapefile_module()
    counts, lost, unplaced = count_animals(module, shapefile)

    placed = sum(counts.values())
    print("animals placed in a county: %d" % placed)
    print("  no town recorded:         %d" % lost["no town"])
    print("  state outside VT/NH/ME:   %d" % lost["no state"])
    if unplaced:
        print()
        print("  %d towns still do not resolve, covering %d animals:"
              % (len(unplaced), sum(unplaced.values())))
        for (state, town), count in sorted(unplaced.items(),
                                           key=lambda kv: -kv[1]):
            print("    %-5s %-28s %3d" % (state, town[:28], count))
        print()
        sys.exit("Add these to %s before I draw anything. A map that silently "
                 "drops animals is worse than no map." % ALIASES)
    print()

    states = set(module.STATES)
    by_code = {code: abbreviation for code, abbreviation
               in module.STATES.items()}
    shapes = list(read_county_shapes(shapefile, set(by_code)))
    if not shapes:
        sys.exit("no counties found for %s in %s"
                 % (", ".join(sorted(states)), COUNTY))
    print("counties drawn: %d across %d states"
          % (len(shapes), len({s for s, _, _, _ in shapes})))

    figure, axis = plt.subplots(figsize=(6.6, 6.2))
    figure.patch.set_facecolor(SURFACE)
    axis.set_facecolor(SURFACE)

    rows, extent = [], [180.0, 90.0, -180.0, -90.0]
    reach = {}
    for statefp, countyfp, name, parts in shapes:
        count = counts.get((statefp, countyfp), 0)
        band = band_of(count)
        colour = EMPTY if band is None else RAMP[band]
        biggest, best = None, -1.0
        for part in parts:
            xs = [p[0] for p in part]
            ys = [p[1] for p in part]
            extent[0] = min(extent[0], min(xs))
            extent[1] = min(extent[1], min(ys))
            extent[2] = max(extent[2], max(xs))
            extent[3] = max(extent[3], max(ys))
            here = reach.setdefault(statefp, [180.0, -180.0, -90.0])
            here[0] = min(here[0], min(xs))
            here[1] = max(here[1], max(xs))
            here[2] = max(here[2], max(ys))
            span = (max(xs) - min(xs)) * (max(ys) - min(ys))
            if span > best:
                best, biggest = span, part
            axis.add_patch(Polygon(part, closed=True, facecolor=colour,
                                   edgecolor=EDGE, linewidth=0.6, zorder=2))
        rows.append({"state": by_code[statefp], "county": name,
                     "county_fips": "%s%s" % (statefp, countyfp),
                     "animals": count})
        if count:
            x, y = centroid(biggest)
            axis.text(x, y, "%d" % count, ha="center", va="center",
                      fontsize=7.4, zorder=4,
                      color="white" if band >= 3 else INK)

    # A state name sits above that state's northern edge, which for all three
    # is open space, rather than on top of its own counties.
    for statefp, (west, east, north) in reach.items():
        axis.text((west + east) / 2.0, north + 0.10, by_code[statefp],
                  ha="center", va="bottom", fontsize=10.5, color=MUTED,
                  zorder=5)

    margin = 0.25
    axis.set_xlim(extent[0] - margin, extent[2] + margin)
    axis.set_ylim(extent[1] - margin, extent[3] + margin * 3)
    middle = math.radians((extent[1] + extent[3]) / 2.0)
    axis.set_aspect(1.0 / max(math.cos(middle), 0.1))
    axis.axis("off")

    patches = [Polygon([(0, 0)], facecolor=EMPTY, edgecolor=PALE,
                       linewidth=0.6)]
    labels = ["none"]
    for index, (low, high) in enumerate(BINS):
        patches.append(Polygon([(0, 0)], facecolor=RAMP[index],
                               edgecolor=EDGE, linewidth=0.6))
        labels.append(band_label(low, high))
    # The bounding box of these three states has open Atlantic in its
    # south-east corner, south of the Maine coast, which is the largest area
    # certain to hold no county. The north-west corner is Quebec but northern
    # Vermont reaches into it.
    legend = axis.legend(patches, labels, loc="lower right", frameon=False,
                         fontsize=8, title="Animals sampled",
                         title_fontsize=8.5, borderpad=0.2, labelspacing=0.5,
                         handlelength=1.3)
    legend.get_title().set_ha("left")
    figure.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)

    table = pd.DataFrame(rows).sort_values(["animals", "state", "county"],
                                           ascending=[False, True, True])
    with_animals = table[table["animals"] > 0]
    print("counties holding animals: %d of %d"
          % (len(with_animals), len(table)))
    print("busiest: %s"
          % ", ".join("%s %s %d" % (r["state"], r["county"], r["animals"])
                      for _, r in with_animals.head(5).iterrows()))
    if not write:
        print()
        print("Nothing written. Rerun with --write to save %s.png, .pdf and %s"
              % (OUT, TABLE))
        return
    for path in ("%s.png" % OUT, "%s.pdf" % OUT, TABLE):
        if os.path.exists(path):
            sys.exit("%s already exists. Move it aside first." % path)
    figure.savefig(OUT + ".png", dpi=300, facecolor=figure.get_facecolor())
    figure.savefig(OUT + ".pdf", facecolor=figure.get_facecolor())
    table.to_csv(TABLE, sep="\t", index=False)
    print("wrote %s.png, %s.pdf and %s" % (OUT, OUT, TABLE))


if __name__ == "__main__":
    main()
