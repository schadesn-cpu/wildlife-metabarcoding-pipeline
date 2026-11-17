#!/usr/bin/env python3
"""extract_its_entomopathogens.py

Writes the six entomopathogen ITS sequences to a small fasta, so the detections
in the VT-12 update can be checked by someone other than us.

WHY
---
its_asv_table_FILTERED.tsv carries feature_id, taxon and read counts but no
sequences, and results/tables/ has asv_sequences_16S / 18S / cytb.fasta with no
ITS equivalent. So the Beauveria, Metarhizium, Cordyceps, Simplicillium and
Leptobacillium calls currently travel as names a reader has to take on trust.
Six short sequences let C. Sullivan, or a mycologist she asks, BLAST them
independently. Several rest on 2 to 19 reads, which is exactly when an
independent look is worth most.

HOW IT FINDS THINGS, rather than assuming
-----------------------------------------
  names + reads + library : parsed from the entomopathogen block of
                            its_build_report.txt
  name -> feature_id      : matched through its_asv_table_FILTERED.tsv on the
                            lowest_named column, with the read count used as a
                            cross-check
  sequences               : every .qza under qiime2/ITS1-2/ is opened as a zip
                            and searched for data/dna-sequences.fasta, so the
                            artifact can be named anything

It refuses rather than guesses: if a name matches no feature, or matches more
than one, or the read counts disagree between the two files, it says so and
writes nothing. If an id has no sequence in any artifact it is reported and
the others are still written.

OUTPUT
  results/tables/its_entomopathogen_seqs.fasta
Headers carry the taxon, read count and library so the file stands alone.
Will not overwrite an existing file.
"""
import glob
import os
import re
import sys
import zipfile

REPORT = "results/tables/its_build_report.txt"
TABLE = "results/tables/its_asv_table_FILTERED.tsv"
QIIME = "qiime2/ITS1-2"
OUT = "results/tables/its_entomopathogen_seqs.fasta"


def parse_report(path):
    """The block that begins '=== entomopathogenic genera', lines like
         Beauveria                            19 reads   BearA-01"""
    hits, inside = [], False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "entomopathogenic genera" in line:
                inside = True
                continue
            if inside:
                m = re.match(r'^\s{2,}(\S.*?)\s{2,}(\d[\d,]*)\s+reads?\s+(\S+)\s*$',
                             line.rstrip())
                if m:
                    hits.append((m.group(1).strip(),
                                 int(m.group(2).replace(",", "")),
                                 m.group(3).strip()))
                elif hits and line.strip().startswith("The "):
                    break
    return hits


def main():
    root = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                              else "~/meedlab/tick_project")
    os.chdir(root)
    for need in (REPORT, TABLE):
        if not os.path.isfile(need):
            sys.exit("missing %s" % need)
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)

    wanted = parse_report(REPORT)
    if not wanted:
        sys.exit("could not parse the entomopathogen block in %s" % REPORT)
    print("from %s:" % REPORT)
    for name, reads, lib in wanted:
        print("  %-32s %5d reads  %s" % (name, reads, lib))
    print()

    rows = []
    with open(TABLE, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ic = head.index("feature_id")
        nc = head.index("lowest_named")
        rc = head.index("total_reads")
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) > max(ic, nc, rc):
                rows.append((r[ic], r[nc].strip(), r[rc]))

    resolved, problems = [], []
    for name, reads, lib in wanted:
        match = [r for r in rows if r[1] == name]
        if not match:
            # the report prints a genus where the table may name a species
            match = [r for r in rows if r[1].startswith(name + " ")
                     or name.startswith(r[1] + " ")]
        if len(match) != 1:
            problems.append("%s matched %d features in %s"
                            % (name, len(match), TABLE))
            continue
        fid, named, got = match[0]
        try:
            if int(got) != reads:
                problems.append("%s: report says %d reads, table says %s"
                                % (name, reads, got))
                continue
        except ValueError:
            pass
        resolved.append((fid, named, reads, lib))

    if problems:
        print("STOPPING. The two files disagree:")
        for p in problems:
            print("  %s" % p)
        sys.exit(1)

    seqs, sources = {}, {}
    artifacts = sorted(glob.glob(os.path.join(QIIME, "**", "*.qza"),
                                 recursive=True))
    if not artifacts:
        sys.exit("no .qza found under %s" % QIIME)
    want_ids = {fid for fid, _, _, _ in resolved}
    for art in artifacts:
        try:
            with zipfile.ZipFile(art) as z:
                members = [n for n in z.namelist()
                           if n.endswith("data/dna-sequences.fasta")]
                if not members:
                    continue
                text = z.read(members[0]).decode("utf-8", "replace")
        except (zipfile.BadZipFile, OSError):
            continue
        fid = None
        for line in text.splitlines():
            if line.startswith(">"):
                fid = line[1:].split()[0]
            elif fid in want_ids and fid not in seqs:
                seqs.setdefault(fid, "")
                seqs[fid] += line.strip()
                sources[fid] = art
            elif fid in seqs:
                seqs[fid] += line.strip()

    print("searched %d artifacts under %s" % (len(artifacts), QIIME))
    missing = [f for f in want_ids if f not in seqs]
    if missing:
        print("  no sequence found for: %s" % ", ".join(missing))
    print()

    written = 0
    with open(OUT, "w", encoding="utf-8") as fh:
        for fid, named, reads, lib in resolved:
            if fid not in seqs:
                continue
            fh.write(">%s %s | %d reads | %s | VT-12 ITS1-2\n"
                     % (fid, named, reads, lib))
            s = seqs[fid]
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")
            written += 1
            print("  %-32s %4d bp  from %s" % (named, len(s), sources[fid]))
    print()
    print("wrote %s  (%d of %d sequences)" % (OUT, written, len(resolved)))


if __name__ == "__main__":
    main()
