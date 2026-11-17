#!/usr/bin/env python3
"""check_18s_primers.py

I ask the reads what primers made them, rather than asking a spreadsheet.
I write nothing.

WHY THIS, AND WHY NOW
---------------------
The 18S taxonomy was produced by pr2-18S-classifier-v4v8.qza, and that
classifier's own provenance says it was extracted with 515F
(GTGYCAGCMGCCGCGGTAA) and 926R (CCGYCAATTYMTTTRAGTTT) -- the 16S primers.
Those are universal SSU primers and in eukaryotes they span 18S V4-V5, so the
classifier is not obviously wrong. It is only wrong if the 18S libraries were
made with something else.

PrimerTable.xlsx says the 18S amplicon is V9, 1391F/1510R. The compiled
workbook's Demux sheet says "18S V4" in all twelve batches. Those disagree,
and two V9 classifiers were built and neither was used. Nothing on disk so far
settles which primers are actually in the 18S data.

The reads do. Whatever trimmed the primers off the 18S reads recorded which
primers it trimmed, in the provenance of the demultiplexed or trimmed
artifacts. And the sequences themselves carry the answer: an 18S V4-V5
fragment and an 18S V9 fragment are different lengths and different sequence.

So I look at three things, in order of how hard they are to argue with:

  1. Every cutadapt or trim action in the 18S pipeline's provenance, and the
     primers it names.
  2. The length distribution of the 18S representative sequences. The V9
     fragment from 1391F/1510R runs about 210 to 310 bp; the eukaryotic
     515F/926R fragment runs about 550 to 650 bp before truncation. Section
     2.3 truncates 18S at 220 single-end, so the TRIMMED length will not
     separate them -- it is the pre-truncation record that does.
  3. The most abundant representative sequences themselves, printed whole, so
     one BLAST settles it if steps 1 and 2 leave it open.

WHAT I AM NOT DOING
-------------------
Deciding. The 16S comparison already showed 41,015 of 89,970 ASVs differ at
genus between the two 16S taxonomies, so the cost of guessing here is not
small. If these three do not agree, the next step is a BLAST of one sequence
against nt, not a judgement call.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_18s_primers.py
"""
import os
import re
import sys
import zipfile

PRIMERS = {
    "GTGYCAGCMGCCGCGGTAA": "515F  (16S V4 fwd / euk SSU)",
    "CCGYCAATTYMTTTRAGTTT": "926R  (16S V5 rev / euk SSU)",
    "GGACTACNVGGGTWTCTAAT": "806R  (16S V4 rev)",
    "GTACACACCGCCCGTC": "1391F (18S V9 fwd)",
    "TGATCCTTCTGCAGGTTCACCTAC": "1510R (18S V9 rev)",
    "CCAGCASCYGCGGTAATTCC": "TAReuk454FWD1 (18S V4 fwd)",
    "ACTTTCGTTCTTGATYRA": "TAReukREV3 (18S V4 rev)",
}
KEYS = ("f_primer", "r_primer", "front", "adapter", "anywhere", "p_trim_left",
        "trunc_len", "trim_left")
PARAM = r'^[^\S\n]*(?:-[^\S\n]*)?%s:[^\S\n]*(\S.*?)[^\S\n]*$'


def actions(path):
    """I yield (member, action name, text) for every action.yaml inside."""
    try:
        with zipfile.ZipFile(path) as z:
            for n in z.namelist():
                if n.endswith("action/action.yaml"):
                    with z.open(n) as fh:
                        text = fh.read().decode("utf-8", "replace")
                    m = re.search(r'^[^\S\n]+action:[^\S\n]*(\S+)', text,
                                  flags=re.M)
                    yield n, (m.group(1) if m else "?"), text
    except Exception as exc:                              # noqa: BLE001
        print("    unreadable: %s" % exc)


def name_primer(value):
    v = value.strip().strip("'\"").upper()
    for seq, label in PRIMERS.items():
        if seq in v or v in seq:
            return "%s  <- %s" % (value, label)
    return value


def sequences(path, limit=3):
    """I return the first few representative sequences, with lengths."""
    try:
        with zipfile.ZipFile(path) as z:
            hit = [n for n in z.namelist()
                   if n.endswith(("data/dna-sequences.fasta",
                                  "data/sequences.fasta"))]
            if not hit:
                return [], []
            with z.open(sorted(hit, key=len)[0]) as fh:
                text = fh.read().decode("utf-8", "replace")
    except Exception:                                     # noqa: BLE001
        return [], []
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs[:limit], [len(s) for s in seqs]


def main():
    if not os.path.isdir("qiime2"):
        sys.exit("Run me from the project root; I see no qiime2/ directory.")

    for marker in ("18S", "16S"):
        root = os.path.join("qiime2", marker)
        if not os.path.isdir(root):
            continue
        arts = []
        for here, dirs, names in os.walk(root):
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for n in sorted(names):
                if n.endswith(".qza") and "classifier" not in n.lower():
                    arts.append(os.path.join(here, n))

        print("=" * 74)
        print("%s  --  %d artifacts" % (marker, len(arts)))
        print("=" * 74)

        seen = set()
        for path in arts:
            for member, act, text in actions(path):
                if not re.search(r'cutadapt|trim|demux|denoise|dada2', act, re.I):
                    continue
                got = {}
                for key in KEYS:
                    for m in re.finditer(PARAM % key, text, flags=re.M):
                        v = m.group(1)
                        if v not in ("null", "0", "[]"):
                            got.setdefault(key, []).append(v)
                if not got:
                    continue
                sig = (act, repr(sorted(got.items())))
                if sig in seen:
                    continue
                seen.add(sig)
                print()
                print("  action %s" % act)
                print("    seen in %s" % path)
                for key, vals in got.items():
                    for v in vals[:4]:
                        print("      %-12s %s" % (key, name_primer(v)))

        reps = [p for p in arts
                if re.search(r'rep.?seq|representative', p, re.I)]
        for path in reps[:2]:
            seqs, lens = sequences(path)
            if not lens:
                continue
            print()
            print("  representative sequences: %s" % path)
            print("    %d sequences, length min %d median %d max %d"
                  % (len(lens), min(lens), sorted(lens)[len(lens) // 2],
                     max(lens)))
            if marker == "18S":
                print("    V9 from 1391F/1510R is ~210-310 bp;")
                print("    eukaryotic 515F/926R is ~550-650 bp untruncated.")
                for s in seqs[:2]:
                    print()
                    print("    >%d bp" % len(s))
                    for i in range(0, min(len(s), 400), 80):
                        print("    %s" % s[i:i + 80])

    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
