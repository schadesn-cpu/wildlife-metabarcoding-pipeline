#!/usr/bin/env python3
"""check_which_verified_table.py  --  READ ONLY. I write nothing.

THE QUESTION
------------
Three BLAST-verified 16S taxonomy tables exist, all 89,970 rows, all tracing
to the full-length classifier:

    taxonomy_16S_blast_verified.tsv               2026-07-08
    taxonomy_16S_blast_verified_PRE_DROPLIST.tsv  2026-07-08
    taxonomy_16S_droplist_applied.tsv             2026-08-13   <- newest

They differ in what was dropped, not in how ASVs were classified. So there is
no taxonomy risk between them -- but the paper has to cite one, and if two
different figures read two different generations then two numbers in the
manuscript were computed under different drop rules. That is worth five
minutes to rule out.

A grep for the filenames returns 14 scripts and tells you nothing about which
way the data flowed. This separates the writers from the readers, and resolves
the case where the path is assigned to a variable first.

HOW I CLASSIFY EACH MENTION
---------------------------
  WRITE    the same line carries a write idiom (to_csv, write.table, open w)
  READ     the same line carries a read idiom (read_csv, read_table, fread)
  VAR->W   the line assigns the path to a variable, and that variable is
           later written
  VAR->R   same, but later read
  COMMENT  the mention is inside a comment -- not counted either way
  ?        none of the above; printed so you can look at it yourself

A script that writes a table is its producer. A script that reads one and
then writes into results/tables or results/figures is feeding the paper, and
those are the ones whose choice of generation matters.

Longest filename wins the match, so blast_verified.tsv is never credited with
a mention of blast_verified_PRE_DROPLIST.tsv.
"""
import os
import re
import sys

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                          else "~/meedlab/tick_project")
SCRIPTDIR = os.path.join(ROOT, "scripts")

# longest first: PRE_DROPLIST must win over the bare blast_verified stem
TOKENS = [
    "taxonomy_16S_blast_verified_PRE_DROPLIST",
    "taxonomy_16S_droplist_applied",
    "taxonomy_16S_blast_verified",
]

WRITE_IDIOM = re.compile(
    r'\bto_csv\b|\bto_excel\b|\bwrite\.table\b|\bwrite\.csv\b|\bwrite_tsv\b'
    r'|\bwrite_csv\b|\bwriteLines\b|\bsavetxt\b|\bfwrite\b|\bto_pickle\b'
    r'|\bopen\s*\([^)]*["\']w["\']|\.write\s*\(|\bsaveRDS\b|>\s*["\']?\$?\w*/')
READ_IDIOM = re.compile(
    r'\bread_csv\b|\bread_table\b|\bread_tsv\b|\bread_excel\b|\bread\.table\b'
    r'|\bread\.csv\b|\bread_delim\b|\bfread\b|\breadRDS\b|\breadLines\b'
    r'|\bpd\.read\w*|\bnp\.loadtxt\b|\bgenfromtxt\b|\bread_biom\b'
    r'|\bopen\s*\(\s*[^,)]*\)|\bopen\s*\([^)]*["\']r["\']')
OUTPUT_IDIOM = re.compile(r'results/(tables|figures|analysis)/')

ASSIGN_PY = re.compile(r'^\s*([A-Za-z_]\w*)\s*(?::[^=]+)?=\s*[^=]')
ASSIGN_R = re.compile(r'^\s*([A-Za-z_.][\w.]*)\s*(?:<-|=)\s*[^=]')


def strip_comment(line, ext):
    """I return the code part of the line, and whether the token was in code."""
    if ext in (".py", ".sh", ".r", ".R"):
        out, q = [], None
        for ch in line:
            if q:
                out.append(ch)
                if ch == q:
                    q = None
                continue
            if ch in "\"'":
                q = ch
                out.append(ch)
                continue
            if ch == "#":
                break
            out.append(ch)
        return "".join(out)
    return line


def classify(lines, idx, ext):
    line = lines[idx]
    code = strip_comment(line, ext)
    hit_in_code = any(t in code for t in TOKENS)
    if not hit_in_code:
        return "COMMENT", None
    if WRITE_IDIOM.search(code):
        return "WRITE", None
    if READ_IDIOM.search(code):
        return "READ", None
    m = (ASSIGN_R if ext in (".r", ".R") else ASSIGN_PY).match(code)
    if m:
        var = m.group(1)
        usage = re.compile(r'\b%s\b' % re.escape(var))
        for j, other in enumerate(lines):
            if j == idx or not usage.search(strip_comment(other, ext)):
                continue
            oc = strip_comment(other, ext)
            if WRITE_IDIOM.search(oc):
                return "VAR->W", (var, j + 1)
            if READ_IDIOM.search(oc):
                return "VAR->R", (var, j + 1)
        return "?", (var, None)
    return "?", None


def token_of(line):
    for t in TOKENS:          # longest first
        if t in line:
            return t
    return None


def main():
    if not os.path.isdir(SCRIPTDIR):
        sys.exit("no scripts/ at %s -- run me from the project root" % ROOT)

    findings = {t: {"writers": [], "readers": [], "unclear": []} for t in TOKENS}
    feeds_output = {}

    for dirpath, dirnames, filenames in os.walk(SCRIPTDIR):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in sorted(filenames):
            ext = os.path.splitext(fn)[1]
            if ext not in (".py", ".R", ".r", ".sh"):
                continue
            path = os.path.join(dirpath, fn)
            rel = os.path.relpath(path, ROOT)
            try:
                with open(path, encoding="utf-8", errors="replace") as fh:
                    lines = fh.read().splitlines()
            except OSError as exc:
                print("  !! cannot read %s: %s" % (rel, exc))
                continue

            writes_output = any(
                OUTPUT_IDIOM.search(strip_comment(ln, ext))
                and WRITE_IDIOM.search(strip_comment(ln, ext))
                for ln in lines)
            feeds_output[rel] = writes_output

            for i, ln in enumerate(lines):
                tok = token_of(ln)
                if not tok:
                    continue
                verdict, detail = classify(lines, i, ext)
                rec = (rel, i + 1, verdict, detail, ln.strip()[:88])
                if verdict in ("WRITE", "VAR->W"):
                    findings[tok]["writers"].append(rec)
                elif verdict in ("READ", "VAR->R"):
                    findings[tok]["readers"].append(rec)
                elif verdict == "?":
                    findings[tok]["unclear"].append(rec)

    bar = "=" * 78
    for tok in TOKENS:
        f = findings[tok]
        print(bar)
        print(tok + ".tsv")
        print(bar)
        for label, recs in (("PRODUCED BY", f["writers"]),
                            ("READ BY", f["readers"]),
                            ("UNRESOLVED", f["unclear"])):
            if not recs:
                print("  %s: none" % label)
                continue
            print("  %s:" % label)
            for rel, lineno, verdict, detail, src in recs:
                star = "  *paper*" if feeds_output.get(rel) else ""
                where = ""
                if detail and detail[1]:
                    where = "  (via %s, line %d)" % detail
                elif detail:
                    where = "  (var %s, never read or written)" % detail[0]
                print("    %-46s :%-4d %-7s%s%s"
                      % (rel, lineno, verdict, star, where))
                print("        %s" % src)
        print()

    print(bar)
    print("VERDICT")
    print(bar)
    # A script that PRODUCES one of the three necessarily reads the generation
    # before it. That is the pipeline, not a disagreement, so producers are
    # excluded from the conflict test -- otherwise the droplist step always
    # looks like a second opinion on the table it was built to supersede.
    producers = {r[0] for tok in TOKENS for r in findings[tok]["writers"]}
    paper_readers = {}
    for tok in TOKENS:
        who = sorted({r[0] for r in findings[tok]["readers"]
                      if feeds_output.get(r[0]) and r[0] not in producers})
        if who:
            paper_readers[tok] = who
    if producers:
        print("  excluded as pipeline producers (they read the generation they")
        print("  supersede, by design):")
        for p in sorted(producers):
            print("    %s" % p)
        print()
    if not paper_readers:
        print("  No script that writes into results/tables, results/figures or")
        print("  results/analysis reads any of the three by name. The genus")
        print("  tables are reaching the taxonomy some other way -- most likely")
        print("  a variable built from a config, or the .qza directly.")
        print("  Grep for taxonomy_verified/ next.")
    elif len(paper_readers) == 1:
        tok, who = next(iter(paper_readers.items()))
        print("  One generation feeds the paper: %s.tsv" % tok)
        print("  Cite that file. Read by:")
        for w in who:
            print("    %s" % w)
    else:
        print("  *** CONFLICT: paper outputs are built from %d different"
              % len(paper_readers))
        print("      generations of the verified taxonomy. ***")
        for tok, who in paper_readers.items():
            print("    %s.tsv" % tok)
            for w in who:
                print("        %s" % w)
        print()
        print("  These differ in what was dropped, so two numbers in the")
        print("  manuscript were computed under different drop rules. Settle")
        print("  which generation is canonical before the numbers are final.")
    print()
    print("  * = script also writes into results/tables, figures or analysis")
    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
