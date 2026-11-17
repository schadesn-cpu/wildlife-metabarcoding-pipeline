#!/usr/bin/env python3
"""I write a MANIFEST.txt at the root of each archive directory.

WHY
---
Two archive directories hold retired output and retired scripts, and neither
says what is in it or why. A directory named 2026-09-24_restructure carries
its reason in its name and nothing else, so anyone who finds a file there has
to guess whether it was superseded, wrong, or simply moved. A manifest turns
that guess into a record.

WHAT I RECORD, AND WHAT I DO NOT
--------------------------------
I record what the filesystem knows: the subdirectory, the date and reason
parsed from its name, and every file with its size and modification time. I
do not invent a reason for an individual file. Where a subdirectory name
carries no reason I say so rather than supplying one, because a wrong reason
in a manifest is worse than an absent one.

I write only MANIFEST.txt, and only inside archive directories that already
exist. I create no directories and I touch no archived file.
"""
import os
import sys
import time

ROOTS = ("results/_archive", "scripts/_archive")
MANIFEST = "MANIFEST.txt"


def human(size):
    """I render a byte count in the largest unit that keeps it readable."""
    for unit in ("B", "K", "M", "G"):
        if size < 1024 or unit == "G":
            return "%6.1f%s" % (size, unit)
        size = size / 1024.0


def describe(name):
    """I split an archive subdirectory name into its date and its reason."""
    parts = name.split("_", 1)
    if len(parts) == 2 and parts[0].count("-") == 2:
        return parts[0], parts[1].replace("_", " ")
    if name.count("-") == 2 and len(name) == 10:
        return name, "no reason recorded in the directory name"
    return "", "no date recorded in the directory name"


def walk(root):
    """I return the subdirectories of an archive root, and its loose files."""
    entries = sorted(os.listdir(root))
    subdirs = [e for e in entries if os.path.isdir(os.path.join(root, e))]
    loose = [e for e in entries
             if os.path.isfile(os.path.join(root, e)) and e != MANIFEST]
    return subdirs, loose


def listing(path, indent):
    """I return one line per file beneath path, with size and mtime."""
    lines, total = [], 0
    for base, _, files in os.walk(path):
        for name in sorted(files):
            if name == MANIFEST:
                continue
            full = os.path.join(base, name)
            try:
                stat = os.stat(full)
            except OSError as problem:
                lines.append("%s%-54s unreadable (%s)"
                             % (indent, name, problem))
                continue
            total += stat.st_size
            rel = os.path.relpath(full, path)
            lines.append("%s%-54s %s  %s"
                         % (indent, rel, human(stat.st_size),
                            time.strftime("%Y-%m-%d",
                                          time.localtime(stat.st_mtime))))
    return lines, total


def write_one(root):
    """I write the manifest for one archive root and return a summary line."""
    subdirs, loose = walk(root)
    out = []
    out.append("MANIFEST for " + root)
    out.append("Written " + time.strftime("%Y-%m-%d") +
               " by scripts/write_archive_manifests.py")
    out.append("")
    out.append("This directory holds retired material. Nothing here is read by")
    out.append("any current script. Sizes and dates are from the filesystem;")
    out.append("reasons are parsed from directory names and are not inferred.")
    out.append("")
    grand, counted = 0, 0
    if loose:
        out.append("=" * 72)
        out.append("Files at the root of " + root)
        out.append("=" * 72)
        for name in loose:
            full = os.path.join(root, name)
            stat = os.stat(full)
            grand += stat.st_size
            counted += 1
            out.append("  %-54s %s  %s"
                       % (name, human(stat.st_size),
                          time.strftime("%Y-%m-%d",
                                        time.localtime(stat.st_mtime))))
        out.append("")
    for name in subdirs:
        date, reason = describe(name)
        path = os.path.join(root, name)
        lines, total = listing(path, "  ")
        grand += total
        counted += len(lines)
        out.append("=" * 72)
        out.append(name)
        out.append("  archived: " + (date if date else "date not in the name"))
        out.append("  reason:   " + reason)
        out.append("  files:    %d, %s total" % (len(lines), human(total).strip()))
        out.append("=" * 72)
        out.extend(lines if lines else ["  (empty)"])
        out.append("")
    out.append("=" * 72)
    out.append("Total: %d files, %s" % (counted, human(grand).strip()))
    out.append("=" * 72)
    target = os.path.join(root, MANIFEST)
    open(target, "w", encoding="utf-8").write("\n".join(out) + "\n")
    return "%s: %d subdirectories, %d files, %s" % (
        target, len(subdirs), counted, human(grand).strip())


def main():
    if not os.path.isdir("results") or not os.path.isdir("scripts"):
        sys.exit("run me from the project root; I cannot see results/ and "
                 "scripts/ here.")
    missing = [r for r in ROOTS if not os.path.isdir(r)]
    if len(missing) == len(ROOTS):
        sys.exit("neither archive directory exists: " + ", ".join(ROOTS))
    for root in ROOTS:
        if root in missing:
            print("skipped, does not exist: " + root)
            continue
        print(write_one(root))


if __name__ == "__main__":
    main()
