#!/usr/bin/env python3
"""
07b_blast_verify_patch.py

Patch to update 07b_blast_verify.py BLAST output format to include
species scientific names (sscinames) and accession (sacc) for
downstream Entrez verification.

Apply this change to line ~305 in 07b_blast_verify.py:

OLD:
    "-outfmt", "6 qseqid sseqid stitle pident length evalue bitscore",

NEW:
    "-outfmt", "6 qseqid sacc sscinames staxids pident length evalue bitscore stitle",

Also update 07c_blast_qc_unclassified.py line ~286:

OLD:
    "-outfmt",  "6 qseqid sseqid staxid pident length evalue bitscore",

NEW:
    "-outfmt",  "6 qseqid sacc sscinames staxids pident length evalue bitscore stitle",

This ensures all three scripts (07b, 07c, 07d) produce consistent output
that includes species names directly in the BLAST results, plus accession
numbers that can be verified via 07e_entrez_verify_accessions.py.

Note: sscinames requires a local BLAST database with taxonomy info.
Remote BLAST (-remote) may not support all format specifiers.
For remote BLAST, use the original format and pass results through
07e_entrez_verify_accessions.py for accession verification.
"""

import subprocess
import sys
from pathlib import Path


def patch_file(filepath: Path, old: str, new: str, label: str) -> bool:
    """Replace old string with new in file. Returns True if patched."""
    content = filepath.read_text()
    if old in content:
        content = content.replace(old, new)
        filepath.write_text(content)
        print(f"  ✅ Patched {label} in {filepath.name}")
        return True
    elif new in content:
        print(f"  ⏭️  {filepath.name} already patched for {label}")
        return False
    else:
        print(f"  ⚠️  Could not find target string for {label} in {filepath.name}")
        return False


def main():
    scripts_dir = Path("scripts")
    if not scripts_dir.exists():
        print("Run from project root (where scripts/ directory is)")
        sys.exit(1)

    print("Patching BLAST output formats...\n")

    # Patch 07b
    f07b = scripts_dir / "07b_blast_verify.py"
    if f07b.exists():
        patch_file(
            f07b,
            '"6 qseqid sseqid stitle pident length evalue bitscore"',
            '"6 qseqid sacc sscinames staxids pident length evalue bitscore stitle"',
            "BLAST outfmt"
        )
    else:
        print(f"  ❌ {f07b} not found")

    # Patch 07c
    f07c = scripts_dir / "07c_blast_qc_unclassified.py"
    if f07c.exists():
        patch_file(
            f07c,
            '"6 qseqid sseqid staxid pident length evalue bitscore"',
            '"6 qseqid sacc sscinames staxids pident length evalue bitscore stitle"',
            "BLAST outfmt"
        )
    else:
        print(f"  ❌ {f07c} not found")

    print("\n07d already uses the correct format (sscinames).")
    print("\nDone. The updated BLAST format is:")
    print('  "6 qseqid sacc sscinames staxids pident length evalue bitscore stitle"')
    print("\nTo verify BLAST results, run:")
    print("  python scripts/07e_entrez_verify_accessions.py verify \\")
    print("    --blast-results <blast_output.tsv> \\")
    print("    --marker <marker> --email your@email.edu \\")
    print("    --outdir <output_dir>/")


if __name__ == "__main__":
    main()
