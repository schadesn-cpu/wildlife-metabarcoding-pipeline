# projects/tick

Analysis code for the tick-borne pathogen surveillance study (Schade et al.).

## Which pipeline version this ran against

The generic pipeline for this study is the state frozen at the tag
`v1-loon-tick` (commit c5e6d84, 15 June 2026, "Frozen loon/tick pipeline
(pre-v2 refactor)"), not the current top level of this repository. The top
level is the post-refactor v2 generic pipeline.

Of the study scripts here, 34 match `v1-loon-tick` byte for byte, 6 differ
from it, and the rest are study-specific and exist only here. So these
scripts, not the repository's top-level copies, are the record of what ran.
`07b_blast_verify.py` in particular differs from the tagged version, and
`07b_blast_verify_patch.py` has no counterpart in the generic pipeline.

## What is not here

Sequence reads go to the NCBI Sequence Read Archive. Feature tables,
taxonomies, the BLAST verification record and the derived genus tables are
deposited with the manuscript.

Sample collection records are withheld. They hold the collection town and
hunter-harvest record for each animal, which together could identify
individual harvest locations and the people who submitted them. Any script
here that needs those values reads them from a file under `metadata/` that
is deliberately absent, and will say so when run. County-level aggregation,
which is what the manuscript reports, is reproducible from the deposited
metadata.

One reconciliation script is also withheld for the same reason.
