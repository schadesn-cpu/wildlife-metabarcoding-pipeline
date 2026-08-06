# Improvements & Sticky Spots — running log

Working notes on things to fix or reconsider before (or shortly after) the public
push. Ordered roughly by impact. Checked items are done; unchecked are open.

---

## Directory / path system  (HIGH — decide before finishing the naming sweep)

### AGREED TARGET LAYOUT
```
project_root/
  reads/                                    project-level input
  classifiers/                              project-level references
  logs/run_manifest.jsonl                   provenance ledger (dates live HERE, not in folders)

  qiime2/{marker}/{dataset_id}/             INTERMEDIATES ONLY (.qza/.qzv)
      imported/  dada2/  taxonomy/  tree/  diversity/

  results/{marker}/{dataset_id}/[{question}/]{analysis}/   DELIVERABLES ONLY
      taxonomy/                  count table + barplot (figure lives here)
      diversity/rarefied_{depth}/  stats + PCoA/alpha (depth scoped HERE, not at top)
      differential_abundance/    ANCOM-BC table + plot
      presence_absence/          detection table + plot

  submission_figures/            FLAT, curated, config-driven, never auto-overwritten
      figure_1_diversity_pcoa.svg (+ .png)
```

Design rules:
- `dataset_id` replaces `all` — a meaningful sample-set name set in config.
- `{marker}` stays its own level (multi-marker grouping; processing is per-marker).
- `{question}` layer is OPTIONAL — only for a marker serving two purposes (e.g. 18S
  diet vs parasites). Single-purpose markers skip it. [DECISION: optional, not always-present]
- Figures live WITH their analysis (no sibling `figures/`); a figure is a property
  of an analysis, not a peer of it.
- Rarefaction depth scopes inside `diversity/` only (taxonomy/ANCOM-BC aren't rarefied).
- `qiime2/` = intermediates only; `results/` = deliverables only. Never cross.

### Figures — two tiers
- WORKING figures: auto-generated every run, live inside their analysis folder, PNG
  for quick looking. Regenerating them never touches submission figures.
- SUBMISSION figures: a curated set assembled from a `submission_figures:` config
  block (each entry = source analysis + structured overrides: size, DPI, palette,
  group order, panel labels, titles on/off, format). Emitted as PNG + editable SVG
  into a FLAT `submission_figures/` with figure-numbered names.
- Polish tiers: pipeline auto-figures → config assembles submission set (structured
  tweaks in code) → OPTIONAL hand-tuning of the SVG in Inkscape/Illustrator for the
  last 5%. Inkscape is NOT on the critical path; code-driven polish gets ~90-95%.

### Implementation stages (do in order, verify each)
1. [ ] Consolidate to one PathBuilder scheme (remove dead old methods; migrate the
   2 real callers). Fixes the `results/` leak. ← DOING FIRST
2. [ ] `dataset_id` from config replaces `all` default; thread through.
3. [ ] Scope depth to `diversity/rarefied_{depth}/`.
4. [ ] Optional `{question}` layer in results paths.
5. [ ] Figures move inside analysis folders; add `submission_figures/` + config block
   + PNG/SVG emission.

- [ ] **Two overlapping PathBuilder schemes.** `config_loader.py` has an older set
  of methods that build `qiime2/{marker}/{stage}/` (no dataset level, ~lines
  516–589) and the `engine_*` set that builds `qiime2/{marker}/{dataset}/{stage}/`
  (612–690). Depending on which a script calls, the *same* artifact can land in
  two different folders. Consolidate to ONE family (the `engine_*` dataset-scoped
  one looks canonical) and delete the other.
- [ ] **QIIME artifacts in `results/` — actually a `.qzv` convention question.**
  Hunt found the qiime-in-results files are mostly `.qzv` VISUALIZATIONS (taxa
  barplots, ANCOM, PERMANOVA/group-significance from group_diversity + 05_run_full),
  not stray `.qza` data. DECISION NEEDED:
    (a) `.qza`→qiime2/, `.qzv`→results/ (viz = deliverable; current behavior), or
    (b) all `.qza`+`.qzv`→qiime2/; results/ holds only QIIME-free outputs
        (exported .tsv tables + .png/.svg figures). Cleaner for non-QIIME users.
  Leaning (b) per "results/ should be openable without QIIME." Still confirm any
  actual stray `.qza` (vs `.qzv`) by looking at a real run tree.
- [ ] **`all` is a non-descriptive dataset slug.** Every path is
  `{marker}/all/...`. Replace the default with a meaningful dataset id set in the
  config (e.g. `dataset_id: loonblood2023`), so paths read
  `qiime2/16S/loonblood2023/...` and `results/16S/loonblood2023/figures/...`.
- [ ] **Inconsistent slug conventions.** Core-metrics dirs use
  `{depth}_{analysis}` (e.g. `r5000_DvT`) while everything else uses `{dataset}`.
  Pick one convention and document it.

### Open design decisions (see chat)
- Keep `{marker}/` as its own path level (recommended — preserves per-marker
  grouping for multi-marker studies) vs. fold marker into one token
  (`16S_loon_2023`).
- Put the run **date** in folder names vs. keep it only in the provenance ledger
  (`logs/run_manifest.jsonl`). Leaning: ledger only, to avoid dated-folder clutter
  and orphaned partial re-runs — unless dated release snapshots are wanted.

---

## Script rules pass  (in progress)

- [x] Rule 5 — no `pass`-in-`except` (config_loader type-probe; primer_advisor
  now counts-and-warns on unreadable reads).
- [x] Rule 7 — no nested functions. All extracted to module level
  (`blast_verify`, `blast_qc`, `presence_absence`, `taxonomy_table` ×3).
  Note: the first nested-def audit missed control-flow-nested defs (e.g. a def
  inside an `else`); corrected to walk all ancestors and re-scanned.
- [ ] Rule 1 — accurate names. ~59 single-char names left across ~17 files.
  **Hold `config_loader.py` until the directory design is settled** (that file
  gets rewritten by the path work). Other scripts can proceed.
- [ ] Rule 5 — review the 12 broad `except Exception` (most are defensible
  top-level `main()` wrappers; narrow the mid-code ones).

## Voice
- [x] Decision: neutral/impersonal; "author" ok for attribution; no first-person
  "I", no "you" in code-logic comments, no personal names. Second person is
  allowed in `--help`/usage text.
