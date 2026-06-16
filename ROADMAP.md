# Pipeline v2 — Roadmap

Remaining work after the v2 refactor, in recommended order. The porting is done
and the run-flow is fixed; nothing below blocks a run, but Tier 1 should be done
before the next *full* end-to-end run on the new dataset. Each item lists why it
matters, rough effort, risk, what it depends on, and how you know it's finished.

Status legend: [ ] not started · [~] partially there · [x] done

---

## Tier 1 — Complete the run-flow (correctness; before the next full run)  ✅ DONE

These are small and high-value: they close the gaps between "runs end-to-end"
and "stops partway." Most were surfaced by the review pass.

### 1. Wire the remaining steps into `pipeline.py run`  ✅
- [x] **Metadata build.** Added a `metadata` stage (after denoise) that builds
  the per-marker QIIME metadata via `make_metadata.py`. Driven by new config keys
  `metadata.source_sheet` / `metadata.source_id_column`; skips with a clear
  warning when they're blank. A full `run` now reaches diversity unaided.
- [x] **Demux QC.** `parse_demux.py` now runs in `step_qc`, after primer
  detection, treated as non-fatal (it warns; it doesn't stop the run).
- [x] **DADA2 advisor.** Left manual by design (trunc-len needs human review).
  `step_denoise` falls back to engine defaults with a logged note pointing at the
  advisor, and the checkpoint hint marks it advisory/manual. Auto-run opt-in not
  added (deferred; not needed).
- **Done when:** `pipeline.py run` goes qc → … → diversity with no manual steps —
  verified by a full `--dry-run` across all eight stages (0 errors).

### 2. Clean two registry / hint mismatches (found in review)  ✅
- [x] **Phantom `rarefaction` entry** removed from the script map.
- [x] **Checkpoint "Next" command** now emits a correct `pipeline.py run --steps
  <stage>` (via a new `stage` field on each registry step).

---

## Tier 2 — Single source of truth (robustness)

Do these once the run-flow is whole. Medium effort, removes whole classes of
"why did it use that value / write there" surprises.

### 3. Reconcile the dual config  ✅ DONE
- [x] **Finding:** the engine's per-marker config was *vestigial* — it loaded
  `config/defaults.yml` + `config/markers/<m>.yml`, stored them, and never read
  them (`get_cfg` had no callers). So it wasn't two configs in conflict; it was
  dead scaffolding that silently ignored anything you set.
- [x] Removed the dead machinery from the engine (`load_combined_config`,
  `load_config_file`, `deep_update`, `Context.config`, `get_cfg`, `--config`).
  `pipeline_config.yml` is now the single source of truth.
- [x] **Bonus bug fixed while tracing it:** cutadapt was getting no primers (a
  silent no-op). `step_denoise` now resolves primers from `markers.<m>.primers`
  or the detected `primers_detected.tsv` and passes them through.
- **Done when:** one config file holds every parameter; the engine reads no
  separate config — verified (engine compiles and parses with `--config` gone,
  no dangling references). You can delete the repo's `config/` tree.

### 4. Route remaining inline paths through `PathBuilder`  ✅ DONE
- [x] Added an `engine_*` block to `PathBuilder` encoding the engine's real
  layout (`qiime2/<marker>/<dataset>/<stage>/<generic-name>` and the matching
  `results/...`), plus `qiime2_root` / `reads_dir` / `primers_detected_tsv` /
  `engine_manifest_tsv`. Routed every inline path f-string in `pipeline.py`
  through them.
- [x] **Finding:** PathBuilder's *existing* methods describe a different (older,
  marker-suffixed, no-dataset) layout than the engine actually writes, so the
  inline literals couldn't just point at them — hence the new `engine_*` methods.
  The old methods are left for the legacy scripts (item 6).
- **Done when:** no marker/dataset path is spelled out inline in `pipeline.py` —
  verified (no `cfg.resolve(f"...")` layout literals remain; each new method
  proven byte-identical to the literal it replaced; full dry-run unchanged).

---

## Tier 3 — Polish (low effort; fold in alongside the above)

### 5. Stale references in docstrings  ✅ DONE
- [x] Updated old numbered-script names in `pipeline.py` docstrings/comments to
  current names.

### 6. (Partially done) Port the remaining legacy analysis/figure scripts
- [x] `taxonomy_table.py` (was `07_taxonomy_table.py`, 854 lines) — ported onto
  the backbone; derives paths from `--marker`, checkpoints, records to the ledger.
- [x] `08_run_diversity_stats.py` + `08b_run_cod_diversity.py` → **unified** into
  `group_diversity.py` (primary + optional confound variables; COD generalized away).
  `diversity_stats.py` and the `cod` step removed. Config-derived, checkpointed,
  fail-loud. Tested (primary + confound dry-run).
- [ ] `run_all_figures.py` (559 lines) — last one; same pattern. Reconcile diversity
  figures with `group_diversity` (which can draw its own) to avoid duplication.
- **Why:** consistency, not function — the remaining three run as-is today via
  `resolve_script`. Each is a mechanical repeat of the proven port pattern
  (header → backbone imports → config-driven defaults/paths via PathBuilder →
  checkpoint + provenance → verbatim domain logic).
- **Done when:** every script the orchestrator calls is on the backbone, at which
  point PathBuilder's two path conventions can collapse into one.

---

## Tier 4 — Optional analysis stages (config-gated; new scope)

Make the standalone analyses that aren't part of the core flow into first-class,
toggleable stages under `analyses:` in the config — same backbone treatment
(config-driven paths, validation, checkpoint + provenance), each skipped loudly
when disabled.

### 7. Rarefaction  ✅ DONE
- [x] `rarefaction.py` (was `06_rarefaction.py`) ported + wired as `step_rarefaction`
  between `taxonomy` and `diversity`, gated by `analyses.rarefaction.enabled`.

### 8. Presence / absence  ✅ DONE
- [x] `11b_presence_absence.py` → `presence_absence.py`, wired as
  `step_presence_absence` after `taxonomy`, gated by `analyses.presence_absence`
  (markers list + thresholds + sample label). Functionally tested end-to-end
  (reads TSV, no QIIME needed).

### 9. BLAST refinement chain  ✅ DONE
- [x] `07d_blast_refine_unresolved` → `blast_refine.py` ported onto the backbone
  (config-driven, derives exported TSV/FASTA/DB/thresholds, checkpoint, fail-loud,
  `--apply` off by default). Candidate-selection + dry-run tested.
- [x] **Refined-taxonomy flow — decided: option (b), built.** `taxonomy_table`
  prefers a refined taxonomy when BLAST is enabled and one exists, loudly, and
  writes `TAXONOMY_SOURCE_<marker>.md` (source, ASVs changed, BLAST settings,
  commit, how to revert) next to the count tables. Nothing auto-rewrites silently;
  classifier runs document themselves too.
- [x] `07c_blast_qc_unclassified` → `blast_qc.py` ported (flags classifier
  conflicts → `confirmed_artefacts` list; conflict taxa config-driven, loon
  defaults genericised). E-utilities network dependency documented as non-fatal.
  `--skip-blast` path tested.
- [x] `07b_blast_verify` → `blast_verify.py` ported (targeted verification;
  config-driven `min_reads` default, advisory only).
- [x] **`step_blast` wired** after `taxonomy`, gated by `analyses.blast.enabled`,
  per-marker over `analyses.blast.markers`, running refine → QC → verify (each
  independent/advisory; one failing doesn't block the others). Disabled by
  default; skips loudly. Dry-run verified (3 tools × markers).
- Entrez accession-verify (`07e`) intentionally dropped (internet + NCBI creds).

---

## Recommended order

~~1 → 2~~ ~~3~~ ~~4~~ ~~5~~ (done) → **6** (taxonomy + diversity done; only `run_all_figures` left) ·
~~7~~ ~~8~~ ~~9~~ (optional-analyses track complete: rarefaction, presence/absence, BLAST).

Rationale: Tier 1 (1–2) is the difference between the pipeline running end-to-end
on the new dataset and stopping partway, and it's small — do it first. Tier 2
(3–4) hardens the foundation: the dual config is the biggest remaining foot-gun,
and centralizing paths is easier once the config story is settled. Tier 3 (5–6)
is cosmetic or optional uniformity — fold 5 into whichever change touches those
files, and treat 6 as a "someday, for consistency" task to tackle only after 3–4
fix the conventions, so the ports aren't done twice.

## Separately: housekeeping (one-time, on your side)
- [ ] `git commit` the current repo state, then `git tag v1-loon-tick` **before**
  dropping the v2 files in, so the loon/tick outputs and old scripts are frozen.
- [ ] Add a `qc:` block to `pipeline_config.yml` if you want to tune thresholds
  (`min_reads`, retention cutoffs) instead of using defaults.
- [ ] One-line note in the README: requires QIIME 2 amplicon 2024.5.
