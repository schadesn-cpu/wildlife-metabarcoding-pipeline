# Quick Start — Wildlife Metabarcoding Pipeline (v2)

A short guide to get you running. The pipeline is **study-agnostic**: it reads
everything from one config file, so you point it at *your* data and markers —
you don't edit any code.

---

## 1. Get access and clone the `v2` branch

You'll get a GitHub invite to the repo — accept it first. Then clone and switch
to the `v2` branch (the current work lives there, **not** on `main`, so don't
skip the checkout):

```bash
git clone https://github.com/schadesn-cpu/wildlife-metabarcoding-pipeline.git
cd wildlife-metabarcoding-pipeline
git checkout v2
```

Tip: note this folder's path (or alias it) — you'll call `pipeline.py` from it.
For example: `alias wmp="python ~/wildlife-metabarcoding-pipeline/pipeline.py"`

---

## 2. Activate the QIIME 2 environment

The pipeline runs inside QIIME 2 2024.5 (amplicon distribution):

```bash
conda activate qiime2-amplicon-2024.5
```

If you don't have QIIME 2 2024.5 installed, install it first per the official
QIIME 2 docs (https://docs.qiime2.org). `environment.yml` in the repo lists the
pipeline's extra Python dependencies.

---

## 3. Set up your project folder

Keep your data and config in **your own folder**, separate from the pipeline
clone. One folder per study:

```
my_study/
  reads/                 # your demultiplexed FASTQs
  metadata/source/       # your master sample sheet
  pipeline_config.yml    # created in the next step
```

Generate a starter config inside that folder:

```bash
cd ~/my_study
python ~/wildlife-metabarcoding-pipeline/pipeline.py init
```

The pipeline auto-discovers `pipeline_config.yml` in your current directory, and
every path resolves **relative to it** — so always run pipeline commands from
your project folder, and your data stays neatly inside it.

---

## 4. Edit the config

Open `pipeline_config.yml` and set, at minimum:

| Setting | What it is |
|---|---|
| `project:` | `name`, `email`, and `root: "."` (your project folder) |
| `active_markers:` | only the markers you ran (e.g. `16S`, `MiFish`) |
| `markers:` | per-marker primers / amplicon length / rarefaction depth |
| `samples.id_regex` | pattern that pulls your sample IDs out of the FASTQ names |
| `metadata.source_sheet` + `source_id_column` | your master sample sheet and its key column |
| `groups.primary.column` | your main variable of interest |

Optional but useful:

- **Confound checks** — add columns to `groups.confounds` to test that your
  primary signal isn't driven by batch/source/season effects.
- **Optional analyses** — `blast`, `presence_absence`, and `rarefaction` are
  **off by default**. Turn them on in the `analyses:` block if you want them.

---

## 5. Check, then run

```bash
# from inside ~/my_study
python ~/wildlife-metabarcoding-pipeline/pipeline.py check                 # validate config + input paths
python ~/wildlife-metabarcoding-pipeline/pipeline.py list                  # show the steps, in order
python ~/wildlife-metabarcoding-pipeline/pipeline.py run --steps all --dry-run  # preview the commands
python ~/wildlife-metabarcoding-pipeline/pipeline.py run --steps all       # run it for real
```

After every step the pipeline prints a **checkpoint**: what it produced, which
QZV to open (and what to look for in it), and the exact next command to run.
Just follow that — you don't need to memorize the order.

---

## Handy to know

- **Run only some steps:** `run --steps taxonomy diversity figures`
- **Re-runs are safe:** finished outputs are skipped automatically; add
  `--force` to redo a step from scratch.
- **Inspect a QIIME result:** `qiime tools view <file>.qzv`, or drag the `.qzv`
  onto https://view.qiime2.org (renders in your browser, nothing is uploaded).
- **Lost?** The checkpoint footer of the last step always tells you what to look
  at and what to run next.

---

*Questions about a specific marker, your metadata layout, or a step that won't
pass `check`? Ask the repo owner — some defaults in the config are examples and
may need adjusting for your study.*
