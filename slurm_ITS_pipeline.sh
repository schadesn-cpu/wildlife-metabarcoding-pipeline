#!/bin/bash
# =============================================================================
# slurm_ITS_pipeline.sh
# SLURM pipeline for ITS1-2 taxonomy classification and table generation
# Loon project — MEED Lab, UNH
#
# Usage:
#   bash slurm_ITS_pipeline.sh          # submits all jobs with dependencies
#   bash slurm_ITS_pipeline.sh --dry-run # prints job commands without submitting
#
# Jobs submitted (in order, with dependencies):
#   1. classify   — UNITE classifier against ITS rep-seqs (~30-60 min)
#   2. tabulate   — taxonomy visualization QZV (~2 min)
#   3. taxtable   — 08_taxonomy_table.py → TSV count table (~5 min)
#   4. diversity  — core-metrics-non-phylogenetic (~10 min)
#
# Adapting for other markers:
#   Change MARKER, CLASSIFIER, TABLE, REPSEQS, METADATA, and OUTDIR
#   variables below. The job dependency chain is marker-agnostic.
#   For 16S add --p-n-jobs to match --cpus-per-task.
#   For MiFish/cytb: taxonomy already done — start from taxtable step only.
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT=/home/users/sne24/meedlab/loon_project
MARKER=ITS1-2
CLASSIFIER=${PROJECT}/classifiers/unite-ver10-99-nb-classifier.qza
REPSEQS=${PROJECT}/qiime2/dada2/rep-seqs_ITS1-2.qza
TABLE=${PROJECT}/qiime2/dada2/table_ITS1-2.qza
METADATA=${PROJECT}/metadata/qiime/metadata_16S.tsv   # same sample set
TAXONOMY_OUT=${PROJECT}/qiime2/taxonomy/taxonomy_ITS1-2.qza
RESULTS_OUT=${PROJECT}/results/ITS/all
LOG_DIR=${PROJECT}/logs/ITS

# SLURM settings — tune to your cluster
PARTITION=shared
CPUS=16           # match --p-n-jobs in classify step
MEM=64G
TIME=04:00:00     # 4 hours — classify-sklearn is the bottleneck
EMAIL=sne24@unh.edu

# Rarefaction depth — check qiime2/dada2/table_ITS1-2.qzv first
# Use a depth that retains most samples. 1000 is a safe starting point for ITS.
RAREFY_DEPTH=1000

# ── Dry run flag ───────────────────────────────────────────────────────────────
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN — no jobs will be submitted ==="
fi

# ── Setup ──────────────────────────────────────────────────────────────────────
mkdir -p ${LOG_DIR} ${RESULTS_OUT}/taxonomy ${RESULTS_OUT}/diversity

submit() {
    # Usage: submit <jobname> <dependency_jid_or_empty> <script_content>
    local name=$1
    local dep=$2
    local script=$3
    local tmpscript=$(mktemp /tmp/slurm_${name}_XXXX.sh)
    echo "${script}" > ${tmpscript}

    if [[ "${DRY_RUN}" == "true" ]]; then
        echo ""
        echo "=== Job: ${name} (dep: ${dep:-none}) ==="
        cat ${tmpscript}
        rm ${tmpscript}
        echo "FAKE_JID_${name}"
        return
    fi

    local dep_flag=""
    if [[ -n "${dep}" ]]; then
        dep_flag="--dependency=afterok:${dep}"
    fi

    local jid
    jid=$(sbatch ${dep_flag} ${tmpscript} | awk '{print $NF}')
    rm ${tmpscript}
    echo "${jid}"
}

# ── Job 1: Taxonomy classification ────────────────────────────────────────────
CLASSIFY_SCRIPT="#!/bin/bash
#SBATCH --job-name=${MARKER}_classify
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --output=${LOG_DIR}/slurm_classify_%j.out
#SBATCH --error=${LOG_DIR}/slurm_classify_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}

source ~/.bashrc
conda activate qiime2-amplicon-2024.5
cd ${PROJECT}

echo \"=== ${MARKER} taxonomy classification ===\"
echo \"Start: \$(date)\"

qiime feature-classifier classify-sklearn \\
    --i-classifier  ${CLASSIFIER} \\
    --i-reads       ${REPSEQS} \\
    --o-classification ${TAXONOMY_OUT} \\
    --p-n-jobs      ${CPUS} \\
    --p-confidence  0.7 \\
    --verbose

echo \"Done: \$(date)\"
"

echo "Submitting Job 1: classify..."
JID1=$(submit "${MARKER}_classify" "" "${CLASSIFY_SCRIPT}")
echo "  Job 1 ID: ${JID1}"

# ── Job 2: Taxonomy tabulate (QZV visualization) ─────────────────────────────
TABULATE_SCRIPT="#!/bin/bash
#SBATCH --job-name=${MARKER}_tabulate
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=${LOG_DIR}/slurm_tabulate_%j.out
#SBATCH --error=${LOG_DIR}/slurm_tabulate_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${EMAIL}

source ~/.bashrc
conda activate qiime2-amplicon-2024.5
cd ${PROJECT}

echo \"=== ${MARKER} taxonomy tabulate ===\"
echo \"Start: \$(date)\"

qiime metadata tabulate \\
    --m-input-file ${TAXONOMY_OUT} \\
    --o-visualization ${RESULTS_OUT}/taxonomy/taxonomy_${MARKER}.qzv

echo \"Done: \$(date)\"
"

echo "Submitting Job 2: tabulate (after job ${JID1})..."
JID2=$(submit "${MARKER}_tabulate" "${JID1}" "${TABULATE_SCRIPT}")
echo "  Job 2 ID: ${JID2}"

# ── Job 3: Taxonomy table (TSV count table via 08_taxonomy_table.py) ──────────
TAXTABLE_SCRIPT="#!/bin/bash
#SBATCH --job-name=${MARKER}_taxtable
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=${LOG_DIR}/slurm_taxtable_%j.out
#SBATCH --error=${LOG_DIR}/slurm_taxtable_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}

source ~/.bashrc
conda activate qiime2-amplicon-2024.5
cd ${PROJECT}

echo \"=== ${MARKER} taxonomy table ===\"
echo \"Start: \$(date)\"

python scripts/08_taxonomy_table.py \\
    --taxonomy  ${TAXONOMY_OUT} \\
    --table     ${TABLE} \\
    --marker    ITS \\
    --include   Fungi \\
    --exclude   Bacteria,Archaea,Viruses,Plantae,Animalia \\
    --outdir    ${RESULTS_OUT}/taxonomy/

echo \"Done: \$(date)\"
"

echo "Submitting Job 3: taxtable (after job ${JID1})..."
JID3=$(submit "${MARKER}_taxtable" "${JID1}" "${TAXTABLE_SCRIPT}")
echo "  Job 3 ID: ${JID3}"

# ── Job 4: Core diversity metrics ─────────────────────────────────────────────
DIVERSITY_SCRIPT="#!/bin/bash
#SBATCH --job-name=${MARKER}_diversity
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=${LOG_DIR}/slurm_diversity_%j.out
#SBATCH --error=${LOG_DIR}/slurm_diversity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}

source ~/.bashrc
conda activate qiime2-amplicon-2024.5
cd ${PROJECT}

echo \"=== ${MARKER} core diversity ===\"
echo \"Start: \$(date)\"

mkdir -p ${RESULTS_OUT}/diversity

qiime diversity core-metrics \\
    --i-table           ${TABLE} \\
    --p-sampling-depth  ${RAREFY_DEPTH} \\
    --m-metadata-file   ${METADATA} \\
    --output-dir        ${RESULTS_OUT}/diversity/core_metrics_depth${RAREFY_DEPTH}/ \\
    --p-n-jobs-or-threads 4

echo \"Done: \$(date)\"
"

echo "Submitting Job 4: diversity (after job ${JID1})..."
JID4=$(submit "${MARKER}_diversity" "${JID1}" "${DIVERSITY_SCRIPT}")
echo "  Job 4 ID: ${JID4}"

# ── Summary ────────────────────────────────────────────────────────────────────
echo ""
echo "=== Pipeline submitted ==="
echo "  Job 1 (classify)  : ${JID1}"
echo "  Job 2 (tabulate)  : ${JID2}  [runs after ${JID1}]"
echo "  Job 3 (taxtable)  : ${JID3}  [runs after ${JID1}]"
echo "  Job 4 (diversity) : ${JID4}  [runs after ${JID1}]"
echo ""
echo "Monitor with:"
echo "  squeue -u sne24"
echo "  tail -f ${LOG_DIR}/slurm_classify_*.out"
echo ""
echo "NOTE: Check rarefaction depth before diversity runs."
echo "  Inspect: qiime2/dada2/table_ITS1-2.qzv"
echo "  Current depth: ${RAREFY_DEPTH} — update RAREFY_DEPTH in this script if needed."
