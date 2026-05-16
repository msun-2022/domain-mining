#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh  —  Viral PLA2-like domain analysis pipeline
#
# Runs all four steps in order.  Pass --dry-run to print commands only.
#
# Usage:
#   ./run_pipeline.sh [--dry-run] [--flank 100] [--threads 4]
#
# Requirements:
#   python >= 3.10  |  cd-hit (bioconda)  |  standard library only for Python
# =============================================================================

set -euo pipefail

# ---------- defaults ---------------------------------------------------------
FLANK=100
THREADS=4
DRY_RUN=false
DEDUP_ID=0.99
CLUSTER_ID=0.70

# ---------- parse arguments --------------------------------------------------
while [[ $# -gt 0 ]]; do
  case $1 in
    --dry-run)   DRY_RUN=true ;;
    --flank)     FLANK="$2";      shift ;;
    --threads)   THREADS="$2";    shift ;;
    --dedup-id)  DEDUP_ID="$2";   shift ;;
    --cluster-id) CLUSTER_ID="$2"; shift ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
  shift
done

run() {
  echo "▶  $*"
  $DRY_RUN && return
  "$@"
}

echo "============================================================"
echo " Viral PLA2-like Domain Pipeline"
echo "  flank=${FLANK}  threads=${THREADS}"
echo "  dedup_id=${DEDUP_ID}  cluster_id=${CLUSTER_ID}"
echo "============================================================"

# Step 1 — Fetch sequences from EBI InterPro
run python scripts/01_fetch_sequences.py \
    --output data/raw/viralexport.json

# Step 2 — Flatten JSON → CSV
run python scripts/02_json_to_csv.py \
    --input  data/raw/viralexport.json \
    --output data/processed/pfam_sequences_flattened.csv

# Step 3 — Extract extended domain sequences
run python scripts/03_extract_domain_seqs.py \
    --boundaries data/processed/pfam_viralpla2_domain_ext.csv \
    --fasta      data/processed/pfam_viralpla2_full.fasta \
    --flank      "${FLANK}" \
    --output     results/pfam_viralpla2_full_domain_ext${FLANK}.fasta

# Step 4 — Cluster sequences
run python scripts/04_cluster_sequences.py \
    --input      results/pfam_viralpla2_full_domain_ext${FLANK}.fasta \
    --outdir     results/clustering \
    --dedup-id   "${DEDUP_ID}" \
    --cluster-id "${CLUSTER_ID}" \
    --threads    "${THREADS}"

echo ""
echo "✓  Pipeline complete."
echo "   Results are in results/"
