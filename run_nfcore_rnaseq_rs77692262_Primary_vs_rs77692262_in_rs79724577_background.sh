#!/usr/bin/env bash
set -euo pipefail

# run_nfcore_rnaseq_rs77692262_Primary_vs_rs77692262_in_rs79724577_background.sh
# Purpose: Run nf-core/rnaseq for rs77692262_Primary_vs_rs77692262_in_rs79724577_background comparison using Singularity
# Reference: GRCh38 (custom GTF + FASTA)
# Requirements: Nextflow ≥22, Singularity/Apptainer

# --------------------------
# Singularity cache settings
# --------------------------
export NXF_SINGULARITY_CACHEDIR="/path/to/singularity_cache"
export SINGULARITY_TMPDIR="/path/to/singularity_tmp"

# --------------------------
# User-defined paths
# --------------------------
SAMPLESHEET="/path/to/samplesheet_rs77692262_Primary_vs_rs77692262_in_rs79724577_background.csv"
OUTDIR="/path/to/results/rs77692262_Primary_vs_rs77692262_in_rs79724577_background"
GTF="/path/to/reference/Homo_sapiens.GRCh38.111.gtf.gz"
FASTA="/path/to/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
WORKDIR="/path/to/nextflow_workdir"

# --------------------------
# Execute pipeline
# --------------------------
nextflow run nf-core/rnaseq \
  -profile singularity \
  --input "${SAMPLESHEET}" \
  --outdir "${OUTDIR}" \
  --gtf "${GTF}" \
  --fasta "${FASTA}" \
  --igenomes_ignore \
  --genome null \
  -w "${WORKDIR}"
