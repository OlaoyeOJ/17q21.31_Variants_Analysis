#!/usr/bin/env bash
set -euo pipefail

# run_nfcore_rnaseq_F_vs_R.sh
# Purpose: Run nf-core/rnaseq for F vs R comparison using Singularity
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
SAMPLESHEET="/path/to/samplesheet_F_vs_R.csv"
OUTDIR="/path/to/results/F_vs_R"
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
