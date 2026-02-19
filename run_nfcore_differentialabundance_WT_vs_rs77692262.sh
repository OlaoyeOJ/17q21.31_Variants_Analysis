#!/usr/bin/env bash
set -euo pipefail

# run_nfcore_differentialabundance_WT_vs_rs77692262.sh
# Purpose: Run nf-core/differentialabundance (DESeq2 workflow) for WT vs T1
# Input: STAR+Salmon gene-level counts from nf-core/rnaseq
# Profile: rnaseq + singularity
# Requirements: Nextflow ≥22, Singularity/Apptainer

# --------------------------
# Singularity cache settings
# --------------------------
export NXF_SINGULARITY_CACHEDIR="/path/to/singularity_cache"
export SINGULARITY_TMPDIR="/path/to/singularity_tmp"

# --------------------------
# User-defined paths
# --------------------------
OUTDIR="/path/to/results/WT_vs_rs77692262/nfcore_DESeq"
SAMPLESHEET="/path/to/samplesheet_WT_vs_rs77692262.csv"
CONTRASTS="/path/to/contrasts_WT_vs_rs77692262.csv"
COUNT_MATRIX="/path/to/star_salmon/salmon.merged.gene_counts.tsv"
LENGTH_MATRIX="/path/to/star_salmon/salmon.merged.gene_lengths.tsv"
GTF="/path/to/reference/Homo_sapiens.GRCh38.111.gtf.gz"
WORKDIR="/path/to/nextflow_workdir"

# --------------------------
# Execute pipeline
# --------------------------
nextflow run nf-core/differentialabundance \
  -profile rnaseq,singularity \
  --outdir "${OUTDIR}" \
  --input "${SAMPLESHEET}" \
  --contrasts "${CONTRASTS}" \
  --matrix "${COUNT_MATRIX}" \
  --transcript_length_matrix "${LENGTH_MATRIX}" \
  --gtf "${GTF}" \
  -w "${WORKDIR}"
