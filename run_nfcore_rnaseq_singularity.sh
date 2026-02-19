#!/usr/bin/env bash

# run_nfcore_rnaseq_singularity.sh
# Purpose: Launch nf-core/rnaseq using Singularity with a custom reference genome
# Assumes: nf-core/rnaseq is available and Singularity is configured on the system

set -euo pipefail

# --------------------------
# Singularity cache settings
# --------------------------
export NXF_SINGULARITY_CACHEDIR=/dev/shm/<user>/singularity_cache
export SINGULARITY_TMPDIR=/dev/shm/<user>/singularity_tmp

# --------------------------
# Run nf-core/rnaseq
# --------------------------
nextflow run nf-core/rnaseq \
    -profile singularity \
    --input /path/to/samplesheet.csv \
    --outdir /path/to/results/rnaseq_output \
    --gtf /path/to/reference/Homo_sapiens.GRCh38.gtf.gz \
    --fasta /path/to/reference/Homo_sapiens.GRCh38.fa.gz \
    --igenomes_ignore \
    --genome null \
    -w /path/to/nextflow_workdir
