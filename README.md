#!/usr/bin/env Rscript

# run_motifbreakr_from_rsid.R
# Purpose: Identify TF motif disruptions for rsIDs using motifbreakR (dbSNP155 GRCh38) and save results + PDFs.
# Outputs:
#   - motifbreakR_results.tsv
#   - one PDF per rsID (strong effects)

# --------------------------
# 0) Install/load packages
# --------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "motifbreakR",
  "MotifDb",
  "BSgenome",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BiocParallel"
)

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(motifbreakR)
  library(MotifDb)
  library(BSgenome)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel)
})

# Increase download timeout for Bioconductor resources
options(timeout = 1e6)

# --------------------------
# 1) Inputs
# --------------------------
# NOTE: Your original script only had rs4698413 but later plotted rs77692262/rs79724577.
# Put ALL rsIDs you want analyzed here.
snps <- c("rs4698413", "rs77692262", "rs79724577")

# Alt alleles to use in plotMB (optional; set to NULL to skip allele-specific plotting)
alt_alleles <- list(
  rs4698413  = "T",
  rs77692262 = c("C", "G"),
  rs79724577 = "C"
)

# --------------------------
# 2) Build SNP GRanges for motifbreakR
# --------------------------
snps.mb <- snps.from.rsid(
  rsid = snps,
  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  search.genome = BSgenome.Hsapiens.UCSC.hg38
)

# --------------------------
# 3) Run motifbreakR
# --------------------------
data("motifbreakR_motif")  # built-in motif set distributed with motifbreakR

results <- motifbreakR(
  snpList   = snps.mb,
  filterp   = TRUE,                 # filter by p-value rather than pct score
  pwmList   = motifbreakR_motif,
  threshold = 1e-4,
  method    = "ic",
  bkg       = c(A=0.25, C=0.25, G=0.25, T=0.25),
  BPPARAM   = SerialParam()
)

# --------------------------
# 4) Save results table
# --------------------------
results_df <- as.data.frame(results, row.names = NULL)

# motifPos is a list-column; convert to comma-separated text for TSV export
if ("motifPos" %in% colnames(results_df)) {
  results_df$motifPos <- vapply(results_df$motifPos, function(x) paste(x, collapse = ", "), character(1))
}

out_tsv <- "motifbreakR_results.tsv"
write.table(results_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote results table: ", out_tsv)
message("Rows: ", nrow(results_df), " | Cols: ", ncol(results_df))

# --------------------------
# 5) Plot PDFs per rsID
# --------------------------
for (rs in snps) {
  rs_hits <- results[results$SNP_id == rs]
  if (length(rs_hits) == 0) {
    message("No hits for ", rs, " at current threshold; skipping plot.")
    next
  }

  # If no alt allele specified, create one generic plot
  if (is.null(alt_alleles[[rs]])) {
    pdf(paste0(rs, ".pdf"))
    plotMB(results = results, rsid = rs, effect = "strong")
    dev.off()
    next
  }

  # Otherwise make allele-specific plots
  for (aa in alt_alleles[[rs]]) {
    pdf(paste0(rs, "_alt_", aa, ".pdf"))
    plotMB(results = results, rsid = rs, effect = "strong", altAllele = aa)
    dev.off()
  }
}

# Print SNP object (useful for logging/debug)
snps.mb
