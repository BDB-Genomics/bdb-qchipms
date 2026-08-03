#!/usr/bin/env Rscript

# Statistical Enrichment & Volcano Plot for qChIP-MS

args <- commandArgs(trailingOnly = TRUE)
input_tsv <- ifelse(length(args) >= 1, args[1], "data/protein_intensities.tsv")
output_tsv <- ifelse(length(args) >= 2, args[2], "results/enrichment/qchip_ms_enriched_proteins.tsv")
output_png <- ifelse(length(args) >= 3, args[3], "results/plots/qchip_ms_volcano.png")
fdr_cutoff <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)
log2fc_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), 1.0)

df <- read.delim(input_tsv, stringsAsFactors = FALSE, check.names = FALSE)

# Locate replicate columns
igg_cols <- grep("IgG", colnames(df), value = TRUE)
trf2_cols <- grep("TRF2|TERF2", colnames(df), value = TRUE)

cat("Processing", nrow(df), "proteins across", length(trf2_cols), "treatment and", length(igg_cols), "control replicates.\n")

# Calculate means and stats if not pre-computed
trf2_mat <- as.matrix(df[, trf2_cols])
igg_mat <- as.matrix(df[, igg_cols])

# Ensure numeric values
storage.mode(trf2_mat) <- "numeric"
storage.mode(igg_mat) <- "numeric"

trf2_mean <- rowMeans(trf2_mat, na.rm = TRUE)
igg_mean <- rowMeans(igg_mat, na.rm = TRUE)

log2FC <- trf2_mean - igg_mean

pvals <- apply(cbind(trf2_mat, igg_mat), 1, function(row) {
  t_vals <- row[1:length(trf2_cols)]
  i_vals <- row[(length(trf2_cols)+1):length(row)]
  t_vals <- t_vals[!is.na(t_vals)]
  i_vals <- i_vals[!is.na(i_vals)]
  if (length(t_vals) >= 2 && length(i_vals) >= 2 && sd(t_vals) + sd(i_vals) > 0) {
    return(t.test(t_vals, i_vals, var.equal = FALSE)$p.value)
  } else {
    return(1.0)
  }
})

neg_log10_pval <- -log10(pvals)
fdr <- p.adjust(pvals, method = "BH")

res <- data.frame(
  Protein_ID = df[["Protein IDs"]],
  Gene_Name = df[["Gene names"]],
  Log2FC = round(log2FC, 3),
  PValue = pvals,
  NegLog10PVal = round(neg_log10_pval, 3),
  FDR = fdr,
  IsEnriched = (log2FC >= log2fc_cutoff & fdr <= fdr_cutoff),
  stringsAsFactors = FALSE
)

# Sort by significance
res <- res[order(res$FDR, -res$Log2FC), ]

write.table(res, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Enrichment statistics saved to:", output_tsv, "\n")

# Generate Volcano Plot
png(output_png, width = 2000, height = 1800, res = 300)
par(mar = c(5, 5, 4, 2))

cols <- ifelse(res$IsEnriched, "#E41A1C", "#999999")
plot(res$Log2FC, res$NegLog10PVal, col = cols, pch = 19, cex = 0.8,
     xlab = "Log2 Fold Change (TRF2 / IgG)",
     ylab = "-Log10 (p-value)",
     main = "qChIP-MS TERF2 Enrichment Volcano Plot")

abline(v = log2fc_cutoff, lty = 2, col = "blue")
abline(h = -log10(0.05), lty = 2, col = "blue")

# Label top hits
top_hits <- head(res[res$IsEnriched, ], 10)
if (nrow(top_hits) > 0) {
  text(top_hits$Log2FC, top_hits$NegLog10PVal, labels = top_hits$Gene_Name, pos = 4, cex = 0.7, col = "#800000")
}

dev.off()
cat("Volcano plot saved to:", output_png, "\n")
