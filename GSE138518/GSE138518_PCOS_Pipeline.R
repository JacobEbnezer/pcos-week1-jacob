# Clear env and set working directory
rm(list = ls())
setwd("D:/3rd Year Internship/GSE138518")

# Load DESeq2
library(DESeq2)

# Load raw counts
raw_counts <- read.delim("GSE138518_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)

# Filter out low counts
cleaned_counts <- raw_counts[rowSums(raw_counts) > 10, ]

# Setup sample metadata
sample_metadata <- data.frame(
  condition = factor(c("PCOS", "PCOS", "PCOS", "Control", "Control", "Control")),
  row.names = colnames(cleaned_counts)
)

# Initialize DESeq2 
dds <- DESeqDataSetFromMatrix(countData = cleaned_counts,
                              colData = sample_metadata,
                              design = ~ condition)

# Run pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "PCOS", "Control"))

# Filter for significance (padj < 0.05)
resOrdered <- res[order(res$padj), ]
sig_genes <- subset(resOrdered, padj < 0.05)

# Save to CSV
write.csv(as.data.frame(sig_genes), file = "PCOS_Significant_Genes_Final.csv")

# Generate Volcano Plot
plot(res$log2FoldChange, -log10(res$padj), 
     main = "Volcano Plot: PCOS vs Control", 
     xlab = "Log2 Fold Change", 
     ylab = "-Log10 Adjusted P-value", 
     pch = 20, col = "darkgray")

with(subset(res, padj < 0.05 & abs(log2FoldChange) > 1), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "red"))

abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)
