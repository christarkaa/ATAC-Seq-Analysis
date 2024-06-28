# Set working directory
setwd("/Users/christophertarkaa/ATACSeq/Peaks")

# Load necessary packages
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(rtracklayer)
library(BiocParallel)
library(EnhancedVolcano)

# Load peak data
treat1.peaks <- read.table("SRR24135553_peaks1.bed", sep="\t")[,1:3]
treat2.peaks <- read.table("SRR24135554_peaks1.bed", sep="\t")[,1:3]
control1.peaks <- read.table("SRR24135555_peaks1.bed", sep="\t")[,1:3]
control2.peaks <- read.table("SRR24135556_peaks1.bed", sep="\t")[,1:3]

# Assign column names
colnames(treat1.peaks) <- c("chrom", "start", "end")
colnames(treat2.peaks) <- c("chrom", "start", "end")
colnames(control1.peaks) <- c("chrom", "start", "end")
colnames(control2.peaks) <- c("chrom", "start", "end")

# Convert dataframes to GRanges objects
treat1.gr <- makeGRangesFromDataFrame(treat1.peaks)
treat2.gr <- makeGRangesFromDataFrame(treat2.peaks)
control1.gr <- makeGRangesFromDataFrame(control1.peaks)
control2.gr <- makeGRangesFromDataFrame(control2.peaks)

# Define the parallel function to fetch raw counts from BigWig files
get_raw_counts_parallel <- function(peaks, bigwig_file, cores = 6) {
  # Import BigWig file
  coverage <- import(bigwig_file, format = "BigWig", which = peaks)
  
  # Use bplapply to process in parallel with the specified number of cores
  counts <- bplapply(seq_along(peaks), function(i) {
    region <- peaks[i]
    overlapping_coverage <- subsetByOverlaps(coverage, region)
    sum(score(overlapping_coverage) * width(overlapping_coverage))
  }, BPPARAM = MulticoreParam(cores))
  
  return(unlist(counts))
}

# Define the merged peaks by union of all peak sets
all.peaks <- reduce(c(treat1.gr, treat2.gr, control1.gr, control2.gr))

# Paths to BigWig files
treat1_bigwig <- "/Users/christophertarkaa/ATACSeq/Peaks/filtered.SRR24135553.bw"
treat2_bigwig <- "/Users/christophertarkaa/ATACSeq/Peaks/filtered.SRR24135554.bw"
control1_bigwig <- "/Users/christophertarkaa/ATACSeq/Peaks/filtered.SRR24135555.bw"
control2_bigwig <- "/Users/christophertarkaa/ATACSeq/Peaks/filtered.SRR24135556.bw"

# Fetch raw counts for treat1_bigwig and measure time
cat("Fetching counts for treat1_bigwig...\n")  # Inform the start of the task
treat1_start_time <- Sys.time()
treat1_counts <- get_raw_counts_parallel(all.peaks, treat1_bigwig, cores = 6)
treat1_end_time <- Sys.time()
treat1_elapsed_time <- treat1_end_time - treat1_start_time
cat("Time taken for treat1_counts:", treat1_elapsed_time, "\n")  # Print elapsed time and move to a new line

# Fetch raw counts for treat2_bigwig and measure time
cat("Fetching counts for treat2_bigwig...\n")  # Inform the start of the task
treat2_start_time <- Sys.time()
treat2_counts <- get_raw_counts_parallel(all.peaks, treat2_bigwig, cores = 6)
treat2_end_time <- Sys.time()
treat2_elapsed_time <- treat2_end_time - treat2_start_time
cat("Time taken for treat2_counts:", treat2_elapsed_time, "\n")  # Print elapsed time and move to a new line

# Fetch raw counts for control1_bigwig and measure time
cat("Fetching counts for control1_bigwig...\n")  # Inform the start of the task
control1_start_time <- Sys.time()
control1_counts <- get_raw_counts_parallel(all.peaks, control1_bigwig, cores = 6)
control1_end_time <- Sys.time()
control1_elapsed_time <- control1_end_time - control1_start_time
cat("Time taken for control1_counts:", control1_elapsed_time, "\n")  # Print elapsed time and move to a new line

# Fetch raw counts for control2_bigwig and measure time
cat("Fetching counts for control2_bigwig...\n")  # Inform the start of the task
control2_start_time <- Sys.time()
control2_counts <- get_raw_counts_parallel(all.peaks, control2_bigwig, cores = 6)
control2_end_time <- Sys.time()
control2_elapsed_time <- control2_end_time - control2_start_time
cat("Time taken for control2_counts:", control2_elapsed_time, "\n")  # Print elapsed time and move to a new line

treat1_counts

# Combine counts into a matrix ensuring all are the same length
counts_matrix <- cbind(treat1_counts, treat2_counts, control1_counts, control2_counts)
colnames(counts_matrix) <- c("treat1", "treat2", "control1", "control2")

# Convert counts to integers
counts_matrix <- round(counts_matrix)

# Create a DESeq2 dataset
condition <- factor(c("treat", "treat", "control", "control"))
colData <- data.frame(row.names = colnames(counts_matrix), condition)
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = colData, design = ~ condition)

# Normalize counts and perform differential analysis
dds <- DESeq(dds)
res <- results(dds)

# Order results by adjusted p-value
res <- res[order(res$padj), ]

# Summarize results
summary(res)

# Perform PCA analysis
vsd <- vst(dds, blind=FALSE)  # or use rlog: vsd <- rlog(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Plot PCA
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  ggtitle("PCA Plot of ATAC-Seq Data") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# Save PCA plot
ggsave("PCA_plot.png", plot = pca_plot)

# Plot MA-plot
png("MA_plot.png")
ma_plot <- plotMA(res, main = "MA-Plot of Differential Peaks", ylim = c(-5, 5))
dev.off()


# Visualize significant peaks (padj < 0.05)
png("Significant_MA_plot.png")
sig_peaks <- res[which(res$padj < 0.05), ]
sig_ma_plot <- plotMA(sig_peaks, main = "Significant Differential Peaks (padj < 0.05)", ylim = c(-5, 5))
dev.off()


# Export significant peaks to a file
sig_peaks_gr <- all.peaks[which(rownames(res) %in% rownames(sig_peaks))]
export(sig_peaks_gr, "significant_peaks.bed")

# Save the counts and results
write.csv(as.data.frame(counts_matrix), "counts_matrix.csv")
write.csv(as.data.frame(res), "differential_peaks_results.csv")

# Generate Volcano Plot
png("Volcano_plot.png")
EnhancedVolcano(res,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano Plot of Differential Peaks',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'grey30',
                lab = NA)  
dev.off()
