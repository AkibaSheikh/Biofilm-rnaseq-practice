# ---- Step 2: Load the data ----
# Dataset: GSE302777 (P. aeruginosa, biofilm vs planktonic, +/- midazolam)

getwd()
list.files("data")

library(DESeq2)
library(tidyverse)

# Read a file whose separator I don't know in advance
read_any <- function(f) {
  first <- readLines(f, n = 1)
  seps <- c(",", ";", "\t", " ", "|")
  n <- sapply(seps, function(s) length(regmatches(first, gregexpr(s, first, fixed = TRUE))[[1]]))
  s <- seps[which.max(n)]
  cat("Using separator: [", s, "]\n")
  read.table(f, sep = if (s == " ") "" else s, header = TRUE,
             check.names = FALSE, quote = "\"", comment.char = "")
}

counts_raw <- read_any("data/GSE302777_counts_GEO_benzodiazepines_PAO1DK.csv.gz")

# The real column names are in the first row (the file has generic Column1... labels)
stopifnot(counts_raw[1, 1] == "Geneid")
colnames(counts_raw) <- as.character(unlist(counts_raw[1, ]))
counts_raw <- counts_raw[-1, ]

# Gene IDs become row names; keep only the 12 sample columns as numbers
rownames(counts_raw) <- counts_raw$Geneid
counts <- counts_raw[, colnames(counts_raw) != "Geneid"]
counts[] <- lapply(counts, as.numeric)

# Sample table, built from the sample names
# DKC = untreated control, DKM = midazolam; BF = biofilm, PL = planktonic
coldata <- data.frame(
  growth    = ifelse(grepl("BF$", colnames(counts)), "biofilm", "planktonic"),
  treatment = ifelse(grepl("^DKC", colnames(counts)), "control", "midazolam"),
  row.names = colnames(counts)
)
coldata$growth    <- factor(coldata$growth,    levels = c("planktonic", "biofilm"))
coldata$treatment <- factor(coldata$treatment, levels = c("control", "midazolam"))

# ---- Checks ----
dim(counts)                            # expect 5678 x 12
head(counts)
coldata
table(coldata$growth, coldata$treatment)   # expect 3 in every cell
all(counts == round(counts))           # expect TRUE
stopifnot(all(colnames(counts) == rownames(coldata)))

ls()

# ---- Step 3: Quality control and filtering ----

# Total reads counted per sample
lib_sizes <- colSums(counts)
lib_sizes
summary(lib_sizes)

qc_df <- data.frame(sample = names(lib_sizes),
                    total = as.numeric(lib_sizes),
                    growth = coldata$growth,
                    treatment = coldata$treatment)

p_qc <- ggplot(qc_df, aes(x = sample, y = total / 1e6, fill = growth)) +
  geom_col() +
  labs(x = "Sample", y = "Total counts (millions)",
       title = "Total counts per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_qc
ggsave("figures/total_counts_per_sample.png", p_qc, width = 7, height = 4, dpi = 300)


genes_before <- nrow(counts)

keep <- rowSums(counts) >= 10
counts_filtered <- counts[keep, ]

genes_after <- nrow(counts_filtered)
c(before = genes_before, after = genes_after, removed = genes_before - genes_after)

qc_summary <- data.frame(genes_before, genes_after, removed = genes_before - genes_after)
write.csv(qc_summary, "results/filtering_summary.csv", row.names = FALSE)
write.csv(qc_df, "results/total_counts_per_sample.csv", row.names = FALSE)

# ---- Step 4: Differential expression with DESeq2 ----
# Main comparison: biofilm vs planktonic, untreated controls only

ctrl <- coldata$treatment == "control"
counts_ctrl  <- round(as.matrix(counts_filtered[, ctrl]))
coldata_ctrl <- droplevels(coldata[ctrl, ])

stopifnot(all(colnames(counts_ctrl) == rownames(coldata_ctrl)))
table(coldata_ctrl$growth)      # expect 3 planktonic, 3 biofilm

dds <- DESeqDataSetFromMatrix(countData = counts_ctrl,
                              colData   = coldata_ctrl,
                              design    = ~ growth)
dds <- DESeq(dds)
resultsNames(dds)

# Unshrunk results: use these for p-values and counting significant genes
res_unshrunk <- results(dds, contrast = c("growth", "biofilm", "planktonic"))

# Shrunk fold changes: use these for plots (volcano, ranking)
res <- lfcShrink(dds, coef = "growth_biofilm_vs_planktonic", type = "apeglm")

resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), "results/DE_results.csv")

summary(res_unshrunk, alpha = 0.05)

up   <- sum(res_unshrunk$padj < 0.05 & res_unshrunk$log2FoldChange >  1, na.rm = TRUE)
down <- sum(res_unshrunk$padj < 0.05 & res_unshrunk$log2FoldChange < -1, na.rm = TRUE)
c(up_in_biofilm = up, down_in_biofilm = down)

sum(res_unshrunk$padj < 0.05 & abs(res_unshrunk$log2FoldChange) > 2, na.rm = TRUE)

# ---- Step 5: PCA ----

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "growth", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = growth, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "PCA of control samples (biofilm vs planktonic)") +
  theme_bw()

p_pca
ggsave("figures/pca.png", p_pca, width = 6, height = 5, dpi = 300)

library(ggrepel)
geom_text_repel(size = 3, show.legend = FALSE)

# ---- Step 6: Volcano plot ----
library(EnhancedVolcano)

top10 <- head(rownames(resOrdered), 10)

p_volcano <- EnhancedVolcano(res,
                             lab = rownames(res),
                             x = "log2FoldChange",
                             y = "padj",
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             selectLab = top10,
                             title = "Biofilm vs planktonic (untreated controls)",
                             subtitle = "Positive log2FC = higher in biofilm",
                             drawConnectors = TRUE,
                             max.overlaps = 20)

p_volcano
ggsave("figures/volcano.png", p_volcano, width = 7, height = 6, dpi = 300)

# ---- Step 7: Heatmap of top 30 genes ----
library(pheatmap)

resOrdered_clean <- resOrdered[!is.na(resOrdered$padj), ]
top <- head(rownames(resOrdered_clean), 30)

mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat)

pheatmap(mat,
         annotation_col = coldata_ctrl["growth"],
         main = "Top 30 genes: biofilm vs planktonic")

png("figures/heatmap_top30.png", width = 2000, height = 1800, res = 300)
pheatmap(mat,
         annotation_col = coldata_ctrl["growth"],
         main = "Top 30 genes: biofilm vs planktonic")
dev.off()
