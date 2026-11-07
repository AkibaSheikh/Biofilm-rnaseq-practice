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
