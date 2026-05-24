library(DESeq2)
library(GEOquery)
library(org.Hs.eg.db)
library(AnnotationDbi)

dir.create("data/rds", recursive = TRUE, showWarnings = FALSE)

TARGET_GENES <- c("PCMT1", "PCMTD1", "PCMTD2")


# GSE234297 — Peripheral Blood (Entrez IDs)
TARGET_BLOOD <- c("PCMT1" = "5110", "PCMTD1" = "115294", "PCMTD2" = "55251")

gse_blood  <- getGEO("GSE234297", GSEMatrix = TRUE, getGPL = FALSE)
meta_blood <- pData(gse_blood[[1]])
meta_blood$case_id   <- sub("PeripheralBlood_", "", meta_blood$title)
meta_blood$condition <- factor(meta_blood$`disease state:ch1`,
                               levels = c("Healthy control", "sALS"))
rownames(meta_blood) <- meta_blood$case_id

counts_blood <- read.table(
  gzfile("data/raw/GSE234297_gene_raw_counts.txt.gz"),
  header = TRUE, sep = "\t", row.names = 1
)

common       <- intersect(colnames(counts_blood), rownames(meta_blood))
counts_sub   <- counts_blood[, common]
meta_sub     <- meta_blood[common, , drop = FALSE]

counts_mat <- as.matrix(counts_sub)
storage.mode(counts_mat) <- "integer"

dds_blood <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta_sub,
  design    = ~ condition
)
dds_blood <- dds_blood[rowSums(counts(dds_blood) >= 10) >= 5, ]
dds_blood <- DESeq(dds_blood)
saveRDS(dds_blood, "data/rds/dds_GSE234297.rds")
cat("dds_GSE234297 saved.\n")

# GSE124439 — Cortex (Gene Symbols)
gse_motor  <- getGEO("GSE124439", GSEMatrix = TRUE, getGPL = FALSE)
meta_motor <- pData(gse_motor[[1]])
meta_motor$sample_id <- meta_motor$title

raw_dir <- "data/raw/GSE124439_RAW"
files   <- list.files(raw_dir, pattern = "_counts\\.txt\\.gz$", full.names = TRUE)

count_list <- lapply(files, function(f) {
  df <- read.table(gzfile(f), header = FALSE, sep = "\t", skip = 1)
  colnames(df) <- c("gene", "count")
  setNames(df$count, df$gene)
})

sample_ids   <- gsub("GSM[0-9]+_(.+)_counts\\.txt\\.gz", "\\1", basename(files))
count_matrix <- do.call(cbind, count_list)
colnames(count_matrix) <- sample_ids

common_m <- intersect(colnames(count_matrix), meta_motor$sample_id)
count_m2 <- count_matrix[, common_m, drop = FALSE]
meta_m2  <- meta_motor[match(common_m, meta_motor$sample_id), , drop = FALSE]

meta_m3 <- meta_m2[
  meta_m2$`sample group:ch1` %in% c("ALS Spectrum MND", "Non-Neurological Control"),
  , drop = FALSE
]
meta_m3$condition <- factor(
  meta_m3$`sample group:ch1`,
  levels = c("Non-Neurological Control", "ALS Spectrum MND")
)
rownames(meta_m3) <- meta_m3$sample_id
count_m3 <- count_m2[, rownames(meta_m3), drop = FALSE]

dds_motor <- DESeqDataSetFromMatrix(
  countData = round(count_m3),
  colData   = meta_m3,
  design    = ~ condition
)
dds_motor <- dds_motor[rowSums(counts(dds_motor) >= 10) >= 5, ]
dds_motor <- DESeq(dds_motor)
saveRDS(dds_motor, "data/rds/dds_GSE124439.rds")
cat("dds_GSE124439 saved.\n")

# GSE153960 — ALS Spectrum MND (Ensembl IDs with version)
TARGET_153960 <- c(
  "PCMT1"  = "ENSG00000120265.17",
  "PCMTD1" = "ENSG00000168300.14",
  "PCMTD2" = "ENSG00000203880.12"
)

gse_153960 <- getGEO("GSE153960", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
meta_153960 <- pData(gse_153960)

counts_153960 <- read.table(
  gzfile("data/raw/GSE153960_counts.tsv.gz"),
  header = TRUE, sep = "\t", row.names = 1
)

common_153 <- intersect(colnames(counts_153960), rownames(meta_153960))
counts_sub_153 <- counts_153960[, common_153]
meta_sub_153   <- meta_153960[common_153, , drop = FALSE]

keep_153 <- meta_sub_153$`group:ch1` %in%
  c("ALS Spectrum MND", "Non-Neurological Control")
meta_sub_153  <- meta_sub_153[keep_153, , drop = FALSE]
counts_sub_153 <- counts_sub_153[, rownames(meta_sub_153)]

meta_sub_153$condition <- factor(
  meta_sub_153$`group:ch1`,
  levels = c("Non-Neurological Control", "ALS Spectrum MND")
)

count_int <- as.matrix(counts_sub_153)
storage.mode(count_int) <- "integer"

dds_153960 <- DESeqDataSetFromMatrix(
  countData = count_int,
  colData   = meta_sub_153,
  design    = ~ condition
)
dds_153960 <- dds_153960[rowSums(counts(dds_153960) >= 10) >= 5, ]
dds_153960 <- DESeq(dds_153960)
saveRDS(dds_153960, "data/rds/dds_GSE153960.rds")
cat("dds_GSE153960 saved.\n")