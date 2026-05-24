library(DESeq2)
library(dplyr)

dir.create("results", recursive = TRUE, showWarnings = FALSE)

TARGET_BLOOD  <- c("PCMT1" = "5110", "PCMTD1" = "115294", "PCMTD2" = "55251")
TARGET_MOTOR  <- c("PCMT1" = "PCMT1", "PCMTD1" = "PCMTD1", "PCMTD2" = "PCMTD2")
TARGET_153960 <- c(
  "PCMT1"  = "ENSG00000120265.17",
  "PCMTD1" = "ENSG00000168300.14",
  "PCMTD2" = "ENSG00000203880.12"
)

dds_blood  <- readRDS("data/rds/dds_GSE234297.rds")
dds_motor  <- readRDS("data/rds/dds_GSE124439.rds")
dds_153960 <- readRDS("data/rds/dds_GSE153960.rds")

run_wilcox <- function(dds, target_genes, dataset_name,
                       control_label, als_label, output_file) {
  
  vst_mat <- assay(vst(dds, blind = FALSE))
  meta    <- as.data.frame(colData(dds))
  
  common  <- intersect(colnames(vst_mat), rownames(meta))
  vst_mat <- vst_mat[, common, drop = FALSE]
  meta    <- meta[common, , drop = FALSE]
  meta    <- meta[meta$condition %in% c(control_label, als_label), , drop = FALSE]
  vst_mat <- vst_mat[, rownames(meta), drop = FALSE]
  meta$condition <- factor(meta$condition, levels = c(control_label, als_label))
  
  cat("Dataset:", dataset_name, "| Groups:\n")
  print(table(meta$condition))
  
  results_list <- lapply(names(target_genes), function(sym) {
    gid <- target_genes[[sym]]
    
    if (!(gid %in% rownames(vst_mat))) {
      return(data.frame(dataset = dataset_name, gene = sym, gene_id = gid,
                        control_n = NA, als_n = NA,
                        control_median = NA, als_median = NA,
                        median_difference_ALS_minus_control = NA,
                        wilcox_pvalue = NA, note = "Gene not found"))
    }
    
    expr_df <- data.frame(
      expression = as.numeric(vst_mat[gid, ]),
      condition  = meta[colnames(vst_mat), "condition"]
    )
    expr_df <- expr_df[!is.na(expr_df$expression) & !is.na(expr_df$condition), ]
    
    ctrl_val <- expr_df$expression[expr_df$condition == control_label]
    als_val  <- expr_df$expression[expr_df$condition == als_label]
    wt       <- wilcox.test(expression ~ condition, data = expr_df,
                            exact = FALSE, alternative = "two.sided")
    
    data.frame(
      dataset = dataset_name, gene = sym, gene_id = gid,
      control_n = length(ctrl_val), als_n = length(als_val),
      control_median = median(ctrl_val, na.rm = TRUE),
      als_median     = median(als_val,  na.rm = TRUE),
      median_difference_ALS_minus_control =
        median(als_val, na.rm = TRUE) - median(ctrl_val, na.rm = TRUE),
      wilcox_pvalue = wt$p.value, note = "OK"
    )
  })
  
  res_df <- do.call(rbind, results_list)
  res_df$wilcox_padj_BH_within_dataset <- p.adjust(res_df$wilcox_pvalue, method = "BH")
  
  write.csv(res_df, output_file, row.names = FALSE)
  cat("Saved:", output_file, "\n")
  return(res_df)
}

wilcox_blood <- run_wilcox(
  dds_blood, TARGET_BLOOD, "GSE234297_Peripheral_Blood",
  "Healthy control", "sALS", "results/wilcoxon_GSE234297_blood.csv"
)

wilcox_motor <- run_wilcox(
  dds_motor, TARGET_MOTOR, "GSE124439_Cortex",
  "Non-Neurological Control", "ALS Spectrum MND",
  "results/wilcoxon_GSE124439_cortex.csv"
)

wilcox_153960 <- run_wilcox(
  dds_153960, TARGET_153960, "GSE153960_ALS_Spectrum_MND",
  "Non-Neurological Control", "ALS Spectrum MND",
  "results/wilcoxon_GSE153960.csv"
)

wilcox_all <- rbind(wilcox_blood, wilcox_motor, wilcox_153960)
wilcox_all$wilcox_padj_BH_all_tests <- p.adjust(wilcox_all$wilcox_pvalue, method = "BH")

write.csv(wilcox_all, "results/wilcoxon_all_three_datasets_PCMT_genes.csv",
          row.names = FALSE)
cat("Combined Wilcoxon results saved.\n")
print(wilcox_all)