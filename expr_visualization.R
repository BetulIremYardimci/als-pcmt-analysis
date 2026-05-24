#Updated 24 May 2026
#Betül İrem YARDIMCI

library(DESeq2)
library(ggplot2)
library(dplyr)

dir.create("plots",   recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

dds_blood  <- readRDS("data/rds/dds_GSE234297.rds")
dds_motor  <- readRDS("data/rds/dds_GSE124439.rds")
dds_153960 <- readRDS("data/rds/dds_GSE153960.rds")

wilcox_all <- read.csv("results/wilcoxon_all_three_datasets_PCMT_genes.csv")

TARGET_BLOOD  <- data.frame(gene = c("PCMT1","PCMTD1","PCMTD2"),
                            gene_id = c("5110","115294","55251"),
                            stringsAsFactors = FALSE)
TARGET_MOTOR  <- data.frame(gene = c("PCMT1","PCMTD1","PCMTD2"),
                            gene_id = c("PCMT1","PCMTD1","PCMTD2"),
                            stringsAsFactors = FALSE)
TARGET_153960 <- data.frame(gene = c("PCMT1","PCMTD1","PCMTD2"),
                            gene_id = c("ENSG00000120265.17",
                                        "ENSG00000168300.14",
                                        "ENSG00000203880.12"),
                            stringsAsFactors = FALSE)

get_deseq_results <- function(dds, target_df, contrast) {
  res    <- results(dds, contrast = contrast, alpha = 0.05,
                    pAdjustMethod = "BH", independentFiltering = FALSE)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  merge(target_df, res_df, by = "gene_id", all.x = TRUE)
}

res_blood  <- get_deseq_results(dds_blood,  TARGET_BLOOD,
                                c("condition","sALS","Healthy control"))
res_motor  <- get_deseq_results(dds_motor,  TARGET_MOTOR,
                                c("condition","ALS Spectrum MND","Non-Neurological Control"))
res_153960 <- get_deseq_results(dds_153960, TARGET_153960,
                                c("condition","ALS Spectrum MND","Non-Neurological Control"))

plot_boxplot <- function(dds, target_df, deseq_res, wilcox_df,
                         dataset_name, dataset_key, dataset_label,
                         control_label, als_label,
                         control_display, als_display,
                         output_dir) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  vst_mat <- assay(vst(dds, blind = FALSE))
  saveRDS(vst_mat, file.path("data/rds", paste0("vst_mat_", dataset_name, ".rds")))
  
  meta    <- as.data.frame(colData(dds))
  common  <- intersect(colnames(vst_mat), rownames(meta))
  vst_mat <- vst_mat[, common, drop = FALSE]
  meta    <- meta[common, , drop = FALSE]
  meta$condition <- factor(meta$condition,
                           levels = c(control_label, als_label),
                           labels = c(control_display, als_display))
  
  expr_long <- do.call(rbind, lapply(seq_len(nrow(target_df)), function(i) {
    sym <- target_df$gene[i]
    gid <- target_df$gene_id[i]
    if (!(gid %in% rownames(vst_mat))) return(NULL)
    data.frame(
      vst_expr  = as.numeric(vst_mat[gid, ]),
      condition = meta[colnames(vst_mat), "condition"],
      gene      = sym, stringsAsFactors = FALSE
    )
  }))
  expr_long <- expr_long[!is.na(expr_long$condition), ]
  expr_long$condition <- factor(expr_long$condition,
                                levels = c(control_display, als_display))
  expr_long$gene <- factor(expr_long$gene, levels = c("PCMT1","PCMTD1","PCMTD2"))
  
  wlx_sub  <- wilcox_df[wilcox_df$dataset == dataset_key, ]
  label_df <- merge(
    deseq_res[, c("gene","log2FoldChange","padj")],
    wlx_sub[, c("gene","wilcox_pvalue","wilcox_padj_BH_within_dataset")],
    by = "gene", all.x = TRUE
  )
  label_df$label <- paste0(
    "log2FC = ", round(label_df$log2FoldChange, 2), "\n",
    "DESeq2 padj = ", ifelse(is.na(label_df$padj), "NA",
                             formatC(label_df$padj, format = "e", digits = 2)), "\n",
    "Wilcox p = ",   ifelse(is.na(label_df$wilcox_pvalue), "NA",
                            formatC(label_df$wilcox_pvalue, format = "e", digits = 2)), "\n",
    "Wilcox padj = ", ifelse(is.na(label_df$wilcox_padj_BH_within_dataset), "NA",
                             formatC(label_df$wilcox_padj_BH_within_dataset,
                                     format = "e", digits = 2))
  )
  ymax     <- aggregate(vst_expr ~ gene, data = expr_long, max)
  label_df <- merge(label_df, ymax, by = "gene", all.x = TRUE)
  label_df$y    <- label_df$vst_expr + 0.35
  label_df$gene <- factor(label_df$gene, levels = c("PCMT1","PCMTD1","PCMTD2"))
  
  p <- ggplot(expr_long, aes(x = condition, y = vst_expr, fill = condition)) +
    geom_boxplot(outlier.size = 0.6, width = 0.5, alpha = 0.85) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    geom_text(data = label_df,
              aes(x = 1.5, y = y, label = label),
              inherit.aes = FALSE, size = 3.0, color = "gray30") +
    scale_fill_manual(values = setNames(c("#2980B9","#C0392B"),
                                        c(control_display, als_display))) +
    labs(title    = "PCMT1 / PCMTD1 / PCMTD2 Expression in ALS",
         subtitle = dataset_label,
         x = NULL, y = "VST Expression", fill = "Condition") +
    theme_bw(base_size = 13) +
    theme(plot.title      = element_text(face = "bold", size = 14),
          plot.subtitle   = element_text(size = 12),
          strip.text      = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          axis.text.x     = element_text(size = 10))
  
  ggsave(file.path(output_dir, paste0(dataset_name, "_boxplot.pdf")),
         p, width = 12, height = 5)
  ggsave(file.path(output_dir, paste0(dataset_name, "_boxplot.png")),
         p, width = 12, height = 5, dpi = 300)
  cat("Saved:", output_dir, "\n")
  return(p)
}

plot_boxplot(dds_blood, TARGET_BLOOD, res_blood, wilcox_all,
             "GSE234297_blood", "GSE234297_Peripheral_Blood",
             "GSE234297 | Peripheral Blood | VST | sALS n=96, Control n=48",
             "Healthy control", "sALS", "Healthy control", "sALS",
             "plots/GSE234297")

plot_boxplot(dds_motor, TARGET_MOTOR, res_motor, wilcox_all,
             "GSE124439_cortex", "GSE124439_Cortex",
             "GSE124439 | Cortex | VST | ALS n=145, Control n=17",
             "Non-Neurological Control", "ALS Spectrum MND",
             "Control", "ALS", "plots/GSE124439")

plot_boxplot(dds_153960, TARGET_153960, res_153960, wilcox_all,
             "GSE153960", "GSE153960_ALS_Spectrum_MND",
             "GSE153960 | VST | ALS n=442, Control n=90",
             "Non-Neurological Control", "ALS Spectrum MND",
             "Non-Neurological Control", "ALS Spectrum MND",
             "plots/GSE153960")

cat("All boxplots complete.\n")