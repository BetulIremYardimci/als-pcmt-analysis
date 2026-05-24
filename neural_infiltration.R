#Updated 24 May 2026
#Betül İrem YARDIMCI

library(GSVA)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)

dir.create("results/neural", recursive = TRUE, showWarnings = FALSE)

PCMT_GENES <- c("PCMT1", "PCMTD1", "PCMTD2")

sets_kegg <- c(
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"
)
sets_reactome <- c(
  "REACTOME_NEURONAL_SYSTEM",
  "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
  "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "REACTOME_PROTEIN_PROTEIN_INTERACTIONS_AT_SYNAPSES",
  "REACTOME_NEUREXINS_AND_NEUROLIGINS"
)
sets_go <- c(
  "GOBP_AXON_DEVELOPMENT",
  "GOBP_AXONAL_TRANSPORT",
  "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"
)

neural_genes <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "C2",
          subcategory = "CP:KEGG_LEGACY") %>% filter(gs_name %in% sets_kegg),
  msigdbr(species = "Homo sapiens", category = "C2",
          subcategory = "CP:REACTOME")   %>% filter(gs_name %in% sets_reactome),
  msigdbr(species = "Homo sapiens", category = "C5",
          subcategory = "GO:BP")         %>% filter(gs_name %in% sets_go)
) %>% dplyr::select(gs_name, gene_symbol) %>% distinct()

neural_sets <- lapply(split(neural_genes$gene_symbol, neural_genes$gs_name), unique)
cat("Gene sets loaded:", length(neural_sets), "\n")

run_neural_correlation <- function(vst_rds, dataset_name, id_type = "SYMBOL") {
  
  cat("\n--- Processing:", dataset_name, "---\n")
  vst_raw <- readRDS(vst_rds)
  
  if (id_type == "ENSEMBL_VERSION") {
    rownames(vst_raw) <- gsub("\\..*", "", rownames(vst_raw))
    id_type <- "ENSEMBL"
  }
  
  if (id_type != "SYMBOL") {
    gmap <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = rownames(vst_raw),
                                  columns = "SYMBOL",
                                  keytype = id_type) %>%
      filter(!is.na(SYMBOL)) %>%
      distinct(!!sym(id_type), .keep_all = TRUE)
    colname_key <- if (id_type == "ENSEMBL") "ENSEMBL" else "ENTREZID"
    vst_raw     <- vst_raw[gmap[[colname_key]], ]
    rownames(vst_raw) <- gmap$SYMBOL
    vst_raw <- vst_raw[!duplicated(rownames(vst_raw)), ]
  }
  
  cat("Matrix dimensions:", dim(vst_raw), "\n")
  
  present <- intersect(PCMT_GENES, rownames(vst_raw))
  if (length(present) < length(PCMT_GENES))
    warning("Missing genes: ", paste(setdiff(PCMT_GENES, rownames(vst_raw)), collapse = ", "))
  
  gsva_param   <- gsvaParam(vst_raw, neural_sets, kcdf = "Gaussian")
  neural_scores <- gsva(gsva_param, verbose = FALSE)
  saveRDS(neural_scores, file.path("data/rds", paste0("neural_scores_", dataset_name, ".rds")))
  
  pcmt_expr  <- vst_raw[present, ]
  common     <- intersect(colnames(pcmt_expr), colnames(neural_scores))
  pcmt_sub   <- pcmt_expr[, common]
  neural_sub <- neural_scores[, common]
  
  cor_mat  <- matrix(NA, nrow = length(present), ncol = nrow(neural_sub),
                     dimnames = list(present, rownames(neural_sub)))
  pval_mat <- cor_mat
  
  for (g in present) {
    for (s in rownames(neural_sub)) {
      tst <- cor.test(as.numeric(pcmt_sub[g, ]),
                      as.numeric(neural_sub[s, ]),
                      method = "spearman", exact = FALSE)
      cor_mat[g, s]  <- tst$estimate
      pval_mat[g, s] <- tst$p.value
    }
  }
  
  write.csv(cor_mat,  paste0("results/neural/correlation_matrix_", dataset_name, ".csv"))
  write.csv(pval_mat, paste0("results/neural/pvalue_matrix_", dataset_name, ".csv"))
  cat("Correlation matrices saved for", dataset_name, "\n")
  
  list(cor = cor_mat, pval = pval_mat)
}

res_153960 <- run_neural_correlation("data/rds/vst_mat_GSE153960.rds",
                                     "GSE153960", id_type = "ENSEMBL_VERSION")
res_124439 <- run_neural_correlation("data/rds/vst_mat_GSE124439_cortex.rds",
                                     "GSE124439_cortex", id_type = "SYMBOL")
res_234297 <- run_neural_correlation("data/rds/vst_mat_GSE234297_blood.rds",
                                     "GSE234297_blood", id_type = "ENTREZID")

read_cor <- function(cor_file, pval_file, label) {
  cor_m  <- as.matrix(read.csv(cor_file,  row.names = 1))
  pval_m <- as.matrix(read.csv(pval_file, row.names = 1))
  cor_df  <- as.data.frame(cor_m) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "GeneSet", values_to = "r")
  pval_df <- as.data.frame(pval_m) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "GeneSet", values_to = "pval")
  left_join(cor_df, pval_df, by = c("Gene","GeneSet")) %>%
    mutate(dataset = label)
}

df_153960 <- read_cor("results/neural/correlation_matrix_GSE153960.csv",
                      "results/neural/pvalue_matrix_GSE153960.csv",
                      "GSE153960\n(ALS Spectrum, n=532)")
df_124439 <- read_cor("results/neural/correlation_matrix_GSE124439_cortex.csv",
                      "results/neural/pvalue_matrix_GSE124439_cortex.csv",
                      "GSE124439\n(Cortex, n=162)")
df_234297 <- read_cor("results/neural/correlation_matrix_GSE234297_blood.csv",
                      "results/neural/pvalue_matrix_GSE234297_blood.csv",
                      "GSE234297\n(Blood, n=144)")

plot_df <- bind_rows(df_153960, df_124439, df_234297) %>%
  mutate(
    sig_label = ifelse(pval < 0.001, "***",
                       ifelse(pval < 0.01,  "**",
                              ifelse(pval < 0.05,  "*", ""))),
    GeneSet_short = GeneSet %>%
      gsub("GOBP_|KEGG_|REACTOME_", "", .) %>%
      gsub("_", " ", .) %>%
      str_to_title() %>%
      str_wrap(width = 22),
    Gene    = factor(Gene, levels = c("PCMT1","PCMTD1","PCMTD2")),
    dataset = factor(dataset, levels = c(
      "GSE153960\n(ALS Spectrum, n=532)",
      "GSE124439\n(Cortex, n=162)",
      "GSE234297\n(Blood, n=144)"))
  )

gs_order <- df_153960 %>%
  filter(Gene == "PCMT1") %>%
  mutate(GeneSet_short = GeneSet %>%
           gsub("GOBP_|KEGG_|REACTOME_", "", .) %>%
           gsub("_", " ", .) %>%
           str_to_title() %>%
           str_wrap(width = 22)) %>%
  arrange(r) %>%
  pull(GeneSet_short)

plot_df$GeneSet_short <- factor(plot_df$GeneSet_short, levels = gs_order)

p <- ggplot(plot_df, aes(x = GeneSet_short, y = Gene, fill = r)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sig_label), size = 4.5, fontface = "bold",
            color = "white", vjust = 0.8) +
  geom_text(aes(label = round(r, 2)), size = 3, color = "white", vjust = 2.2) +
  scale_fill_gradientn(
    colors = c("#053061","#2166AC","#4393C3","#92C5DE","#F7F7F7",
               "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"),
    limits = c(-0.65, 0.95),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9),
    name   = "Spearman\nr",
    guide  = guide_colorbar(barwidth = 1, barheight = 10,
                            ticks.colour = "grey30", frame.colour = "grey30")
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_wrap(~ dataset, ncol = 1) +
  labs(
    title    = "PCMT Gene Family × Neural Gene Set Activity in ALS",
    subtitle = "Spearman Correlation | Three Independent Datasets",
    x = NULL, y = NULL,
    caption  = "* p < 0.05   ** p < 0.01   *** p < 0.001"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 10, color = "grey40", hjust = 0),
    plot.caption     = element_text(size = 9,  color = "grey50", hjust = 0),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 9, color = "grey20"),
    axis.text.y      = element_text(size = 11, face = "bold.italic", color = "grey10"),
    axis.ticks       = element_blank(),
    strip.text       = element_text(size = 10, face = "bold", color = "grey20"),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

ggsave("results/neural/PCMT_neural_correlation_3datasets.jpeg",
       p, width = 14, height = 11, dpi = 300, quality = 97)

cat("Combined heatmap saved.\n")