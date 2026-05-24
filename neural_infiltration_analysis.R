vst_153960 <- readRDS("data/rds/vst_mat_GSE153960.rds")

cat("Boyut:", dim(vst_153960), "\n")
cat("İlk gen isimleri:\n")
print(rownames(vst_153960)[1:5])
cat("İlk örnek isimleri:\n")
print(colnames(vst_153960)[1:5])

# Metadata da kontrol et
dds_153960 <- readRDS("data/rds/dds_GSE153960.rds")
cat("\nSample grupları:\n")
print(table(dds_153960$condition))

library(GSVA)
library(msigdbr)
library(dplyr)

dir.create("results/ALS/neural", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Gen isimlerini temizle (versiyon numarasını kaldır) ----
vst_clean <- vst_153960
rownames(vst_clean) <- gsub("\\..*", "", rownames(vst_clean))

cat("Temizlendi. İlk 3 gen:\n")
print(rownames(vst_clean)[1:3])

# ---- 2. ENSEMBL → SYMBOL dönüşümü ----
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = rownames(vst_clean),
                                  columns = "SYMBOL",
                                  keytype = "ENSEMBL") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cat("Eşleşen gen sayısı:", nrow(gene_map), "\n")

# Matrisi sembol ile yeniden indexle
vst_symbol <- vst_clean[gene_map$ENSEMBL, ]
rownames(vst_symbol) <- gene_map$SYMBOL
vst_symbol <- vst_symbol[!duplicated(rownames(vst_symbol)), ]

cat("Final matris boyutu:", dim(vst_symbol), "\n")

# PCMT genlerinin matristeki varlığını kontrol et
cat("PCMT genleri:\n")
print(vst_symbol[c("PCMT1", "PCMTD1", "PCMTD2"), 1:3])



# ---- 3. Neural Gen Setleri ----
# BRCA'da kullandığımız setlerin aynısı
sets_kegg <- c("KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
               "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY")

sets_reactome <- c("REACTOME_NEURONAL_SYSTEM",
                   "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
                   "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
                   "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
                   "REACTOME_PROTEIN_PROTEIN_INTERACTIONS_AT_SYNAPSES",
                   "REACTOME_NEUREXINS_AND_NEUROLIGINS")

sets_go <- c("GOBP_AXON_DEVELOPMENT",
             "GOBP_AXONAL_TRANSPORT",
             "GOBP_ANTEROGRADE_AXONAL_TRANSPORT")

neural_genes_combined <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "C2", 
          subcategory = "CP:KEGG_LEGACY") %>% filter(gs_name %in% sets_kegg),
  msigdbr(species = "Homo sapiens", category = "C2", 
          subcategory = "CP:REACTOME") %>% filter(gs_name %in% sets_reactome),
  msigdbr(species = "Homo sapiens", category = "C5", 
          subcategory = "GO:BP") %>% filter(gs_name %in% sets_go)
) %>% dplyr::select(gs_name, gene_symbol) %>% distinct()

neural_geneset_list <- split(neural_genes_combined$gene_symbol,
                             neural_genes_combined$gs_name)
neural_geneset_list <- lapply(neural_geneset_list, unique)

cat("Gen seti sayısı:", length(neural_geneset_list), "\n")

# ---- 4. ssGSEA ----
cat("ssGSEA hesaplanıyor...\n")
gsva_param <- gsvaParam(vst_symbol, neural_geneset_list, kcdf = "Gaussian")
neural_scores_als <- gsva(gsva_param, verbose = TRUE)

cat("ssGSEA tamamlandı. Boyut:", dim(neural_scores_als), "\n")
saveRDS(neural_scores_als, "data/rds/neural_scores_GSE153960.rds")

# ---- 5. Spearman Korelasyon Matrisi ----
cat("Korelasyon hesaplanıyor...\n")

pcmt_genes <- c("PCMT1", "PCMTD1", "PCMTD2")
pcmt_expr  <- vst_symbol[pcmt_genes, ]

# Ortak örnekler
common <- intersect(colnames(pcmt_expr), colnames(neural_scores_als))
cat("Ortak örnek:", length(common), "\n")

pcmt_sub   <- pcmt_expr[, common]
neural_sub <- neural_scores_als[, common]

cor_mat_als  <- matrix(NA, nrow = 3, ncol = nrow(neural_sub),
                       dimnames = list(pcmt_genes, rownames(neural_sub)))
pval_mat_als <- cor_mat_als

for (g in pcmt_genes) {
  for (s in rownames(neural_sub)) {
    test <- cor.test(as.numeric(pcmt_sub[g, ]),
                     as.numeric(neural_sub[s, ]),
                     method = "spearman", exact = FALSE)
    cor_mat_als[g, s]  <- test$estimate
    pval_mat_als[g, s] <- test$p.value
  }
}

print(round(cor_mat_als, 3))

write.csv(cor_mat_als,  "results/ALS/neural/correlation_matrix_GSE153960.csv")
write.csv(pval_mat_als, "results/ALS/neural/pvalue_matrix_GSE153960.csv")
cat("Kaydedildi.\n")



library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

dir.create("results/ALS/neural", recursive = TRUE, showWarnings = FALSE)

cor_mat_als  <- as.matrix(read.csv("results/ALS/neural/correlation_matrix_GSE153960.csv", 
                                   row.names = 1))
pval_mat_als <- as.matrix(read.csv("results/ALS/neural/pvalue_matrix_GSE153960.csv", 
                                   row.names = 1))

# Uzun formata çevir
cor_df <- as.data.frame(cor_mat_als) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "GeneSet", values_to = "r")

pval_df <- as.data.frame(pval_mat_als) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "GeneSet", values_to = "pval")

plot_df <- left_join(cor_df, pval_df, by = c("Gene", "GeneSet")) %>%
  mutate(
    sig_label = ifelse(pval < 0.001, "***",
                       ifelse(pval < 0.01,  "**",
                              ifelse(pval < 0.05,  "*", ""))),
    GeneSet_short = GeneSet %>%
      gsub("GOBP_|KEGG_|REACTOME_", "", .) %>%
      gsub("_", " ", .) %>%
      str_to_title() %>%
      str_wrap(width = 22),
    Gene = factor(Gene, levels = c("PCMT1", "PCMTD1", "PCMTD2"))
  )

# GeneSet sıralama
geneset_order <- plot_df %>%
  group_by(GeneSet_short) %>%
  summarise(mean_r = mean(r)) %>%
  arrange(mean_r) %>%
  pull(GeneSet_short)

plot_df$GeneSet_short <- factor(plot_df$GeneSet_short, levels = geneset_order)

p <- ggplot(plot_df, aes(x = GeneSet_short, y = Gene, fill = r)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_label),
            size = 5, fontface = "bold", color = "white", vjust = 0.8) +
  geom_text(aes(label = round(r, 2)),
            size = 3.2, color = "white", vjust = 2.2) +
  scale_fill_gradientn(
    colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
               "#F7F7F7",
               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    limits = c(-0.5, 0.9),
    breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8),
    name   = "Spearman\nCorrelation (r)",
    guide  = guide_colorbar(barwidth = 1.2, barheight = 12,
                            ticks.colour = "grey30", frame.colour = "grey30")
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title    = "PCMT Gene Family vs. Neural Gene Set Activity",
    subtitle = "Spearman Correlation | GSE153960 ALS Spectrum MND | n = 532",
    x = NULL, y = NULL,
    caption  = "* p < 0.05   ** p < 0.01   *** p < 0.001"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(size = 15, face = "bold", hjust = 0,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0,
                                 margin = margin(b = 12)),
    plot.caption  = element_text(size = 9, color = "grey50", hjust = 0),
    plot.margin   = margin(20, 20, 15, 20),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                               size = 10, color = "grey20"),
    axis.text.y = element_text(size = 12, face = "bold.italic", color = "grey10"),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    legend.position = "right"
  )

ggsave("results/ALS/neural/PCMT_neural_correlation_heatmap_GSE153960.jpeg",
       p, width = 14, height = 4.5, dpi = 300, quality = 97)

cat("Heatmap kaydedildi.\n")
