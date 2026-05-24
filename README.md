# als-pcmt-analysis

Expression and pathway analysis of **PCMT1**, **PCMTD1**, and **PCMTD2** genes in Amyotrophic Lateral Sclerosis (ALS).

## Overview

This project investigates the transcriptional profiles of protein-L-isoaspartate O-methyltransferase (PCMT) family genes across ALS patient datasets to explore their potential role in disease-associated molecular mechanisms.

## Datasets

| GEO Accession | Tissue | Platform | Notes |
|---|---|---|---|
| GSE137810 | Spinal cord | RNA-seq | ALS vs control |
| GSE124439 | Spinal cord | RNA-seq | ALS vs control |
| GSE116622 | Blood | RNA-seq | ALS vs control |
| GSE153960 | iPSC-derived motor neurons | RNA-seq | ALS vs control |

## Analysis Scripts

| Script | Description |
|---|---|
| `data_download.R` | GEO dataset retrieval |
| `extract_dds.R` | DESeq2 object construction |
| `wilcoxon_test_expr.R` | Differential expression (Wilcoxon) |
| `expr_visualization.R` | Expression plots |
| `neural_infiltration.R` | Neural infiltration analysis |
| `neural_infiltration_analysis.R` | Extended infiltration analysis |

## Requirements

- R ≥ 4.2
- GEOquery, DESeq2, ggplot2, limma
