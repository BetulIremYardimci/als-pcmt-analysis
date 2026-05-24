#Updated 24 May 2026
#Betül İrem YARDIMCI

library(GEOquery)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

base_url <- "https://zenodo.org/records/6385747/files"

spinal_files <- list(
  list(
    url      = paste0(base_url, "/Thoracic_Spinal_Cord.tar.gz?download=1"),
    destfile = "data/raw/Thoracic_Spinal_Cord.tar.gz"
  ),
  list(
    url      = paste0(base_url, "/Lumbar_Spinal_Cord.tar.gz?download=1"),
    destfile = "data/raw/Lumbar_Spinal_Cord.tar.gz"
  ),
  list(
    url      = paste0(base_url, "/gencode.v30.gene_meta.tsv.gz?download=1"),
    destfile = "data/raw/gencode.v30.gene_meta.tsv.gz"
  )
)

for (f in spinal_files) {
  if (!file.exists(f$destfile)) {
    cat("Downloading:", basename(f$destfile), "\n")
    download.file(f$url, destfile = f$destfile, mode = "wb")
  } else {
    cat("Already exists, skipping:", basename(f$destfile), "\n")
  }
}

untar("data/raw/Thoracic_Spinal_Cord.tar.gz", exdir = "data/raw/")
untar("data/raw/Lumbar_Spinal_Cord.tar.gz",   exdir = "data/raw/")

geo_counts_gse124439 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124439/suppl/GSE124439_RAW.tar"
geo_counts_gse153960 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153960/suppl/GSE153960_counts.tsv.gz"
geo_counts_gse234297 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234297/suppl/GSE234297_gene_raw_counts.txt.gz"

if (!file.exists("data/raw/GSE124439_RAW.tar")) {
  download.file(geo_counts_gse124439, "data/raw/GSE124439_RAW.tar", mode = "wb")
  untar("data/raw/GSE124439_RAW.tar", exdir = "data/raw/GSE124439_RAW/")
}

if (!file.exists("data/raw/GSE153960_counts.tsv.gz"))
  download.file(geo_counts_gse153960, "data/raw/GSE153960_counts.tsv.gz", mode = "wb")

if (!file.exists("data/raw/GSE234297_gene_raw_counts.txt.gz"))
  download.file(geo_counts_gse234297, "data/raw/GSE234297_gene_raw_counts.txt.gz", mode = "wb")
