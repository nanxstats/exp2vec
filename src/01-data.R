# Get data ---------------------------------------------------------------------

# Get expression reads and metadata links from
# <https://gtexportal.org/home/datasets>
path_reads <- file.path("data", "reads.gz")

if (!file.exists(path_reads)) {
  paste0(
    "https://storage.googleapis.com/",
    "gtex_analysis_v8/rna_seq_data/",
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
  ) |>
    curl::curl_download(destfile = path_reads, quiet = FALSE)
}

path_meta <- file.path("data", "meta.tsv")
if (!file.exists(path_meta)) {
  paste0(
    "https://storage.googleapis.com/",
    "gtex_analysis_v8/annotations/",
    "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
  ) |>
    curl::curl_download(destfile = path_meta, quiet = FALSE)
}

# Load sample name dictionary
meta <- readr::read_tsv(path_meta)

# Create document-term matrix
path_dtm <- file.path("data", "dtm.rds")

if (!file.exists(path_dtm)) {
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)
  dtm <- readr::read_tsv(
    path_reads,
    skip = 2,
    col_types = readr::cols(
      .default = readr::col_integer(),
      Name = readr::col_character(),
      Description = readr::col_character()
    )
  )
  saveRDS(dtm, file = path_dtm)
}

# Load gene read count (DTM)
dtm <- readRDS(path_dtm)

# Preprocess genes -------------------------------------------------------------

# Filter out "useless" genes

# Clean up ensembl id
gene_dtm <- grex::cleanid(dtm[, 1, drop = TRUE])
dtm[, 1] <- gene_dtm

# Map gene id
df_grex <- grex::grex(gene_dtm)
df_entrez <- df_grex[
  !is.na(df_grex$"entrez_id"),
  c("ensembl_id", "entrez_id", "hgnc_symbol", "gene_biotype")
]

# First order by `entrez_id` then de-duplicate `hgnc_symbol`:
# remove the rows with identical `hgnc_symbol` but larger `entrez_id`
df_entrez <- df_entrez[order(as.integer(df_entrez$entrez_id)), ]
idx_hgnc_dup <- which(duplicated(df_entrez$"hgnc_symbol"))
df_entrez[idx_hgnc_dup, ]
df_entrez <- df_entrez[-idx_hgnc_dup, ]

# Keep genes with `ensembl_id` still in the mapping table
dtm <- dtm[which(dtm[, 1, drop = TRUE] %in% df_entrez$"ensembl_id"), ]

# Replace Ensembl IDs with HGNC symbols
dtm[, 1] <- df_entrez[match(dtm[, 1, drop = TRUE], df_entrez$ensembl_id), "hgnc_symbol"]

# Remove gtex annotation
dtm$"Description" <- NULL

# Preprocess samples -----------------------------------------------------------

# List tissue types
unique(meta$SMTS)
sort(table(meta$SMTS), decreasing = TRUE)

# Specify tissue type
tissue_name <- "Pancreas"

get_dtm_tissue <- function(tissue) {
  sample_id_tissue <- as.data.frame(meta[which(meta$SMTS == tissue), "SAMPID"])[, 1]
  dtm_tissue <- dtm[, which(names(dtm) %in% sample_id_tissue)]
  dtm_tissue <- t(as.matrix(dtm_tissue))
  colnames(dtm_tissue) <- dtm$"Name"
  dtm_tissue
}

dtm_tissue <- get_dtm_tissue(tissue_name)
dim(dtm_tissue)

# Save DTM ---------------------------------------------------------------------

path_dtm_tissue <- file.path("output", paste0(tissue_name, "_dtm.rds"))
dtm_tissue |> saveRDS(file = path_dtm_tissue)
