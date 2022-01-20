library("Matrix")

# Generate TCM from DTM --------------------------------------------------------

tissue_name <- "Pancreas"
path_dtm_tissue <- file.path("output", paste0(tissue_name, "_dtm.rds"))
dtm_tissue <- readRDS(path_dtm_tissue)

# Get TCM from DTM by matrix multiplication

# Since the order here is irrelevant, we didn't use sliding windows like word2vec
# or pearson correlation coefficient for gene co-expression inference
# to get the co-occurrence matrix

# Convert document-term matrix to term co-occurrence matrix by inner product
dtm2tcm <- function(x) {
  # Equivalent to t(x) %*% x
  y <- Matrix::crossprod(x)
  # Correct self-cooccurrence, see quanteda:::fcm.dfm
  Matrix::diag(y) <- (Matrix::diag(y) - Matrix::colSums(x)) / 2L
  y
}

tcm_tissue <- dtm2tcm(dtm_tissue)

# Check if there is any NA (possibly due to integer overflow)
range(tcm_tissue)

# Run GloVe on TCM -------------------------------------------------------------

# Vocabulary: gene names
vocab <- colnames(tcm_tissue)

# Quick check if it is symmetrical
identical(tcm_tissue[100, ], tcm_tissue[, 100])

# Make it triangular to save memory
tcm_tissue[lower.tri(tcm_tissue)] <- 0L

# If need to remove diagonal
# mat[!upper.tri(mat)] <- 0

# Check if triangular now
isTriangular(tcm_tissue)

# Convert to dgTMatrix that text2vec needs - this will double the memory usage
tcm_tissue <- as(tcm_tissue, "dgTMatrix")
gc()

# Define model parameters
glove_model <- text2vec::GlobalVectors$new(
  word_vectors_size = 100,
  vocabulary = vocab,
  x_max = 100,
  learning_rate = 0.05
)

# Start training

# On parallelization, reproducibility, and stability:
#
# - text2vec uses all cores by default but the embedding
#   will be different from two runs due to the SGD implementation.
# - To make it fully reproducible, use a single core:
#   <https://github.com/dselivanov/text2vec/issues/251>
# - The embedding stability between runs is also an interesting problem:
#   <http://doi.org/10.18653/v1/N18-1190>

word_vectors_main <- glove_model$fit_transform(x = tcm_tissue, n_iter = 50)

# Save the embedding -----------------------------------------------------------

word_vectors_context <- glove_model$components
word_vectors <- word_vectors_main + t(word_vectors_context)
glove_dump <- glove_model$dump()

glove_dump |>
  saveRDS(file = file.path("output", paste0(tissue_name, "_glove_dump.rds")))
glove_model |>
  saveRDS(file = file.path("output", paste0(tissue_name, "_glove_model.rds")))
word_vectors_main |>
  saveRDS(file = file.path("output", paste0(tissue_name, "_embedding_main.rds")))
word_vectors_context |>
  saveRDS(file = file.path("output", paste0(tissue_name, "_embedding_context.rds")))
word_vectors |>
  saveRDS(file = file.path("output", paste0(tissue_name, "_embedding.rds")))
