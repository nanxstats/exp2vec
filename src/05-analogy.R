tissue_name <- "Pancreas"
word_vectors <- readRDS(file.path("output", paste0(tissue_name, "_embedding.rds")))

# Find gene analogies with the embedding

# BRCA1 - BRCA2 = [?] - TP53

gene_unknown <-
  word_vectors["BRCA1", , drop = FALSE] -
  word_vectors["BRCA2", , drop = FALSE] +
  word_vectors["TP53", , drop = FALSE]

text2vec::sim2(x = word_vectors, y = gene_unknown, method = "cosine", norm = "l2")[, 1] |>
  sort(decreasing = TRUE) |>
  head(n = 5) |>
  round(digits = 4)

# SOX9-TP53

# <https://doi.org/10.1016/j.celrep.2020.107742>
# <https://doi.org/10.1074/jbc.M104231200>
# <https://doi.org/10.1158/1078-0432.CCR-19-0098>
