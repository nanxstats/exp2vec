tissue_name <- "Pancreas"
word_vectors <- readRDS(file.path("output", paste0(tissue_name, "_embedding.rds")))

interest <- c("EGFR", "TP53", "PTEN", "KRAS")

# Use the embedding to find nearest-neighbor genes
query_gene <- function(gene, word_vectors, method = "cosine", norm = "l2") {
  vec <- text2vec::sim2(word_vectors, y = word_vectors[gene, , drop = FALSE], method = method, norm = norm)
  df <- data.frame(gene = rownames(vec), sim = vec[, 1], stringsAsFactors = FALSE)
  df <- df[!(df$gene == gene), ]
  df <- df[order(df$sim, decreasing = TRUE), ]
  rownames(df) <- NULL
  df
}

lst <- lapply(interest, FUN = query_gene, word_vectors = word_vectors)
df <- do.call(cbind, args = lst)
names(df) <- as.vector(t(outer(interest, c("gene", "sim"), paste, sep = ".")))

# Print table
df |>
  DT::datatable() |>
  DT::formatRound(columns = seq_len(length(interest)) * 2, digits = 4)

# EGFR: CCT6A <https://doi.org/10.1002/pmic.201800157>
# PTEN: EIF4EBP2 <https://doi.org/10.1158/1541-7786.mcr-17-0696>
# PTEN: BMPR1A <https://doi.org/10.1093/hmg/ddab094>
# KRAS: SINHCAF <https://doi.org/10.1042/BCJ20170945>
# KRAS: ETNK1 <https://doi.org/10.1002/ijc.30509>
