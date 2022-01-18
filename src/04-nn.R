tissue_name <- "Pancreas"
word_vectors <- readRDS(file.path("output", paste0(tissue_name, "_embedding.rds")))

interest <- tibble::tibble(term = c("EGFR", "TP53", "PTEN", "KRAS"))

# Use the embedding to find nearest-neighbor genes
result <- dplyr::mutate(
  interest,
  # Create similarity matrix
  cos_sim = purrr::map(
    term,
    function(x) {
      text2vec::sim2(word_vectors, y = word_vectors[x, , drop = FALSE], method = "cosine", norm = "l2")
    }
  ),
  # Convert matrix to tibble
  cos_sim_tbl = purrr::map(
    cos_sim,
    function(x) {
      tibble::tibble(term = row.names(x), distance = x[, 1]) |>
        dplyr::arrange(desc(distance))
    }
  ), # Sort distance
  # Find neighbors
  n1 = purrr::map_chr(cos_sim_tbl, function(x) x[[2, 1]]),
  n2 = purrr::map_chr(cos_sim_tbl, function(x) x[[3, 1]]),
  n3 = purrr::map_chr(cos_sim_tbl, function(x) x[[4, 1]]),
  n4 = purrr::map_chr(cos_sim_tbl, function(x) x[[5, 1]]),
  n5 = purrr::map_chr(cos_sim_tbl, function(x) x[[6, 1]])
)

# Print table
result |>
  dplyr::select(-starts_with("cos")) |>
  kableExtra::kbl(escape = FALSE) |>
  kableExtra::kable_classic(full_width = FALSE)

# EGFR: CCT6A <https://doi.org/10.1002/pmic.201800157>
# PTEN: EIF4EBP2 <https://doi.org/10.1158/1541-7786.mcr-17-0696>
# PTEN: BMPR1A <https://doi.org/10.1093/hmg/ddab094>
# KRAS: SINHCAF <https://doi.org/10.1042/BCJ20170945>
# KRAS: ETNK1 <https://doi.org/10.1002/ijc.30509>
