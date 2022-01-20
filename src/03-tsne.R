tissue_name <- "Pancreas"
word_vectors <- readRDS(file.path("output", paste0(tissue_name, "_embedding.rds")))

path_tsne <- file.path("output", paste0(tissue_name, "_tsne.rds"))

if (!file.exists(path_tsne)) {
  set.seed(42)
  tsne_out <- Rtsne::Rtsne(word_vectors, verbose = TRUE, pca = FALSE, max_iter = 5000)
  saveRDS(tsne_out, file = path_tsne)
}

tsne_out <- readRDS(path_tsne)

# Show genes in the 2D t-SNE representation
set.seed(42)
cl <- kmeans(word_vectors, centers = 15, iter.max = 20)

df <- cbind(as.data.frame(tsne_out$Y), as.factor(cl$cluster))
names(df) <- c("x", "y", "cluster")

path_tsne_plot <- file.path("output", paste0(tissue_name, "_tsne.png"))
ragg::agg_png(path_tsne_plot, width = 3840, height = 3840 * 0.618, res = 300)
ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point(ggplot2::aes(colour = cluster), alpha = 0.3) +
  cowplot::theme_minimal_grid() +
  ggsci::scale_color_d3(palette = "category20")
invisible(dev.off())

glove_dump <- readRDS(file.path("output", paste0(tissue_name, "_glove_dump.rds")))
gene_all <- rownames(glove_dump$w_i)

# Locate the genes in the distinctive cluster
idx <- intersect(which(df$x > -40 & df$x < 0), which(df$y > 60 & df$y < 80))
# Identify the cluster number
idx_cl <- df$cluster[idx] |>
  table() |>
  which.max() |>
  names() |>
  as.numeric()
# Inspect the cluster: many RNA genes
gene_all[which(df$cluster == idx_cl)] |> sort()
