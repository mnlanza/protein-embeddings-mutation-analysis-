library(reticulate)

np <- import("numpy")
# default layer is 28
src_df <- np$load("input/baf/input_82_S_TCG_P_embeddings_blocks_28_mlp_l3.npy")
mut_df <- np$load("input/baf/input_82_L_TTG_P_embeddings_blocks_28_mlp_l3.npy")
# py_to_r(df)
# str(py_to_r(df))

ref_mat <- py_to_r(src_df)[1, , ]
mut_mat <- py_to_r(mut_df)[1, , ]

vec_pos1_ref <- ref_mat[1,]
vec_pos1_mut <- mut_mat[1,]

#Compute Euclidean distance
euclidean_distances <- sapply(1:nrow(ref_mat), function(i) {
  src_i <- ref_mat[i, ] 
  tgt_i <- mut_mat[i, ] 
  
  sqrt(sum((src_i - tgt_i)^2))
})
euclidean_distances

#Compute cosine similarity
log_invert_cosine_sim <- sapply(1:nrow(ref_mat), function(i) {
  src_i <- ref_mat[i, ] 
  tgt_i <- mut_mat[i, ] 
  
  x <- sum(src_i * tgt_i) / (sqrt(sum(src_i^2)) * sqrt(sum(tgt_i^2)))
  x <- min(x, 0.999999)  # Prevent log2(0)
  log2(1 - x)
})
log_invert_cosine_sim


library(ggplot2)
library(reshape2)

plot_cosine_similarity <- function(log_invert_cosine_sim, highlight = NULL, 
                                xlim_bounds = NULL) {
  df <- data.frame(pos = 1:length(log_invert_cosine_sim), cos_sim = log_invert_cosine_sim)
  p <- ggplot(df, aes(x = pos, y = cos_sim)) +
    geom_point(color = "darkred", alpha = 0.6, size = 1.2) +
    labs(title = "Log Inverse Cosine Similarity between the two vectors in layer 28 at Each Position", x = "Nucleotide Position",
         y = "Log Inverse Cosine Similarity") +
    theme_minimal()
  # Optional: highlight a specific position
  if (!is.null(highlight)) {
    p <- p + 
      geom_point(data = df[highlight, , drop = FALSE], 
                 aes(x = pos, y = cos_sim),
                 color = "blue", size = 2.5)
  }
  # Optional: x-axis limits
  if (!is.null(xlim_bounds)) {
    p <- p + xlim(xlim_bounds[1], xlim_bounds[2])
  }
  return(p)
}

plot_cosine_similarity(log_invert_cosine_sim, highlight = 269, xlim_bounds = c(260, 360))
plot_cosine_similarity(log_invert_cosine_sim, highlight = 269)

plot_euclidean_distance <- function(euclidean_distances, highlight = NULL, 
                                xlim_bounds = NULL) {
  df <- data.frame(pos = 1:length(euclidean_distances), euclidean_dist = euclidean_distances)
  p <- ggplot(df, aes(x = pos, y = euclidean_dist)) +
    geom_point(color = "darkred", alpha = 0.6, size = 1.2) +
    labs(title = "Euclidean Distance between the two vectors in layer 28 at Each Position", x = "Nucleotide Position",
         y = "Euclidean Distance") +
    theme_minimal()
  if (!is.null(highlight)) {
    p <- p + 
      geom_point(data = df[highlight, , drop = FALSE], 
                 aes(x = pos, y = euclidean_dist),
                 color = "blue", size = 2.5)
  }
  if (!is.null(xlim_bounds)) {
    p <- p + xlim(xlim_bounds[1], xlim_bounds[2])
  }
  return(p)
}
plot_euclidean_distance(euclidean_distances, highlight = 269, xlim_bounds = c(260, 360))
plot_euclidean_distance(euclidean_distances, highlight = 269)

plot_layer_diff <- function(layer_data_tsv) {

}
    # desired tsv contains aid, positio, amino_acid, codon, layer