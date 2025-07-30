#!/usr/bin/env Rscript
library(ggplot2)
#Compute Euclidean distance
compute_euclidean_distances <- function(ref_mat, mut_mat) {
  euclidean_distances <- sapply(1:nrow(ref_mat), function(i) {
    src_i <- ref_mat[i, ]
    tgt_i <- mut_mat[i, ]
    sqrt(sum((src_i - tgt_i)^2))
  })
  return(euclidean_distances)
}

compute_log_invert_cosine_sim <- function(ref_mat, mut_mat) {
#Compute cosine similarity
log_invert_cosine_sim <- sapply(1:nrow(ref_mat), function(i) {
  src_i <- ref_mat[i, ] 
  tgt_i <- mut_mat[i, ] 
  x <- sum(src_i * tgt_i) / (sqrt(sum(src_i^2)) * sqrt(sum(tgt_i^2)))
  x <- min(x, 0.999999)  # Prevent log2(0)
  log2(1 - x)
})
  return(log_invert_cosine_sim)
}

# Helper function to get descriptive mutation type
get_mutation_type <- function(var_name) {
  type <- substr(var_name, 1, 3)
  if (type == "mut") {
    return("AMR Actual Mutation")
  } else if (type == "syn") {
    return("Synonymous Mutation")
  } else if (type == "sto") {
    return("Stop Codon Mutation")
  } else {
    return("Unknown")
  }
}

plot_cosine_similarity <- function(log_invert_cosine_sim, layer, highlight = NULL, 
                                xlim_bounds = NULL, var_name = "") {
  mutation_type <- get_mutation_type(var_name)
  df <- data.frame(pos = 1:length(log_invert_cosine_sim), cos_sim = log_invert_cosine_sim)
  p <- ggplot(df, aes(x = pos, y = cos_sim)) +
    geom_point(color = "darkred", alpha = 0.6, size = 1.2) +
    labs(title = paste0("Log Inverse Cosine Similarity: ", mutation_type, " vs Reference (Layer ", layer, ")"),
         x = "Nucleotide Position",
         y = "Log Inverse Cosine Similarity") +
    theme_minimal()
    
  # Optional: highlight positions with different colors and add to legend
  if (!is.null(highlight)) {
    # Create data frames for each highlight point with labels
    mut_pos_df <- data.frame(
      pos = highlight[1],
      cos_sim = df$cos_sim[highlight[1]],
      label = "Mutation Position"
    )
    gene_end_df <- data.frame(
      pos = highlight[2],
      cos_sim = df$cos_sim[highlight[2]],
      label = "Gene End"
    )
    
    # Add points with different colors and include in legend
    p <- p + 
      geom_point(data = mut_pos_df, aes(x = pos, y = cos_sim, color = label), size = 2.5) +
      geom_point(data = gene_end_df, aes(x = pos, y = cos_sim, color = label), size = 2.5) +
      scale_color_manual(values = c("Gene End" = "green", "Mutation Position" = "blue")) +
      theme(legend.position = "bottom",
            legend.title = element_blank())
  }
  
  # Optional: x-axis limits
  if (!is.null(xlim_bounds)) {
    p <- p + coord_cartesian(xlim = xlim_bounds)
  }
  return(p)
}

plot_euclidean_distance <- function(euclidean_distances, layer, highlight = NULL, 
                                xlim_bounds = NULL, var_name = "") {
  mutation_type <- get_mutation_type(var_name)
  df <- data.frame(pos = 1:length(euclidean_distances), euclidean_dist = euclidean_distances)
  p <- ggplot(df, aes(x = pos, y = euclidean_dist)) +
    geom_point(color = "darkred", alpha = 0.6, size = 1.2) +
    labs(title = paste0("Euclidean Distance: ", mutation_type, " vs Reference (Layer ", layer, ")"),
         x = "Nucleotide Position",
         y = "Euclidean Distance") +
    theme_minimal()
    
  # Optional: highlight positions with different colors and add to legend
  if (!is.null(highlight)) {
    # Create data frames for each highlight point with labels
    mut_pos_df <- data.frame(
      pos = highlight[1],
      euclidean_dist = df$euclidean_dist[highlight[1]],
      label = "Mutation Position"
    )
    gene_end_df <- data.frame(
      pos = highlight[2],
      euclidean_dist = df$euclidean_dist[highlight[2]],
      label = "Gene End"
    )
    
    # Add points with different colors and include in legend
    p <- p + 
      geom_point(data = mut_pos_df, aes(x = pos, y = euclidean_dist, color = label), size = 2.5) +
      geom_point(data = gene_end_df, aes(x = pos, y = euclidean_dist, color = label), size = 2.5) +
      scale_color_manual(values = c("Gene End" = "green", "Mutation Position" = "blue")) +
      theme(legend.position = "bottom",
            legend.title = element_blank())
  }
  
  if (!is.null(xlim_bounds)) {
    p <- p + coord_cartesian(xlim = xlim_bounds)
  }
  return(p)
}

plot_layer_diff <- function(layer_data_tsv, output_dir, embed_dir) {
  # Make output_dir absolute if it isn't already
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Read TSV with colClasses to ensure proper types
  df <- read.delim(layer_data_tsv, header = TRUE, 
                   colClasses = c(coord="integer",
                                codon="character",
                                aa="character",
                                left_margin="integer",
                                layer="integer",
                                start="integer",
                                end="integer",
                                aid="character",
                                mut_codon="character",
                                mut_aa="character",
                                synon_cod="character"))
  
  # Debug: print column names and first row
  print("TSV columns:")
  print(colnames(df))
  print("First row values:")
  print(df[1,])
  
  aa_coord <- df$coord[1] # for embedding file name
  codon <- df$codon[1] # for embedding file name
  amino_acid <- df$aa[1] # for embedding file name
  left_margin <- df$left_margin[1] # for highlight/plot_start
  layer <- df$layer[1]
  start <- df$start[1] # for highlight/plot_start
  end <- df$end[1] # for highlight
  aid <- df$aid[1] # for output dir
  mut_codon <- df$mut_codon[1] # for mutated strand
  mut_aa <- df$mut_aa[1] # for mutated strand
  syn_codon <- df$synon_cod[1] # for synonymous strand 
  
  # Debug: print values used for filename
  print("Values for filename construction:")
  print(paste("aa_coord:", aa_coord))
  print(paste("amino_acid:", amino_acid))
  print(paste("codon:", codon))
  print(paste("mut_codon:", mut_codon))
  print(paste("mut_aa:", mut_aa))
  
  library(reticulate)
  np <- import("numpy")

  # Debug: print full path being attempted
  print("Attempting to load file:")

  src_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_", amino_acid, "_", codon, "_embeddings_blocks_", layer, "_mlp_l3.npy"))
  stop_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_Z_TAG_embeddings_blocks_", layer, "_mlp_l3.npy"))
  mut_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_", mut_aa, "_", mut_codon, "_embeddings_blocks_", layer, "_mlp_l3.npy"))
  syn_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_", amino_acid, "_", syn_codon, "_embeddings_blocks_", layer, "_mlp_l3.npy"))

  ref_mat <- py_to_r(src_df)[1, , ]
  stop_mat <- py_to_r(stop_df)[1, , ]
  mut_mat <- py_to_r(mut_df)[1, , ]
  syn_mat <- py_to_r(syn_df)[1, , ]

  # Compute distances
  mut_diff_dist <- compute_euclidean_distances(ref_mat, mut_mat)
  mut_diff_cosine_sim <- compute_log_invert_cosine_sim(ref_mat, mut_mat)

  syn_diff_dist <- compute_euclidean_distances(ref_mat, syn_mat)
  syn_diff_cosine_sim <- compute_log_invert_cosine_sim(ref_mat, syn_mat)

  stop_diff_dist <- compute_euclidean_distances(ref_mat, stop_mat)
  stop_diff_cosine_sim <- compute_log_invert_cosine_sim(ref_mat, stop_mat)

  # Plotting
  print("Plotting cosine similarity and euclidean distance...")

  if (start <= left_margin) {
    left_margin <- start-1
  }
  mut_nuc_pos <- (aa_coord - 1) * 3
  start_mut_hl <- left_margin + mut_nuc_pos
  gene_len <- end-start+1
  end_mut_hl <- gene_len + left_margin
  #print(c(start_mut_hl, end_mut_hl))

  codon_pos_adj <- 1

  suppressWarnings({
  p1 <- plot_cosine_similarity(mut_diff_cosine_sim, layer = layer, 
                             highlight = c(start_mut_hl+codon_pos_adj, end_mut_hl),
                             var_name = "mut_diff_cosine_sim")
  p2 <- plot_euclidean_distance(mut_diff_dist, layer = layer, 
                               highlight = c(start_mut_hl+codon_pos_adj, end_mut_hl),
                               var_name = "mut_diff_dist")

  p3 <- plot_cosine_similarity(syn_diff_cosine_sim, layer = layer, 
                             highlight = c(start_mut_hl+codon_pos_adj, end_mut_hl),
                             var_name = "syn_diff_cosine_sim")
  p4 <- plot_euclidean_distance(syn_diff_dist, layer = layer, 
                               highlight = c(start_mut_hl+codon_pos_adj, end_mut_hl),
                               var_name = "syn_diff_dist")

  p5 <- plot_cosine_similarity(stop_diff_cosine_sim, layer = layer, 
                             highlight = c(start_mut_hl+codon_pos_adj, end_mut_hl),
                             var_name = "stop_diff_cosine_sim")
  p6 <- plot_euclidean_distance(stop_diff_dist, layer = layer, 
                               highlight = c(start_mut_hl+codon_pos_adj, end_mut_hl),
                               var_name = "stop_diff_dist")
  })

  # Use file.path to construct proper paths
  ggsave(file.path(output_dir, "cosine_sim_mut.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(output_dir, "euclidean_dist_mut.pdf"), p2, width = 8, height = 6)

  ggsave(file.path(output_dir, "cosine_sim_syn.pdf"), p3, width = 8, height = 6)
  ggsave(file.path(output_dir, "euclidean_dist_syn.pdf"), p4, width = 8, height = 6)

  ggsave(file.path(output_dir, "cosine_sim_stop.pdf"), p5, width = 8, height = 6)
  ggsave(file.path(output_dir, "euclidean_dist_stop.pdf"), p6, width = 8, height = 6)
  print("Done!")
}