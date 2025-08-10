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
                                xlim_bounds = NULL, var_name = "", mutation_region = NULL,
                                gene_start_pos = NULL, color_by_codon = TRUE) {
  mutation_type <- get_mutation_type(var_name)
  df <- data.frame(pos = 1:length(log_invert_cosine_sim), cos_sim = log_invert_cosine_sim)
  # Start base plot
  p <- ggplot(df, aes(x = pos, y = cos_sim))

  # Optional: shaded region for mutation codon (with legend) FIRST so points draw on top
  if (!is.null(mutation_region) && length(mutation_region) == 2) {
    region_df <- data.frame(
      xmin = mutation_region[1],
      xmax = mutation_region[2],
      ymin = -Inf,
      ymax = Inf,
      label = factor("Mutation Codon", levels = c("Mutation Codon"))
    )
    p <- p +
      geom_rect(
        data = region_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
        inherit.aes = FALSE,
        alpha = 0.22,
        color = "#6baed6",
        linewidth = 0.2
      ) +
      scale_fill_manual(values = c("Mutation Codon" = "#6baed6"), name = NULL)
  }

  # Points: optionally color by codon position relative to gene_start_pos
  if (!is.null(gene_start_pos) && isTRUE(color_by_codon)) {
    df$codon_i <- ((df$pos - gene_start_pos) %% 3) + 1
    p <- p +
      geom_point(data = df, aes(x = pos, y = cos_sim, color = factor(codon_i)), alpha = 0.6, size = 1.2) +
      scale_color_manual(name = "Codon position",
                         values = c("1" = "red", "2" = "orange", "3" = "green"))
  } else {
    p <- p + geom_point(color = "darkred", alpha = 0.6, size = 1.2)
  }

  p <- p +
    labs(title = paste0("Log Inverse Cosine Similarity: ", mutation_type, " vs Reference (Layer ", layer, ")"),
         x = "Nucleotide Position",
         y = "Log Inverse Cosine Similarity") +
    theme_minimal()
    
  # (moved earlier so it draws under points)

  # Optional: vertical lines for gene start/end (fixed colors to avoid legend conflicts)
  if (!is.null(highlight)) {
    pos_vals <- as.numeric(highlight)
    if (length(pos_vals) >= 2) {
      p <- p +
        geom_vline(xintercept = pos_vals[1], linetype = "dashed", linewidth = 0.5, color = "orange", show.legend = FALSE) +
        geom_vline(xintercept = pos_vals[2], linetype = "dashed", linewidth = 0.5, color = "green", show.legend = FALSE)
    }
  }
  
  # Optional: x-axis limits
  if (!is.null(xlim_bounds)) {
    p <- p + coord_cartesian(xlim = xlim_bounds)
  }
  return(p)
}

plot_euclidean_distance <- function(euclidean_distances, layer, highlight = NULL, 
                                xlim_bounds = NULL, var_name = "", mutation_region = NULL,
                                gene_start_pos = NULL, color_by_codon = TRUE) {
  mutation_type <- get_mutation_type(var_name)
  df <- data.frame(pos = 1:length(euclidean_distances), euclidean_dist = euclidean_distances)
  # Start base plot
  p <- ggplot(df, aes(x = pos, y = euclidean_dist))

  # Optional: shaded region for mutation codon (with legend) FIRST so points draw on top
  if (!is.null(mutation_region) && length(mutation_region) == 2) {
    region_df <- data.frame(
      xmin = mutation_region[1],
      xmax = mutation_region[2],
      ymin = -Inf,
      ymax = Inf,
      label = factor("Mutation Codon", levels = c("Mutation Codon"))
    )
    p <- p +
      geom_rect(
        data = region_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
        inherit.aes = FALSE,
        alpha = 0.22,
        color = "#6baed6",
        linewidth = 0.2
      ) +
      scale_fill_manual(values = c("Mutation Codon" = "#6baed6"), name = NULL)
  }

  # Points: optionally color by codon position relative to gene_start_pos
  if (!is.null(gene_start_pos) && isTRUE(color_by_codon)) {
    df$codon_i <- ((df$pos - gene_start_pos) %% 3) + 1
    p <- p +
      geom_point(data = df, aes(x = pos, y = euclidean_dist, color = factor(codon_i)), alpha = 0.6, size = 1.2) +
      scale_color_manual(name = "Codon position",
                         values = c("1" = "red", "2" = "orange", "3" = "green"))
  } else {
    p <- p + geom_point(color = "darkred", alpha = 0.6, size = 1.2)
  }

  p <- p +
    labs(title = paste0("Euclidean Distance: ", mutation_type, " vs Reference (Layer ", layer, ")"),
         x = "Nucleotide Position",
         y = "Euclidean Distance") +
    theme_minimal()
    
  # (moved earlier so it draws under points)

  # Optional: vertical lines for gene start/end (fixed colors to avoid legend conflicts)
  if (!is.null(highlight)) {
    pos_vals <- as.numeric(highlight)
    if (length(pos_vals) >= 2) {
      p <- p +
        geom_vline(xintercept = pos_vals[1], linetype = "dashed", linewidth = 0.5, color = "orange", show.legend = FALSE) +
        geom_vline(xintercept = pos_vals[2], linetype = "dashed", linewidth = 0.5, color = "green", show.legend = FALSE)
    }
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

  # Plotting for multiple embedding types
  print("Plotting cosine similarity and euclidean distance for embedding types...")

  # Compute positional helpers once
  if (start <= left_margin) {
    left_margin <- start-1
  }
  mut_nuc_pos <- (aa_coord - 1) * 3
  start_mut_hl <- left_margin + mut_nuc_pos
  gene_len <- end-start+1
  end_mut_hl <- gene_len + left_margin
  gene_start_hl <- left_margin + 1
  mut_codon_xmin <- start_mut_hl
  mut_codon_xmax <- start_mut_hl + 2

  for (embed_type in c("mlp_l3", "pre_norm")) {
    message(paste("Processing embed type:", embed_type))
    # Load numpy arrays with suffix based on embed_type
    src_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_", amino_acid, "_", codon, "_embeddings_blocks_", layer, "_", embed_type, ".npy"))
    stop_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_Z_TAG_embeddings_blocks_", layer, "_", embed_type, ".npy"))
    mut_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_", mut_aa, "_", mut_codon, "_embeddings_blocks_", layer, "_", embed_type, ".npy"))
    syn_df <- np$load(paste0(embed_dir, "/output/input_", aa_coord, "_", amino_acid, "_", syn_codon, "_embeddings_blocks_", layer, "_", embed_type, ".npy"))

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

    suppressWarnings({
      p1 <- plot_cosine_similarity(
        mut_diff_cosine_sim, layer = layer,
        highlight = c(gene_start_hl, end_mut_hl),
        var_name = "mut_diff_cosine_sim",
        mutation_region = c(mut_codon_xmin, mut_codon_xmax),
        gene_start_pos = gene_start_hl,
        color_by_codon = TRUE
      )
      p2 <- plot_euclidean_distance(
        mut_diff_dist, layer = layer,
        highlight = c(gene_start_hl, end_mut_hl),
        var_name = "mut_diff_dist",
        mutation_region = c(mut_codon_xmin, mut_codon_xmax),
        gene_start_pos = gene_start_hl,
        color_by_codon = TRUE
      )

      p3 <- plot_cosine_similarity(
        syn_diff_cosine_sim, layer = layer,
        highlight = c(gene_start_hl, end_mut_hl),
        var_name = "syn_diff_cosine_sim",
        mutation_region = c(mut_codon_xmin, mut_codon_xmax),
        gene_start_pos = gene_start_hl,
        color_by_codon = TRUE
      )
      p4 <- plot_euclidean_distance(
        syn_diff_dist, layer = layer,
        highlight = c(gene_start_hl, end_mut_hl),
        var_name = "syn_diff_dist",
        mutation_region = c(mut_codon_xmin, mut_codon_xmax),
        gene_start_pos = gene_start_hl,
        color_by_codon = TRUE
      )

      p5 <- plot_cosine_similarity(
        stop_diff_cosine_sim, layer = layer,
        highlight = c(gene_start_hl, end_mut_hl),
        var_name = "stop_diff_cosine_sim",
        mutation_region = c(mut_codon_xmin, mut_codon_xmax),
        gene_start_pos = gene_start_hl,
        color_by_codon = TRUE
      )
      p6 <- plot_euclidean_distance(
        stop_diff_dist, layer = layer,
        highlight = c(gene_start_hl, end_mut_hl),
        var_name = "stop_diff_dist",
        mutation_region = c(mut_codon_xmin, mut_codon_xmax),
        gene_start_pos = gene_start_hl,
        color_by_codon = TRUE
      )
    })

    # Save into subdirectory per embedding type
    out_dir <- file.path(output_dir, embed_type)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(out_dir, "cosine_sim_mut.pdf"), p1, width = 8, height = 6)
    ggsave(file.path(out_dir, "euclidean_dist_mut.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(out_dir, "cosine_sim_syn.pdf"), p3, width = 8, height = 6)
    ggsave(file.path(out_dir, "euclidean_dist_syn.pdf"), p4, width = 8, height = 6)
    ggsave(file.path(out_dir, "cosine_sim_stop.pdf"), p5, width = 8, height = 6)
    ggsave(file.path(out_dir, "euclidean_dist_stop.pdf"), p6, width = 8, height = 6)
  }
  print("Done!")
}