#!/bin/bash

# Exit on error, undefined var, or pipeline fail
set -euo pipefail

# === Usage Check ===
if [ $# -lt 8 ]; then
  echo "Usage: $0 POS GENE_START GENE_END SEQ_ID AID INPUT_FASTA LAYER TGT_CODON [left_margin=2000] [right_margin=1000]"
  exit 1
fi

# === Positional and Optional Arguments ===
pos="$1"
gene_start="$2"
gene_end="$3"
seq_id="$4"
aid="$5"
input_fasta="$6"
layer="$7"
tgt_codon="$8"
left_margin="${9:-2000}"
right_margin="${10:-1000}"

# Debug output
echo "Debug: Checking parameters..."
echo "pos: $pos"
echo "gene_start: $gene_start"
echo "gene_end: $gene_end"
echo "seq_id: $seq_id"
echo "aid: $aid"
echo "input_fasta: $input_fasta"
echo "layer: $layer"
echo "tgt_codon: $tgt_codon"
echo "left_margin: $left_margin"
echo "right_margin: $right_margin"

# Check required files exist
if [ ! -f "$input_fasta" ]; then
  echo "Error: Input FASTA file $input_fasta not found"
  exit 1
fi

if [ ! -f "input/codon_table" ]; then
  echo "Error: input/codon_table not found"
  exit 1
fi

# Create required directories
mkdir -p output jobs figures

# Define job ID first (moved up)
job_id="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
job_version="$(echo "${seq_id//_/-}-${pos}-layer${layer}" | tr '[:upper:]' '[:lower:]')"

# Create aid-specific directories
mkdir -p "output/${job_id}/${layer}" "figures/${job_id}/${layer}"

# File naming
fasta_out="output/${job_id}/${layer}/query_${seq_id}_${pos}.fasta"
plot_info_table="output/${job_id}/${layer}/query_${seq_id}_${pos}.tab"
job_dir="jobs/${job_id}-${job_version}"


# Create job directory before running anything
mkdir -p "$job_dir"

# Remove old files if they exist
echo "Removing old files if they exist..."
rm -f "$plot_info_table"
rm -f "$fasta_out"

# generate all codon variants at position 83
# change seq-id for each population
python3 scripts/gen_for_codons_variants.py \
  --fasta "$input_fasta" \
  --codon-table input/codon_table \
  --aa-coord "$pos" \
  --seq-id "$seq_id" \
  --output-fasta "$fasta_out" \
  --plot_info_table "$plot_info_table" \
  --left-margin "$left_margin" \
  --right-margin "$right_margin" \
  --gene-start "$gene_start" \
  --gene-end "$gene_end" \
  --aid "$aid" \
  --layer "$layer" \
  --mut-codon "$tgt_codon"

# submit job pre_norm test
evo_gcp submit --job "$job_id" \
  --output_type embedding \
  --input_fasta "$(pwd)/$fasta_out" \
  --job_version "$job_version" \
  --embedding_layers "blocks.${layer}.pre_norm" \
  --wait

# download results
evo_gcp download --job "$job_id" \
  --job_version "$job_version" \
  --jobs_dir "$(pwd)/jobs"
  
# Plot embedding differences (use heredoc to avoid line-continuation/quoting issues)
Rscript - <<RSCRIPT
source('scripts/plot_layer.R')
plot_layer_diff(
  layer_data_tsv='${plot_info_table}',
  output_dir='figures/${aid}/${layer}',
  embed_dir='$(pwd)/jobs/${job_id}-${job_version}'
)
RSCRIPT

# Clean up job directory to free up space
echo "Cleaning up job directory: $job_dir"
if [ -d "$job_dir" ]; then
    echo "  → Removing directory and all its contents..."
    rm -rf "$job_dir"
    if [ $? -eq 0 ]; then
        echo "  ✓ Cleanup successful"
    else
        echo "  ! Warning: Could not fully remove directory"
    fi
else
    echo "  ! Directory does not exist, nothing to clean"
fi

# Note: Plots and data are preserved in figures/${aid}/${layer} and output/${aid}/${layer}