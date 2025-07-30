#!/bin/bash

# Exit on error, undefined var, or pipeline fail
set -euo pipefail

# === Usage Check ===
if [ $# -lt 8 ]; then
  echo "Usage: $0 POS GENE_START GENE_END SEQ_ID AID INPUT_FASTA LAYER [left_margin=2000] [right_margin=1000]"
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
compare_out="output/${job_id}/${layer}/compare_strands_${seq_id}_${pos}.tab"

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

# submit job
evo_gcp submit --job "$job_id" \
  --output_type embedding \
  --input_fasta "$(pwd)/$fasta_out" \
  --job_version "$job_version" \
  --embedding_layers "blocks.${layer}.mlp.l3" \
  --wait
#   --embedding_layers "$layer" \

# download results
evo_gcp download --job "$job_id" \
  --job_version "$job_version" \
  --jobs_dir "$(pwd)/jobs"
  
# Plot embedding differences
Rscript -e "
source('scripts/plot_layer.R')
plot_layer_diff(
  layer_data_tsv='$plot_info_table',
  output_dir='figures/${aid}/${layer}',
  embed_dir='$(pwd)/jobs/${job_id}-${job_version}'
)
"

# jobs are at $(pwd)/jobs/${job_id}-${job_version} job/aid(lowercase)/layer/"${seq_id//_/-}-${pos}-layer${layer}