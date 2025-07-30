#!/bin/bash
set -euo pipefail

# Arguments: layer (required), aid (optional)
layer="${1:-}"
aid_filter="${2:-}"
# Convert aid_filter to lowercase if provided
if [ -n "$aid_filter" ]; then
    aid_filter=$(echo "$aid_filter" | tr '[:upper:]' '[:lower:]')
fi

# Check if layer parameter is provided
if [ -z "$layer" ]; then
    echo "Usage: $0 LAYER [AID]"
    echo "Examples:"
    echo "  $0 28           # Run all aids with layer 28"
    echo "  $0 28 baa      # Run only baa with layer 28"
    exit 1
fi

# Debug output
echo "Debug: In run_embed_mut_analysis.sh"
echo "layer: $layer"
echo "aid_filter: $aid_filter"

# Create all required directories upfront
mkdir -p input output jobs figures

# Create aid-specific directories from TSV
while IFS=$'\t' read -r aid _ _ _ _ _ _ _ _ _ || [ -n "$aid" ]; do
    if [ "$aid" = "aid" ]; then continue; fi
    aid_lower="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
    mkdir -p "output/${aid_lower}" "figures/${aid_lower}"
done < "input/updated_data.tsv"

# Shared parameters
left_margin=2000
right_margin=1000
input_fasta="input/human_contigs_src.fasta"

# Process mutations
while IFS=$'\t' read -r aid gene contig start end strand flipped src_codon tgt_codon mut_pos || [ -n "$contig" ]; do
    # Skip header line
    if [ "$aid" = "aid" ]; then continue; fi
    
    # Convert aid to lowercase for comparison
    aid_lower="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
    
    # Skip if aid filter is set and doesn't match
    if [ -n "$aid_filter" ] && [ "$aid_lower" != "$aid_filter" ]; then continue; fi

    echo "Running $contig with position $mut_pos and gene range $start-$end (aid: $aid_lower, layer: $layer)"
    
    echo "  → Running position $mut_pos"
    ./scripts/analyze_codon_pos_embed.sh "$mut_pos" "$start" "$end" "$contig" "$aid_lower" "$input_fasta" "$layer" "$tgt_codon" "$left_margin" "$right_margin"

done < "input/updated_data.tsv"