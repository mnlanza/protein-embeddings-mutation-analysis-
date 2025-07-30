# Layer Embedding Analysis Pipeline

This repository contains scripts and tools for analyzing evo2 embeddings across different layers, with a focus on comparing embeddings between wild-type and mutated sequences.

## Data Preparation

### Input Data Requirements

1. Human Contigs Source FASTA (`input/human_contigs_src.fasta`)
   - Contains the original (non-mutated) sequences from human contigs
   - Each sequence includes the source/wild-type codons
   - Generated from the original human genomic data

2. Codon Table (`input/codon_table`)
   - Reference table (codon table 11) for codon-amino acid mappings for bacteria

3. Updated Data TSV (`input/updated_data.tsv`)
   - Contains contig information with columns:
     - aid: Subject ID
     - gene: Gene identifier
     - contig: Contig identifier
     - start: Gene start position within the contig
     - end: Gene end position within the contig
     - strand: DNA strand (+ or -)
     - src_codon: Source/wild-type codon
     - tgt_codon: Target/mutated codon
     - flipped: Strand orientation flag (True if the source codon and target codon got flipped when constructing contig)
     - mut_pos: Mutation amino acid position inside the gene

## Pipeline Overview

### 1. Running the Analysis (`run_embed_mut_analysis.sh`)

The main script orchestrates the analysis for multiple sequences:

```bash
# Run analysis for all AIDs at layer 28
./scripts/run_embed_mut_analysis.sh 28

# Run analysis for specific AID at layer 28
./scripts/run_embed_mut_analysis.sh 28 baa
```

Parameters:
- Layer number (required)
- AID filter (optional)

### 2. Codon Analysis Process (`analyze_codon_pos_embed.sh`)

For each mutation in the dataset:

1. **Sequence Generation**
   - Uses `gen_for_codons_variants.py` to generate sequence variants
   - Extracts sequence context around mutation site
   - Parameters:
     - Left margin: 2000bp (default)
     - Right margin: 1000bp (default)

2. **Embedding Generation**
   - Submits sequences to language model via `evo_gcp`
   - Generates embeddings for specified layer
   - Downloads results to `jobs/` directory

3. **Visualization**
   - Uses R scripts to generate plots
   - Compares embeddings between variants
   - Outputs plots to `figures/` directory

## Directory Structure

```
.
├── input/
│   ├── codon_table
│   ├── human_contigs_src.fasta
│   └── updated_data.tsv
├── jobs/
│   └── {aid}-{seq_id}-{pos}-layer{layer}/
│       └── output/
├── output/
│   └── {aid}/
│       └── {layer}/
│           ├── query_{seq_id}_{pos}.fasta
│           └── query_{seq_id}_{pos}.tab
└── figures/
    └── {aid}/
        └── {layer}/
            ├── cosine_sim_*.pdf
            └── euclidean_dist_*.pdf
```

## Output Files

### 1. FASTA Files
- Location: `output/{aid}/{layer}/query_{seq_id}_{pos}.fasta`
- Contains sequence variants for analysis (64 total for each codon position)

### 2. Plot Information Tables
- Location: `output/{aid}/{layer}/query_{seq_id}_{pos}.tab`
- Contains metadata for plotting

### 3. Embedding Files
- Location: `jobs/{aid}-{seq_id}-{pos}-layer{layer}/output/`
- Contains evo2 embeddings

### 4. Visualization Plots
- Location: `figures/{aid}/{layer}/`
- Types:
  - Cosine similarity plots (`cosine_sim_*.pdf`)
  - Euclidean distance plots (`euclidean_dist_*.pdf`)

## Requirements

Dependencies are listed in `requirements.txt`. The pipeline requires:
- Python 3
- R with required plotting libraries
- Access to `evo_gcp` service for embeddings generation
