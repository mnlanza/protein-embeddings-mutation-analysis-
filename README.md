## Evaluating EVO2 (single‑nucleotide LLM) on AMR mutation retention

**Goal**: Evaluate the state‑of‑the‑art, single‑nucleotide‑resolution protein LLM EVO2 on its ability to understand and retain information about key antimicrobial resistance (AMR) mutations along the sequence—so the model can, in principle, help predict functional protein differences under mutation.

**Summary**: Using a reproducible pipeline on gyrA, we compare wild‑type vs mutated sequence embeddings across layers. We observe that the mutation signal peaks at the codon and decays within ~200–500 nt (cosine angle) and ~200 nt (distance), shorter than expected for downstream functional impact. This suggests exploring mutation‑aware “steering” toward biologically relevant positions to extend the model’s retention window.

### Roles this showcases
- **ML/Applied**: embedding analysis, similarity metrics, representation diagnostics
- **SWE/Data**: CLI pipeline, job orchestration, reproducible outputs, clear I/O contracts
- **Comp Bio**: variant generation, codon models, hypothesis‑driven evaluation

---

## Results at a glance
- **Finding**: Mutation vs wild‑type difference vector direction (cosine angle) spikes at the mutation and decays within ~500 nt. Distance magnitude decays faster (~200 nt).
- **Scope**: 8 gyrA genes × 4 variants each; consistent trend across samples.
- **Implication**: Middle layers (e.g., layer 14) capture local effects; targeted steering may help the model retain biologically meaningful mutation information longer.

Key figures (example, layer 14 pre‑norm):
- `figures/baa/14/pre_norm/cosine_sim_mut.pdf` (angle between difference vectors)
- `figures/baa/14/pre_norm/euclidean_dist_mut.pdf` (distance between difference vectors)

---

## Quickstart

### Prereqs
- Python 3, R (with `ggplot2`), and the EVO2 GCP CLI
  - EVO2 GCP: `https://github.com/eitanyaffe/evo2_gcp`
- Python deps: `pip install -r requirements.txt`
- R deps (in R): `install.packages("ggplot2"); install.packages("reticulate")`

### Minimal demo (replicates example figures)
```bash
# From repo root
./run_embed_mut_analysis.sh 14 baa
```
Outputs:
- Embeddings are fetched to `jobs/` (cleaned after plotting)
- Plots saved under `figures/baa/14/pre_norm/`

Runtime: ~10 minutes on GCP per run (varies by inputs).

---

## How it works
1. **Variant generation**: `scripts/gen_for_codons_variants.py` creates 64 codon variants at a target aa position with user‑defined margins and writes a TSV descriptor.
2. **Embeddings**: `scripts/analyze_codon_pos_embed.sh` submits sequences to EVO2 via `evo_gcp` for the requested layer(s) and downloads `.npy` outputs.
3. **Visualization**: `scripts/plot_layer.R` loads `.npy` via `reticulate`, computes cosine angle and Euclidean distance per nucleotide, and saves mutation/synonymous/stop controls.
4. **Batching**: `run_embed_mut_analysis.sh` iterates rows in `input/updated_data.tsv` and orchestrates jobs per AID.

### Inputs
1. `input/human_contigs_src.fasta`: source contigs (wild‑type)
2. `input/codon_table`: codon→aa map (table 11)
3. `input/updated_data.tsv`: columns include `aid, gene, contig, start, end, strand, src_codon, tgt_codon, flipped, mut_pos`

### Outputs
- `output/{aid}/{layer}/query_{seq_id}_{pos}.fasta` (64 variants with margins)
- `output/{aid}/{layer}/query_{seq_id}_{pos}.tab` (plot metadata)
- `figures/{aid}/{layer}/{embed_type}/` PDFs for mutation, synonymous, and stop controls

---

## Repo structure
```
scripts/
  analyze_codon_pos_embed.sh   # per‑position job: variants → evo_gcp → plots
  gen_for_codons_variants.py   # builds 64 codon variants with margins
  plot_layer.R                 # R/ggplot2 visualizations via reticulate
run_embed_mut_analysis.sh      # batch runner over updated_data.tsv (layer, optional AID)
input/                         # FASTA, codon table, and TSV
output/                        # generated FASTA + metadata
figures/                       # saved plots
```

---

## What I contributed
- Designed and implemented the end‑to‑end pipeline (Python, Bash, R)
- Built the variant generator and plotting components; integrated `evo_gcp`
- Conducted experiments on gyrA AMR mutations; analyzed and summarized results
- Mentored by and collaborated with Eitan Yaffe

---

## Requirements
- Python deps in `requirements.txt`
- R: `ggplot2`, `reticulate`
- EVO2 GCP CLI configured per `evo2_gcp` documentation

---

## Public availability / usage notice
This repository is publicly visible (clone, inspect, fork). Running the full system—especially Google Cloud operations—requires:

- A Google Cloud Platform (GCP) project with billing enabled
- Appropriate IAM permissions/service account credentials to create/manage resources (VMs, containers, storage buckets, batch jobs, etc.)
- Enabling required Google APIs for your project
- Awareness that compute/storage usage will incur costs to the GCP project owner

If using in a shared or organizational environment, ensure credentials and IAM policies are properly configured, and do not expose secret keys in public.

## Next directions
- Evaluate layer‑wise retention across more genes and conditions
- Implement and test mutation‑aware steering vectors; measure retention extension
- Add quantitative significance tests and confidence bands around decay curves

---

## License and contact
- License: No license. All rights reserved (contact for permissions beyond personal review).
- Contact: marcolanza@berkeley.edu • LinkedIn: [marconlanza](https://www.linkedin.com/in/marconlanza/)

If you’re hiring for ML, SWE/Data, or comp‑bio roles and this is relevant, I’d love to chat.
