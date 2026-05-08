# Panscan `complex` Module — Tutorial

## Overview

The `complex` module identifies and visualizes complex structural variant loci in a pangenome graph. It detects regions with large SV bubbles, extracts subgraphs using `gfabase`, colors nodes by gene identity using `GraphAligner` alignments, and produces per-sample haplotype walk visualizations using `Bandage`.

---

## Defining a Complex SV Locus

As per the APR ([Arab Pangenome Reference](https://www.nature.com/articles/s41467-025-61645-w)) paper definition, a complex SV locus is defined as a genomic window containing **at least one top-level graph bubble (snarl) with a minimum allele size of 10,000bp and a minimum of 5 alternate alleles**. This reflects regions where multiple large, structurally distinct haplotypes co-exist in the population.

For reproducibility and direct comparison with HPRC-style analyses (following the methodology of the [HPRC draft pangenome paper](https://www.nature.com/articles/s41586-023-05896-x)), users can run with relaxed settings that lower the allele size threshold to 5,000bp (selected by the HPRC authors based on manual inspection of Bandage plots) and reduce the minimum allele count to 1. This will identify a broader set of structurally complex regions and allow direct comparison with HPRC published loci. See the parameter table and example commands below.

---

## Prerequisites

All tools (`gfabase`, `GraphAligner`, `Bandage`) are bundled in the panscan Singularity image (`panscan.sif`). You do not need to install them separately.

Required inputs:

| File | Description |
|---|---|
| `pangenome_complete.gfa(.gz)` | Pangenome graph in GFA format (from Minigraph-Cactus) |
| `pangenome_complete.gfab` | Indexed GFA database (built from GFA, see Step 2) |
| `pangenome_complete.vcf(.gz)` | Pangenome VCF |
| `genes.fa` | Gene sequences in FASTA (from Ensembl) |
| `annotation.gff3` | Gene annotation file (from Ensembl/Gencode) |

---

## A Note on `--sep_pattern`

The separator pattern used in your pangenome's path/walk names must match the `--sep_pattern` argument. This varies depending on how the pangenome was constructed. You can inspect it directly from your GFA:

```bash
grep "^W" pangenome_complete.gfa | head -3 | awk '{print $2, $3, $4}'
```

Common patterns:

| Pangenome | Example path name | `--sep_pattern` |
|---|---|---|
| HPRC v1.1 (CHM13-based) | `CHM13#1#chr1` | `#` |
| UPR (CHM13-based) | `CHM13#0#chr1` | `#0#` |
| GRCh38-based | `GRCh38#0#chr1` | `#0#` |

Using the wrong separator will cause gene coloring and walk extraction to fail silently. Always verify before running.

---
## Trying panscan complex on UPR Phase 1 Data

To run a trial analysis, the APR Phase 1 dataset is publicly available and provides a well-characterized pangenome suitable for testing the `complex` module end-to-end.

**Primary download:** [MBRU Arab Pangenome Reference](https://www.mbru.ac.ae/the-arab-pangenome-reference/)

**Backup / Zenodo archive:** [https://zenodo.org/records/17587133](https://zenodo.org/records/17587133)

Download the pangenome GFA, VCF, and use `--sep_pattern "#0#"` and `--ref_name CHM13` for this dataset.


## Step-by-Step

### Step 1 — Decompress GFA (if needed)

```bash
if [[ ! -f "pangenome_complete.gfa" ]]; then
    gunzip -c pangenome_complete.gfa.gz > pangenome_complete.gfa
fi
```

### Step 2 — Build GFAB index

The `.gfab` file is a SQLite-based indexed version of the GFA that enables fast coordinate-based subgraph extraction.

```bash
singularity exec -B /mnt:/mnt --no-home panscan.sif \
    gfabase load \
    -o pangenome_complete.gfab \
    pangenome_complete.gfa
```

This only needs to be done once per pangenome. The `.gfab` file can be reused for all subsequent `complex` runs.

### Step 3 — Align genes to pangenome graph

`GraphAligner` aligns gene sequences to the pangenome graph to identify which graph nodes correspond to which genes. This enables gene-colored Bandage visualizations.

```bash
singularity exec -B /mnt:/mnt --no-home panscan.sif \
    GraphAligner \
    -g pangenome_complete.gfa \
    -f genes.fa \
    -a genes_vs_pangenome.gaf \
    -x vg \
    --multimap-score-fraction 0.1 \
    --threads 32
```

The output `.gaf` file maps gene names to pangenome segment IDs and is reused across all regions. This step can be slow for whole-genome pangenomes — run it once and cache the result.

### Step 4 — Decompress VCF (if needed)

```bash
if [[ ! -f "pangenome_complete.vcf" ]]; then
    gunzip -c pangenome_complete.vcf.gz > pangenome_complete.vcf
fi
```

### Step 5 — Run complex module

**UPR-definition (default, stricter):**

```bash
panscan complex \
    pangenome_complete.vcf \
    pangenome_complete.gfab \
    --gff3 annotation.gff3 \
    --gaf_file genes_vs_pangenome.gaf \
    --ref_name CHM13 \
    --sep_pattern "#0#" \
    --threads 32 \
    --plot_threads 8 \
    --plot_min_free_gb 30 \
    --min_free_gb 20
```

**HPRC-style (relaxed, for comparison with published loci):**

```bash
panscan complex \
    pangenome_complete.vcf \
    pangenome_complete.gfab \
    --gff3 annotation.gff3 \
    --gaf_file genes_vs_pangenome.gaf \
    --ref_name CHM13 \
    --sep_pattern "#0#" \
    -s 5000 \
    -a 1 \
    -n 1 \
    -l 100000 \
    --sites 1 \
    --threads 32 \
    --plot_threads 8 \
    --plot_min_free_gb 30 \
    --min_free_gb 20
```

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `-s` | 10000 | Minimum SV size (bp) to define a complex site (UPR: 10000, HPRC-style: 5000) |
| `-a` | 5 | Minimum number of ALT alleles at a site (UPR: 5, HPRC-style: 1) |
| `-n` | 1 | Minimum number of large variants per site |
| `-l` | 100000 | Window size (bp) for region definition |
| `--sites` | 1 | Minimum complex sites per window |
| `--ref_name` | CHM13 | Reference sample name as it appears in the pangenome W lines |
| `--sep_pattern` | `#0#` | Path name separator — must match your pangenome (see note above) |
| `--min_free_gb` | 20 | Minimum free RAM (GB) before launching region workers |
| `--plot_min_free_gb` | 0 | Minimum free RAM (GB) before launching each Bandage plot job |
| `--plot_threads` | auto | Number of parallel Bandage workers (see warning below) |
| `--threads` | auto | Total thread budget |

---

## ⚠️ Bandage Plotting — Important Warning

Bandage is a single-threaded application that loads the entire graph into memory for layout computation. Running many Bandage jobs in parallel causes significant I/O and memory overhead on shared HPC clusters and **will cause job crashes** if not properly controlled.

**Always set the following:**

- `--plot_min_free_gb 20` to `40` — ensures sufficient free RAM before each Bandage job launches. Use 40GB for HPRC-scale pangenomes.
- `--plot_threads` to **no more than 8–10** regardless of how many CPUs your job has allocated. More parallel Bandage workers does not improve throughput — it causes memory contention, I/O saturation, and crashes.

```bash
panscan complex \
    ... \
    --threads 64 \
    --plot_threads 8 \
    --plot_min_free_gb 30
```

The graph extraction and walk generation steps scale well with `--threads`. Only the final Bandage rendering step needs to be throttled. Keeping `--plot_threads` low will **not meaningfully affect total runtime** — panscan pipelines the produce and plot steps concurrently, so Bandage rendering overlaps with graph extraction for subsequent regions.

---

## Outputs

```
./
├── complex_regions.csv                   # Detected complex regions with coordinates
├── complex_sites_summary.csv             # Per-site summary with gene annotations
└── all_plottables/
    └── pangenome_complete.<region>.wd/
        ├── <region>.gfa                  # Reference-anchored regional subgraph
        ├── color_csvs/
        │   ├── all_colors.csv            # Node → gene color mapping (all genes)
        │   ├── pc_colors.csv             # Node → gene color mapping (protein-coding only)
        │   ├── gene_colors/              # Per-gene node color CSVs
        │   └── <sample>.hap<N>.walks.csv # Per-sample haplotype walk color CSVs
        ├── colored_gfas/
        │   ├── all_colors_<region>.gfa         # Full graph colored by all genes
        │   ├── pc_colors_<region>.gfa          # Protein-coding genes only
        │   └── <sample>.hap<N>_<region>.gfa    # Per-sample haplotype walk GFAs
        ├── bandage_images/               # Raw Bandage PNG outputs
        └── bandage_images_with_legend/   # Bandage images with gene color legend
```

---

## HPC Example (SLURM)

```bash
#!/bin/bash
#SBATCH -p your-partition
#SBATCH --cpus-per-task 64
#SBATCH --mem 300G
#SBATCH --output=run.txt
#SBATCH --error=run.err
#SBATCH --job-name=panscan-complex

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-64}"
SIF="/path/to/panscan.sif"
BIND="/mnt"

# Step 1 — decompress GFA
if [[ ! -f "pangenome_complete.gfa" ]]; then
    gunzip -c pangenome_complete.gfa.gz > pangenome_complete.gfa
fi

# Step 2 — build gfab index (once per pangenome)
if [[ ! -f "pangenome_complete.gfab" ]]; then
    singularity exec -B ${BIND}:${BIND} --no-home ${SIF} \
        gfabase load -o pangenome_complete.gfab pangenome_complete.gfa
fi

# Step 3 — align genes to graph (once per pangenome)
if [[ ! -f "genes_vs_pangenome.gaf" ]]; then
    singularity exec -B ${BIND}:${BIND} --no-home ${SIF} \
        GraphAligner \
        -g pangenome_complete.gfa \
        -f genes.fa \
        -a genes_vs_pangenome.gaf \
        -x vg --multimap-score-fraction 0.1 \
        --threads ${THREADS}
fi

# Step 4 — decompress VCF
if [[ ! -f "pangenome_complete.vcf" ]]; then
    gunzip -c pangenome_complete.vcf.gz > pangenome_complete.vcf
fi

# Step 5 — run complex module (HPRC-style relaxed settings)
panscan complex \
    pangenome_complete.vcf \
    pangenome_complete.gfab \
    --gff3 annotation.gff3 \
    --gaf_file genes_vs_pangenome.gaf \
    --ref_name CHM13 \
    --sep_pattern "#0#" \
    -s 5000 -a 1 -n 1 -l 100000 --sites 1 \
    --threads ${THREADS} \
    --plot_threads 8 \
    --plot_min_free_gb 30 \
    --min_free_gb 50

echo "Done."
```

---

## Notes

**Reference name (`--ref_name`)** must match exactly how the reference sample appears in your pangenome W lines — this is the sample name (e.g. `CHM13`), not the chromosome name.

**Sample haplotype visualization** — each sample's walk through the region is extracted and colored by gene identity. Samples with missing genotypes (`GT: .|.`) at large SVs will show minimal traversal through the reference graph. This is biologically correct — it indicates the sample carries a large private structural variant not represented in the reference-anchored graph, and should be interpreted as such rather than as a pipeline failure.

**Memory** — `--min_free_gb` controls region-level parallelism (graph extraction workers), while `--plot_min_free_gb` controls Bandage-level parallelism independently. Set both appropriately for your cluster. For HPRC-scale pangenomes, `--min_free_gb 50` and `--plot_min_free_gb 30` are recommended starting points.

**⏱ GFAB Conversion Time Reference**
Building the `.gfab` index is a one-time operation per pangenome. Approximate runtimes on a standard HPC node:
- HPRC v1.1: ~15 minutes
- HPRC v2: ~1.5 hours

Plan accordingly and do not include this step in time-sensitive job allocations.
