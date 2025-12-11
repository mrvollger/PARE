# PARE - Paralog And allelic Regulatory Elements 

PARE is a Snakemake workflow for identifying paralogs and allelic variants of regulatory elements (REs) across genome assemblies. The pipeline extracts regulatory element sequences from haplotype-resolved assemblies, performs all-vs-all alignments to find paralogs and allelic variants, and classifies alignments based on their genomic context.

## Features

- **Extract and align regulatory elements** across multiple haplotype-resolved assemblies
- **Classify alignments** as self, paralog, allelic, ortholog, or unknown
- **Build RE networks** with Union-Find clustering for community detection
- **Export to GraphML** for visualization in Cytoscape, Gephi, etc.
- **Flexible configuration** with sensible defaults

## Installation

Install [pixi](https://pixi.sh/latest/) for environment management:

```bash
pixi install
```

## Quick Start

```bash
# Run with default config
pixi run snakemake -c4

# Run with custom config
pixi run snakemake -c4 --configfile test_data/config.yaml

# Dry run to see planned jobs
pixi run snakemake -n
```

## Configuration

### Minimal Configuration

Only one parameter is required in `config/config.yaml`:

```yaml
sample_table: "config/samples.tsv"
```

### Sample Table Format

`config/samples.tsv` should be a tab-separated file with these columns:

| sample | haplotype | assembly | re_bed |
|--------|-----------|----------|--------|
| HG002 | hap1 | path/to/HG002.1.fasta.gz | path/to/HG002.1.re.bed.gz |
| HG002 | hap2 | path/to/HG002.2.fasta.gz | path/to/HG002.2.re.bed.gz |

**Required columns:**
- `sample`: Sample identifier
- `haplotype`: Haplotype identifier (e.g., hap1, hap2)
- `assembly`: Path to assembly FASTA (can be gzipped)
- `re_bed`: Path to regulatory element BED file (can be gzipped)

### Optional Configuration

All parameters have defaults in `workflow/rules/common.smk`:

```yaml
output_dir: "results"          # Default: results
temp_dir: "temp"               # Default: temp

alignment:
  min_identity: 80             # Minimum alignment identity (%)
  min_coverage: 80             # Minimum query coverage (%)
  tool: "minimap2"             # Alignment tool
  threads: 8                   # Threads per alignment job
  minimap_preset: "asm20"      # Minimap2 preset
  max_secondary: 100           # Maximum secondary alignments
  secondary_threshold: 0.8     # Secondary alignment score threshold
  kmer_size: 15                # K-mer size
  window_size: 10              # Window size

clustering:
  allelic_distance: 10000      # Max distance for allelic calls (bp)
  reciprocal_overlap: 0.5      # Reciprocal overlap threshold
  method: "graph"              # Clustering method

filtering:
  min_re_length: 50            # Minimum RE length (bp)
  max_re_length: 10000         # Maximum RE length (bp)
  skip_repeats: false          # Filter against segmental duplications
  slop: 2000                   # Flanking sequence to add (bp)
```

## Pipeline Overview

### 1. Extract Sequences

- Filters RE BED files by length
- Adds flanking sequence (slop) for better alignment
- Extracts sequences from assemblies

### 2. All-vs-All Alignment

- Aligns all REs to each assembly using minimap2
- Reports multiple secondary alignments to find paralogs
- Adjusts coordinates back to original assembly space
- Trims alignments to exact RE boundaries

### 3. Filter and Classify

Alignments are classified into 5 categories:

| Classification | Description |
|----------------|-------------|
| **self** | RE aligns to its own genomic location |
| **paralog** | RE aligns elsewhere within the same haplotype |
| **allelic** | RE aligns between haplotypes of the same sample |
| **ortholog** | RE aligns between different samples |
| **unknown** | Sample/haplotype could not be determined |

### 4. Build Graph

- Creates an RE-centric network (nodes = REs, edges = alignments)
- Clusters related REs using Union-Find
- Exports GraphML for visualization

## Output Files

```
results/
├── filtered_alignments/
│   └── all_REs_vs_{sample_id}.paf    # Filtered alignments per assembly
└── graphs/
    ├── re_graph.graphml              # Network graph for visualization
    ├── re_clusters.tsv               # Cluster assignments
    ├── annotated.paf                 # PAF with cluster annotations
    └── annotated_res.bed             # REs with cluster assignments
```

## PAF Output Tags

The pipeline adds custom tags to PAF output:

| Tag | Description |
|-----|-------------|
| `id:Z:` | Original RE location (chrom_start_end) |
| `ct:Z:` | Classification (self/paralog/allelic/ortholog/unknown) |
| `ov:i:` | Overlap between query and target positions |
| `qs:Z:` / `ts:Z:` | Query/target sample ID |
| `qh:Z:` / `th:Z:` | Query/target haplotype |

## Supported Contig Naming

The pipeline automatically parses sample and haplotype from contig names:

- `HG00097#1#CM094075.1` → sample: HG00097, haplotype: hap1
- `HG00097#2#CM094075.1` → sample: HG00097, haplotype: hap2
- `chr1_PATERNAL` → sample: HG002, haplotype: hap1
- `chr1_MATERNAL` → sample: HG002, haplotype: hap2

## Test Data

The `test_data/` directory contains example data with two samples (HG002, HG01123) to test all classification types.

```bash
pixi run snakemake -c4 --configfile test_data/config.yaml
```

## Repository Structure

```
PARE/
├── config/
│   ├── config.yaml           # Main configuration
│   └── samples.tsv           # Sample table
├── workflow/
│   ├── Snakefile             # Main workflow entry
│   ├── rules/
│   │   ├── common.smk        # Common functions, defaults
│   │   ├── extract_sequences.smk
│   │   ├── align.smk
│   │   └── graph.smk
│   ├── scripts/
│   │   ├── paf.py            # PAF parser class
│   │   ├── filter_paf.py     # Filter and classify
│   │   ├── adjust_paf_for_slop.py
│   │   ├── paf_to_graph.py   # Graph generation
│   │   └── split_fasta.py
│   ├── envs/
│   │   ├── python.yml
│   │   └── bfx.yml
│   └── profiles/default/
│       └── config.yaml
└── test_data/
```

## Dependencies

Managed via pixi:

- **Python**: ≥3.9, pandas, numpy, networkx
- **Bioinformatics**: minimap2 2.30, samtools 1.21, bedtools, rustybam 0.1.34

## Troubleshooting

**"Chromosome not found in FAI files"**
- Ensure assemblies have `.fai` index files. Generate with `samtools faidx`.

**All alignments classified as "unknown"**
- Check contig naming conventions match supported patterns
- Add support for your naming scheme in `workflow/scripts/paf.py`

**Low alignment counts**
- Adjust `min_identity` and `min_coverage` thresholds
- Check input RE quality and length filters

## Citation

If you use PARE in your research, please cite:
(Citation information to be added)

## License

(License information to be added)