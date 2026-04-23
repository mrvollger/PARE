# Union Peaks Tables

## Thresholds used
- MIN_FIRE_COV: 4 (minimum FIRE coverage to call a peak)
- MIN_FRAC_ACC: 0.2 (minimum fraction accessibility to call a peak)
- DIFF_THRESHOLD: 0.25 (minimum H1-H2 accessibility difference for haplotype-selective)
- MIN_HAP_COV: 10 (minimum coverage on both haplotypes to test selectivity)

## Output files
All files are bgzipped BED files with a `#chrom` header line.
Each file contains the same columns but with different leading coordinates
(renamed to chrom/start/end), filtered to rows where those coordinates exist.

### All sites (every sample x consensus peak that passes in at least one sample)
- `union-peaks-cons.bed.gz`: All rows, leading coords are consensus (graph) coordinates
- `union-peaks-chm13.bed.gz`: Subset where consensus coords are bare chr* names (T2T-CHM13 reference paths)
- `union-peaks-asm.bed.gz`: Rows with assembly coordinates (not NA)
- `union-peaks-hg38.bed.gz`: Rows with hg38 liftover coordinates (not NA)

### Peak-only (is_peak == TRUE)
Filtered to rows where the consensus peak passes thresholds in that specific sample
(fire_coverage >= 4 AND fire_cov_OV_cov >= 0.2).
Each row represents a called peak in a specific sample/haplotype.
- `peaks-cons.bed.gz`: Peaks only, consensus coordinates
- `peaks-chm13.bed.gz`: Peaks only, CHM13 coordinates
- `peaks-asm.bed.gz`: Peaks only, assembly coordinates
- `peaks-hg38.bed.gz`: Peaks only, hg38 coordinates

## Columns (in output order)

| # | Column | Description |
|---|--------|-------------|
| 1 | chrom | Leading coordinate chromosome (renamed from cons_chr, asm_chr, or hg38_chr depending on file) |
| 2 | start | Leading coordinate start |
| 3 | end | Leading coordinate end |
| 4 | ind_key | Individual+platform key (sample_id without haplotype suffix) |
| 5 | consensus_peak_id | Unique ID for the consensus peak in the graph |
| 6 | sample_id | Full sample identifier (e.g. HG00438_PS00991_ONT_HG00438_1) |
| 7 | coverage | Total fiber-seq coverage at this site in this sample |
| 8 | fire_coverage | Number of fibers showing accessible chromatin (FIRE) |
| 9 | score | Peak score from the caller |
| 10 | fire_cov_OV_cov | Fraction of fibers that are accessible (fire_coverage / coverage) |
| 11 | asm | Assembly identifier |
| 12 | asm_chr | Assembly contig name (NA if no assembly mapping) |
| 13 | asm_start | Assembly start coordinate |
| 14 | asm_end | Assembly end coordinate |
| 15 | asm_peak_id | Peak ID in assembly coordinates |
| 16 | hg38_chr | hg38 chromosome (NA if no liftover) |
| 17 | hg38_start | hg38 start coordinate |
| 18 | hg38_end | hg38 end coordinate |
| 19 | exists_in_hg38 | Whether the peak has hg38 coordinates (logical) |
| 20 | overlaps_original_called_peak | Whether this site overlaps a peak from the original per-sample analysis (logical) |
| 21 | sample_id_peak_id_TEMP_COL | Temporary join column from peak calling pipeline |
| 22 | num_overlapping_orig_peaks | Number of original per-sample peaks overlapping this consensus peak |
| 23 | orig_peaks_overlapping_consensus | IDs of original peaks overlapping this consensus peak |
| 24 | seq | DNA sequence at the peak |
| 25 | is_peak | Whether this site passes peak thresholds (fire_coverage >= 4 AND fire_cov_OV_cov >= 0.2) |
| 26 | Individual_ID | Individual identifier (e.g. HG00438), parsed from sample_id |
| 27 | Haplotype | Haplotype number (1 or 2), parsed from sample_id |
| 28 | Platform | Sequencing platform (PacBio or ONT; HG002_HG002 samples are PacBio) |
| 29 | PS | Sample preparation identifier, parsed from sample_id (NA for HG002) |
| 30 | is_SD_cons | Whether the consensus coordinates overlap a segmental duplication |
| 31 | is_SD_asm | Whether the assembly coordinates overlap a segmental duplication |
| 32 | is_SD | Whether either coordinate set overlaps a segmental duplication |
| 33 | chrom_hg38 | hg38 chromosome name for the assembly contig (from chromosome alias table) |
| 34 | autosome | "Autosomes" or "Sex chromosome", based on chrom_hg38/asm_chr |
| 35 | is_complex | Multiple original asm peaks merged into this consensus peak |
| 36 | cons_assembly | Assembly of consensus peak (e.g. HG00438_1, HG002_1, HG002_2, CHM13) |
| 37 | haplotype_distinct | Haplotype classification (see below) |
| 38 | diff | Accessibility difference between haplotypes (fire_cov_OV_cov H1 - H2) |
| 39 | has_het_variants | Whether the peak sequences differ between haplotypes (logical) |
| 40 | local_mismatches | Number of substitutions in best local alignment of haplotype sequences |
| 41 | pct_seq_identity | Percent sequence identity between haplotype sequences (0-1) |
| 42 | fisher_p | Fisher exact test p-value for accessibility difference between haplotypes |
| 43 | fire_coverage_1 | FIRE coverage on haplotype 1 for this individual |
| 44 | fire_coverage_2 | FIRE coverage on haplotype 2 for this individual |
| 45 | coverage_1 | Total coverage on haplotype 1 for this individual |
| 46 | coverage_2 | Total coverage on haplotype 2 for this individual |

### haplotype_distinct categories
- **chrX**: Peak on chrX (excluded from autosomal haplotype analysis)
- **chrY**: Peak on chrY (excluded from autosomal haplotype analysis)
- **Complex**: Multiple original assembly peaks merged into this consensus peak on either haplotype
- **Genetically haplotype-specific**: Peak on one haplotype, no assembly mapping on the other
- **Haplotype-selective**: |diff| >= 0.25 AND Fisher p <= 0.05, with het variants between haplotype sequences
- **Epigenetically haplotype-selective**: |diff| >= 0.25 AND Fisher p <= 0.05, with identical haplotype sequences
- **Shared**: Sufficient coverage but does not meet both diff and Fisher thresholds
- **Untestable**: Insufficient coverage on one or both haplotypes to determine
- **NA**: Peak not called on either haplotype for this individual
