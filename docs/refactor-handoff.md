# PARE union-BED refactor — cluster handoff

Written 2026-04-23. Use this to resume the refactor on Redwood (CHPC) or locally.

## What changed and why

PARE used to read a per-sample-haplotype `samples.tsv` with paths to assembly
FASTAs and RE BED files, extract slopped sequences per-sample, and align each
sample's REs against its assembly to find paralogs.

It now consumes a single **augmented union peaks BED** — a bgzipped TSV with
one row per `(sample_id × consensus_peak_id)` and pre-extracted ±2 kb flanking
sequence. No assembly FASTAs are needed at run time.

### Upstream: `docs/add_slop_to_union_bed.spec.md`

Spec for the helper that adds the flanking sequence. The output BED has the
same columns as the union BED plus three new ones appended:

- `asm_slop_start` (int)
- `asm_slop_end` (int)
- `slop_seq` (string, uppercase, `peak_length + 4kb` typical, clipped at contig ends)

Existing local file: `chr19-union-peaks-cons.slop2kb.bed.gz` (337 MB, produced
by the user on cluster).

## Input files the pipeline needs

1. `config/config.yaml` — points at:
   - `union_bed:`    the augmented union BED (bgzipped)
   - `samples_tsv:`  small sidecar (header + 90 rows for the chr19 test). See
     the one-liner in `workflow/rules/common.smk` docstring to generate it.
     Base Snakemake parses this, not the BED, so workflow-parse stays fast.
   - optional `ft_bam_table:` enables the pileup rules; columns
     `sample_id`, `ft_bam`.

Already have a valid `config/samples.tsv` for chr19 (generated from the slop
BED, 91 lines including header).

## Architecture

### Alignment unit: haplotype_sample = (Individual_ID, Haplotype)

- 78 unique haplotype_samples across the chr19 union (39 individuals × 2 hap).
- `is_primary_sample == TRUE` picks exactly one `sample_id` per
  (Individual_ID, Haplotype) for alignment. Non-primary sample_ids share the
  same genome, same slop_seq, so they're kept only as graph nodes.

### RE_ID format

Every sequence is named `{sample_id}__{consensus_peak_id}` (double-underscore
separator). `sample_id` can contain single underscores (e.g.
`HG002_PS01015_PacBio_HG002_2`), so splitters use the last occurrence of `__`.
See `workflow/scripts/union_bed.py` (`split_re_id`) and
`workflow/scripts/paf.py` (`split_re_id`).

### Pipeline stages

| Stage                           | Rule(s)                            | Jobs | Notes |
|---------------------------------|------------------------------------|------|-------|
| Parse union BED                 | `extract_union_sequences`          | 1    | polars; emits per-haplotype FASTAs + BEDs + global `re_index.tsv.gz` |
| Per-haplotype all-vs-all align  | `align_within_haplotype_sample`    | 78   | `minimap2 -X -c --eqx -k15 -w10 -n2 -m20 -s0` |
| Filter + classify               | `filter_alignments`                | 78   | identity + peak-coverage thresholds; attaches tags |
| Merge per-haplotype PAFs        | `merge_filtered_pafs`              | 1    |
| Graph + cluster                 | `paf_to_graphml`                   | 1    | Union-Find at `consensus_peak_id` level; expands to every `(sample_id, consensus_peak_id)` node |
| Cross-sample identity matrix    | `align_re_to_re` + `annotate_re_to_re_alignments` | 1+1 | unslopped (peak-only) sequences all-vs-all |
| (optional) pileup               | `ft_pileup` etc.                   | n    | only if `ft_bam_table` is configured |

Total jobs for chr19: 163 (confirmed by dry-run).

### PAF tag schema on filtered records

- `id:Z:` query RE_ID (backward compat)
- `cp:Z:` query consensus_peak_id
- `tp:Z:` target consensus_peak_id
- `qs:Z:` / `ts:Z:` query / target sample_id
- `qh:Z:` / `th:Z:` query / target haplotype_sample (`Individual_ID_Haplotype`)
- `ct:Z:` classification: `self` / `paralog` / `allelic` / `allelic-paralog` /
  `ortholog` / `ortholog-paralog` / `unknown`
- `bi:f:` BLAST identity %
- `qc:f:` query peak coverage %
- `tc:f:` target peak coverage %

(With per-haplotype `-X` alignment, records end up as `paralog`. The classifier
is written to handle future cross-haplotype merges too.)

### Primary node vs collapsed node

User wanted both views. Outputs:

- `results/graphs/re_clusters.tsv`: `re_id / consensus_peak_id / sample_id / cluster_id`
  (sample-level nodes)
- `results/graphs/consensus_peak_clusters.tsv`: `consensus_peak_id / cluster_id`
  (collapsed to consensus peak)

`annotated_res.bed` carries every (sample_id × consensus_peak_id) row with
metadata.

## Running on the cluster

```bash
# adjust paths; make sure union BED + samples TSV exist
pixi run snakemake \
    --profile workflow/profiles/default \
    --config \
        union_bed=/path/to/union.slop2kb.bed.gz \
        samples_tsv=/path/to/samples.tsv \
    -j 64
```

Regenerate `samples.tsv` from the slop BED with the one-liner in
`workflow/rules/common.smk` (header comment block above `_load_samples_tsv`).

## Environments

- **Base env**: `pixi.toml` — just snakemake + ruff + snakefmt. No polars.
  Base env must only read `samples_tsv` (small), not the union BED.
- **Rule envs**: `workflow/envs/python.yml` (polars, pandas) and
  `workflow/envs/bfx.yml` (minimap2, samtools, bedtools, rustybam).

## Known follow-ups / gotchas

1. **Test data** (`test_data/`) is still in the old per-sample-assembly format
   and will not run on the new pipeline. The user indicated they'll test on
   real data on the cluster. If we want a regression test, we'd need a small
   slop-BED fixture.
2. **Analysis Rmd files** (`analysis/*.Rmd`) are edited but unstaged. They
   probably consume the old BED schema; verify after a real run.
3. **Pileup sample→ft_bam mapping** was previously inferred from the old
   `samples.tsv`. Now it's an explicit `ft_bam_table` config key (TSV of
   `sample_id / ft_bam`). Pileup is off by default.
4. **Helper script `add_slop_to_union_bed.py`** — spec is at
   `docs/add_slop_to_union_bed.spec.md`. User ran it on cluster; output file
   lives locally as `chr19-union-peaks-cons.slop2kb.bed.gz`.
5. **polars peak coverage computation** in `extract_union_sequences.py`
   computes the inner unslopped seq as
   `slop_seq[asm_start-asm_slop_start : +(asm_end-asm_start)]`. Verify on
   first run — if any row has `asm_slop_start > asm_start` that's a bug in
   the slop script, not here.
6. **Union BED column 37** (`is_primary_sample`) is present in the local file
   but not listed in `union-peaks-README.md` — README is stale. Data is
   authoritative.

## File inventory for commit

### New files
- `workflow/scripts/union_bed.py` — polars helpers (parsing, filtering, RE_ID builders)
- `workflow/scripts/extract_union_sequences.py` — Snakemake script: BED → FASTAs + BEDs + re_index
- `docs/add_slop_to_union_bed.spec.md` — spec for the upstream slop helper
- `docs/refactor-handoff.md` — this file
- `config/samples.tsv` — generated sidecar (91 lines for chr19)

### Modified files
- `workflow/Snakefile` — new wildcard constraints, new targets
- `workflow/rules/common.smk` — drops assembly helpers; loads samples_tsv via stdlib
- `workflow/rules/extract_sequences.smk` — single extract rule + merged-BED helper
- `workflow/rules/align.smk` — per-haplotype `-X` all-vs-all + filter + merge
- `workflow/rules/graph.smk` — one rule (`paf_to_graphml`)
- `workflow/rules/all_by_all.smk` — cross-sample identity matrix on unslopped seqs
- `workflow/rules/pileup.smk` — column-7 sample_id split (no more contig-name parsing)
- `workflow/scripts/paf.py` — adds `RE_ID_SEP`, `split_re_id` (rest untouched)
- `workflow/scripts/filter_paf.py` — re_index lookup, peak coverage, tag schema above
- `workflow/scripts/paf_to_graph.py` — clean rewrite; clusters at consensus_peak_id, expands via re_index
- `workflow/scripts/annotate_re_to_re.py` — consumes new `re_clusters.tsv` and `re_index.tsv.gz`
- `workflow/scripts/merge_peaks_with_annotations.py` — reads new annotated_res.bed header; `Haplotype == "1"` instead of `hap1`
- `config/config.yaml` — new schema (`union_bed`, `samples_tsv`, `ft_bam_table`)

### Intentionally NOT committed (local-only)
- `chr19-union-peaks-cons.bed.gz` / `chr19-union-peaks-cons.slop2kb.bed.gz` — data, large
- `test_data/`, `data/`, `figures/`, `tmp.pdf`, `analysis/tmp.pdf`, `paralog.paf`, `paralog_graph.graphml`
- `claude.md` — user's scratch

### Analysis Rmd files (modified but deferred)
- `analysis/00_data_loader.Rmd` / `01_paralog_variance.Rmd` / `02_paralog_vs_allelic.Rmd` / `03_fire_diff_by_identity.Rmd` / `paralog_accessibility.Rmd`

These predate the rewrite; defer committing until we re-run the pipeline on
real data and confirm the analyses still work against the new output schema.
