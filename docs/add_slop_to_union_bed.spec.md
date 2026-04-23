# `add_slop_to_union_bed.py` — Spec

## Purpose
Augment a union peak BED (like `chr19-union-peaks-cons.bed.gz`) with ±N bp of flanking assembly sequence, so downstream tools (PARE) don't need to touch the assemblies.

## Inputs

1. **Union BED** (`--bed`) — bgzip'd TSV with header. Must contain columns: `Individual_ID`, `asm_chr`, `asm_start`, `asm_end`. Preserve all other columns untouched.
2. **Assembly list** (`--fastas`) — a text file with one FASTA path per line, e.g.:
   ```
   /mmfs1/gscratch/stergachislab/HPRC/assemblies/diploid/NA20503.diploid.fa.gz
   /mmfs1/gscratch/stergachislab/HPRC/assemblies/diploid/NA20752.diploid.fa.gz
   /mmfs1/gscratch/stergachislab/HPRC/assemblies/diploid/NA20762.diploid.fa.gz
   ...
   ```
   Each FASTA must be bgzipped and have both `.fai` and `.gzi` indexes alongside (error out with a clear message if missing).
3. **Slop size** (`--slop`, default `2000`).
4. **Output path** (`--out`) — write bgzipped TSV.

## Individual_ID → assembly mapping

- Derive `Individual_ID` from each FASTA as the filename basename before the first `.` (e.g. `NA20503.diploid.fa.gz` → `NA20503`).
- Build a dict `{Individual_ID: fasta_path}`. If two FASTAs map to the same ID, error out.
- If a BED row's `Individual_ID` isn't in the dict, **error out by default**, but accept a `--skip-missing` flag that instead drops those rows and logs a count.
- HG002-style `Individual_ID=HG002` needs its own diploid FASTA in the list (user will add); script should not hardcode anything.

## Per-row logic

For each row:

1. Skip rows where `asm_chr`, `asm_start`, or `asm_end` is `NA` / empty — write them through unchanged with empty new columns, and log the count.
2. Load the assembly's chrom length from its `.fai` (cache per-assembly).
3. Compute:
   - `asm_slop_start = max(0, asm_start - slop)`
   - `asm_slop_end   = min(chrom_len, asm_end + slop)`
4. Fetch `slop_seq` from the FASTA for `asm_chr:asm_slop_start-asm_slop_end` (0-based half-open). Uppercase the result.
5. Sanity-check: `len(slop_seq) == asm_slop_end - asm_slop_start`. Error (not warn) on mismatch — indicates index/fasta skew.

## Output

Same columns as input, in the same order, with three new columns appended:

- `asm_slop_start` (int)
- `asm_slop_end` (int)
- `slop_seq` (string, uppercase)

Keep the leading `#` on the header. Write with bgzip.

## Implementation notes

- Use **pysam** `FastaFile` — handles `.fa.gz` + `.fai` + `.gzi` natively and is fast. (`conda install -c bioconda pysam`.)
- **Group rows by `Individual_ID`** before fetching, so you open each FASTA once. Many rows share the same assembly (e.g. all `HG002_*` rows point to the same HG002 FASTA), and `FastaFile` construction is slow.
- Stream row-by-row on write — don't hold the output in memory. Reading the input fully into a DataFrame is fine; it's ~1.4M rows.
- Preserve column dtypes as strings on write to avoid `float` coercion of NA-containing int columns.
- `#!/usr/bin/env python3`. Any wrapper shell: `set -euo pipefail`.

## CLI

```
add_slop_to_union_bed.py \
    --bed chr19-union-peaks-cons.bed.gz \
    --fastas assemblies.txt \
    --slop 2000 \
    --out chr19-union-peaks-cons.slop2kb.bed.gz \
    [--skip-missing] \
    [--threads 4]
```

## Sanity tests to run after

```bash
# Row count unchanged
zcat in.bed.gz  | wc -l
zcat out.bed.gz | wc -l

# New columns present at the end of the header
zcat out.bed.gz | head -1 | tr '\t' '\n' | tail -5

# slop_seq length matches slopped coords
zcat out.bed.gz \
  | awk -F'\t' 'NR>1 && $NF != "" {
      n=length($NF); w=$(NF-1)-$(NF-2);
      if (n!=w) { print "MISMATCH", NR, n, w; exit 1 }
    } END { print "OK" }'

# Inner peak sequence matches the original `seq` column:
# slop_seq[ (asm_start - asm_slop_start) : (asm_start - asm_slop_start) + (asm_end - asm_start) ] == seq
```

## Handoff

Once `chr19-union-peaks-cons.slop2kb.bed.gz` exists, hand it back to PARE. PARE will be rewritten to:

- Consume this file directly (no `samples.tsv`, no assembly rules).
- Use `slop_seq` as the alignment query.
- Use `asm_slop_start` / `asm_slop_end` + `asm_chr` to map PAF hits back to real assembly coords.
- Classify alignments using `sample_id × consensus_peak_id` as the primary node, with `consensus_peak_id` preserved on every record for downstream collapsing.
