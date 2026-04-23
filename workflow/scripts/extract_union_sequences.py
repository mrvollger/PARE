#!/usr/bin/env python3
"""
Parse the augmented union peaks BED and emit:
  1. one slopped FASTA per haplotype_sample (alignment queries/targets)
  2. one unslopped FASTA per haplotype_sample (for all-vs-all RE identity)
  3. one unslopped BED per haplotype_sample (asm_chr, asm_start, asm_end, re_id)
  4. one slopped BED per haplotype_sample (asm_chr, asm_slop_start, asm_slop_end, re_id)
  5. a single re_index.tsv.gz covering every (sample_id × consensus_peak_id) row
     (primary and non-primary) with metadata needed to expand paralog edges
     from the haplotype_sample level back to per-sample nodes.

Only is_primary_sample == TRUE rows end up in the FASTAs/BEDs (alignment input).
The re_index includes every row.

Invoked as a Snakemake script with:
    input.bed              path to union.bed.gz
    output.fastas          list of slopped FASTAs, one per haplotype_sample
    output.fastas_unslop   list of unslopped FASTAs, one per haplotype_sample
    output.beds            list of unslopped BEDs, one per haplotype_sample
    output.slop_beds       list of slopped BEDs, one per haplotype_sample
    output.re_index        gzipped TSV
    params.haplotype_samples   ordered list matching the FASTA/BED lists
    params.min_len / params.max_len / params.drop_sd  (optional filters on alignment input)
    log[0]                 log file path
"""

from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path

# snakemake script: add this dir to path so union_bed is importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

import polars as pl  # noqa: E402
from union_bed import (  # noqa: E402
    RE_ID_SEP,
    build_re_index,
    filter_for_alignment,
    read_union_bed,
    with_re_id,
)


def _index_map(paths: list[str], haplotype_samples: list[str]) -> dict[str, str]:
    """Map haplotype_sample -> file path, matching list positions from the rule."""
    if len(paths) != len(haplotype_samples):
        raise ValueError(
            f"path count ({len(paths)}) != haplotype_sample count ({len(haplotype_samples)})"
        )
    return dict(zip(haplotype_samples, paths))


def main() -> None:
    sm = snakemake  # noqa: F821 (injected)

    bed_path: str = sm.input.bed
    hs_list: list[str] = list(sm.params.haplotype_samples)

    fasta_map = _index_map(list(sm.output.fastas), hs_list)
    fasta_unslop_map = _index_map(list(sm.output.fastas_unslop), hs_list)
    bed_map = _index_map(list(sm.output.beds), hs_list)
    slop_bed_map = _index_map(list(sm.output.slop_beds), hs_list)
    re_index_path: str = sm.output.re_index

    min_len = getattr(sm.params, "min_len", None)
    max_len = getattr(sm.params, "max_len", None)
    drop_sd = bool(getattr(sm.params, "drop_sd", False))

    log = open(sm.log[0], "w")

    def say(msg: str) -> None:
        print(msg, file=log, flush=True)

    say(f"reading union BED: {bed_path}")
    df = read_union_bed(bed_path)
    say(f"  rows: {df.height}")
    say(f"  unique sample_id: {df.get_column('sample_id').n_unique()}")
    say(f"  unique consensus_peak_id: {df.get_column('consensus_peak_id').n_unique()}")

    # Emit the global RE index first (covers primary + non-primary rows).
    say("building RE index (all rows)")
    re_idx = build_re_index(df)
    # Write gzipped TSV
    with gzip.open(re_index_path, "wt") as fh:
        fh.write("\t".join(re_idx.columns) + "\n")
        for row in re_idx.iter_rows():
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")
    say(f"  wrote {re_idx.height} rows to {re_index_path}")

    # Filter to alignment input: primary rows with slop_seq.
    aln = filter_for_alignment(df, drop_sd=drop_sd, min_len=min_len, max_len=max_len)
    aln = with_re_id(aln)
    say(f"alignment input rows after filters: {aln.height}")

    # Validate haplotype_sample coverage
    present_hs = set(aln.get_column("haplotype_sample").unique().to_list())
    missing_hs = [h for h in hs_list if h not in present_hs]
    if missing_hs:
        say(f"WARNING: haplotype_samples with no primary rows: {missing_hs}")

    # Split and write per-haplotype_sample files.
    for hs in hs_list:
        sub = aln.filter(pl.col("haplotype_sample") == hs)
        n = sub.height

        # slopped FASTA
        with open(fasta_map[hs], "w") as fh:
            for re_id, seq in sub.select(["re_id", "slop_seq"]).iter_rows():
                fh.write(f">{re_id}\n{seq}\n")

        # unslopped FASTA (inner peak sequence = slop_seq[asm_start-asm_slop_start : that + peak_len])
        # We compute it on the fly so we don't depend on a separate `seq` column.
        inner = sub.with_columns(
            (pl.col("asm_start") - pl.col("asm_slop_start")).alias("_inner_off"),
            (pl.col("asm_end") - pl.col("asm_start")).alias("_inner_len"),
        )
        with open(fasta_unslop_map[hs], "w") as fh:
            for re_id, seq, off, length in inner.select(
                ["re_id", "slop_seq", "_inner_off", "_inner_len"]
            ).iter_rows():
                if seq is None or off is None or length is None:
                    continue
                fh.write(f">{re_id}\n{seq[off : off + length]}\n")

        # unslopped BED (asm_chr, asm_start, asm_end, re_id)
        with open(bed_map[hs], "w") as fh:
            for chrom, start, end, re_id in sub.select(
                ["asm_chr", "asm_start", "asm_end", "re_id"]
            ).iter_rows():
                fh.write(f"{chrom}\t{start}\t{end}\t{re_id}\n")

        # slopped BED (asm_chr, asm_slop_start, asm_slop_end, re_id)
        with open(slop_bed_map[hs], "w") as fh:
            for chrom, start, end, re_id in sub.select(
                ["asm_chr", "asm_slop_start", "asm_slop_end", "re_id"]
            ).iter_rows():
                fh.write(f"{chrom}\t{start}\t{end}\t{re_id}\n")

        say(f"  {hs}: {n} primary rows → FASTAs + BEDs")

    log.close()


if __name__ == "__main__":
    main()
