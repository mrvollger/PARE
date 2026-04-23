#!/usr/bin/env python3
"""
Parser helpers for the augmented union peaks BED.

The input is a bgzipped TSV (header starts with '#') produced by combining the
union peak table with `slop_seq`, `asm_slop_start`, `asm_slop_end` (see
docs/add_slop_to_union_bed.spec.md).

Key columns used downstream:
  - consensus_peak_id   : cross-sample consensus peak ID
  - sample_id           : full composite sample id (e.g. HG002_PS01015_PacBio_HG002_2)
  - Individual_ID       : e.g. HG002
  - Haplotype           : 1 or 2
  - asm_chr / asm_start / asm_end              : coords in that sample's assembly
  - asm_slop_start / asm_slop_end              : slopped coords
  - slop_seq                                    : ±slop flanking sequence
  - is_peak             : peak called in this sample/hap (TRUE/FALSE)
  - is_primary_sample   : used to dedupe down to one sample_id per (Individual_ID, Haplotype)
  - is_SD_cons / is_SD_asm / is_SD             : SD overlap flags

The RE_ID used everywhere downstream in PARE is:
    {sample_id}__{consensus_peak_id}         (double underscore separator)

Alignment is done once per haplotype_sample = (Individual_ID, Haplotype), using
only rows where is_primary_sample == TRUE. Non-primary rows share the same
assembly sequence and are only needed as graph nodes, not alignment queries.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Iterable, Iterator, Optional

import polars as pl


RE_ID_SEP = "__"

# Columns we actually need — read only these for speed
REQUIRED_COLS = [
    "consensus_peak_id",
    "sample_id",
    "Individual_ID",
    "Haplotype",
    "asm_chr",
    "asm_start",
    "asm_end",
    "asm_slop_start",
    "asm_slop_end",
    "slop_seq",
    "is_peak",
    "is_primary_sample",
    "is_SD_cons",
    "is_SD_asm",
    "is_SD",
]


def make_re_id(sample_id: str, consensus_peak_id: str) -> str:
    return f"{sample_id}{RE_ID_SEP}{consensus_peak_id}"


def split_re_id(re_id: str) -> tuple[str, str]:
    """Split RE_ID back into (sample_id, consensus_peak_id). Raises on malformed."""
    idx = re_id.rfind(RE_ID_SEP)
    if idx < 0:
        raise ValueError(f"RE_ID missing '{RE_ID_SEP}' separator: {re_id!r}")
    return re_id[:idx], re_id[idx + len(RE_ID_SEP):]


def haplotype_sample_key(individual_id: str, haplotype: str | int) -> str:
    """Unique key for a haplotype sample: e.g. HG002_1."""
    return f"{individual_id}_{haplotype}"


def read_union_bed(path: str, columns: Optional[Iterable[str]] = None) -> pl.DataFrame:
    """
    Read the augmented union BED (bgzipped or plain) into a polars DataFrame.

    The header begins with '#' (e.g. '#chrom'); polars strips only literal
    matches, so we rename the first column after read.

    By default reads only REQUIRED_COLS; pass columns=... to override.
    """
    use_cols = list(columns) if columns is not None else REQUIRED_COLS

    # Peek header so we can request only the columns we need by name.
    # polars handles .gz via compression='gzip' but bgzip files are gzip-compatible.
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=True,
        comment_prefix=None,
        null_values=["NA", ""],
        columns=use_cols,
        infer_schema_length=10000,
    )

    # Validate required columns are present
    missing = [c for c in use_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"union BED at {path} missing required columns: {missing}"
        )

    # Coerce types we rely on (polars may read these as ints/floats/bools depending on content)
    df = df.with_columns(
        pl.col("Haplotype").cast(pl.Utf8),
        pl.col("asm_start").cast(pl.Int64, strict=False),
        pl.col("asm_end").cast(pl.Int64, strict=False),
        pl.col("asm_slop_start").cast(pl.Int64, strict=False),
        pl.col("asm_slop_end").cast(pl.Int64, strict=False),
    )

    return df


def with_re_id(df: pl.DataFrame) -> pl.DataFrame:
    """Add re_id and haplotype_sample columns."""
    return df.with_columns(
        (pl.col("sample_id") + pl.lit(RE_ID_SEP) + pl.col("consensus_peak_id")).alias("re_id"),
        (pl.col("Individual_ID") + pl.lit("_") + pl.col("Haplotype")).alias("haplotype_sample"),
    )


def filter_for_alignment(
    df: pl.DataFrame,
    *,
    drop_sd: bool = False,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
) -> pl.DataFrame:
    """
    Keep only rows that should be aligned:
      - is_primary_sample == TRUE (one sample_id per haplotype_sample)
      - has slop_seq and asm coords
      - optional length and SD filters
    """
    out = df.filter(
        (pl.col("is_primary_sample") == True)  # noqa: E712
        & pl.col("slop_seq").is_not_null()
        & (pl.col("slop_seq").str.len_chars() > 0)
        & pl.col("asm_chr").is_not_null()
        & pl.col("asm_start").is_not_null()
        & pl.col("asm_end").is_not_null()
    )

    if drop_sd:
        out = out.filter(pl.col("is_SD") != True)  # noqa: E712

    if min_len is not None:
        out = out.filter((pl.col("asm_end") - pl.col("asm_start")) >= min_len)
    if max_len is not None:
        out = out.filter((pl.col("asm_end") - pl.col("asm_start")) <= max_len)

    return out


def haplotype_samples(df: pl.DataFrame) -> list[str]:
    """Ordered list of unique haplotype_sample keys in the table (primary-sample rows only)."""
    prim = df.filter(pl.col("is_primary_sample") == True)  # noqa: E712
    return (
        prim.select(
            (pl.col("Individual_ID") + pl.lit("_") + pl.col("Haplotype")).alias("hs")
        )
        .unique()
        .sort("hs")
        .get_column("hs")
        .to_list()
    )


def sample_ids(df: pl.DataFrame) -> list[str]:
    """All unique sample_ids (primary and non-primary)."""
    return (
        df.select("sample_id")
        .unique()
        .sort("sample_id")
        .get_column("sample_id")
        .to_list()
    )


def write_fasta(records: Iterable[tuple[str, str]], path: str) -> int:
    """Write (name, seq) pairs as FASTA. Returns record count."""
    n = 0
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
            n += 1
    return n


@dataclass
class REIndexRow:
    """One row of the RE index: one per (sample_id × consensus_peak_id)."""
    re_id: str
    sample_id: str
    consensus_peak_id: str
    haplotype_sample: str
    individual_id: str
    haplotype: str
    asm_chr: Optional[str]
    asm_start: Optional[int]
    asm_end: Optional[int]
    asm_slop_start: Optional[int]
    asm_slop_end: Optional[int]
    is_peak: Optional[bool]
    is_primary_sample: Optional[bool]
    is_sd: Optional[bool]


RE_INDEX_COLS = [
    "re_id",
    "sample_id",
    "consensus_peak_id",
    "haplotype_sample",
    "Individual_ID",
    "Haplotype",
    "asm_chr",
    "asm_start",
    "asm_end",
    "asm_slop_start",
    "asm_slop_end",
    "is_peak",
    "is_primary_sample",
    "is_SD",
]


def build_re_index(df: pl.DataFrame) -> pl.DataFrame:
    """Build the per-(sample_id × consensus_peak_id) index with RE_ID + metadata."""
    indexed = with_re_id(df)
    return indexed.select(
        [
            "re_id",
            "sample_id",
            "consensus_peak_id",
            "haplotype_sample",
            "Individual_ID",
            "Haplotype",
            "asm_chr",
            "asm_start",
            "asm_end",
            "asm_slop_start",
            "asm_slop_end",
            "is_peak",
            "is_primary_sample",
            "is_SD",
        ]
    )


if __name__ == "__main__":
    # Quick smoke test: `python union_bed.py <bed.gz>` prints summary stats
    if len(sys.argv) != 2:
        print("usage: union_bed.py <union.bed.gz>", file=sys.stderr)
        sys.exit(2)
    df = read_union_bed(sys.argv[1])
    print(f"rows: {df.height}")
    print(f"unique sample_id: {df.get_column('sample_id').n_unique()}")
    print(f"unique consensus_peak_id: {df.get_column('consensus_peak_id').n_unique()}")
    print(f"haplotype_samples: {len(haplotype_samples(df))}")
    print(f"is_primary_sample=TRUE rows: {df.filter(pl.col('is_primary_sample')==True).height}")
