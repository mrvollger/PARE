#!/usr/bin/env python3
"""
Filter ft pileup output to find peak positions within each RE.

For each continuous run of the same 'name' (RE), find the position with:
1. Maximum score
2. If tied, maximum fire_coverage
3. If still tied, most central position
"""

import sys
import polars as pl


def filter_peaks(df: pl.DataFrame) -> pl.DataFrame:
    """Filter to best peak position per continuous RE group."""
    # Create group IDs for continuous runs of the same name
    df = df.with_columns(
        (pl.col("name") != pl.col("name").shift(1)).cum_sum().alias("group_id")
    )

    # Calculate group boundaries for central position calculation
    group_stats = df.group_by("group_id").agg(
        pl.col("start").min().alias("group_start"),
        pl.col("end").max().alias("group_end"),
    )

    df = df.join(group_stats, on="group_id")

    # Calculate center of each row and distance to group center
    df = df.with_columns(
        ((pl.col("start") + pl.col("end")) / 2).alias("row_center"),
        ((pl.col("group_start") + pl.col("group_end")) / 2).alias("group_center"),
    )
    df = df.with_columns(
        (pl.col("row_center") - pl.col("group_center")).abs().alias("dist_to_center")
    )

    # Sort by group_id, then by score (desc), fire_coverage (desc), dist_to_center (asc)
    df = df.sort(
        ["group_id", "score", "fire_coverage", "dist_to_center"],
        descending=[False, True, True, False],
    )

    # Take first row per group (best by our criteria)
    df = df.group_by("group_id", maintain_order=True).first()

    # Drop helper columns
    df = df.drop(
        ["group_id", "group_start", "group_end", "row_center", "group_center", "dist_to_center"]
    )

    return df


def snakemake_main():
    input_file = snakemake.input.pileup
    output_file = snakemake.output.peaks
    log_file = snakemake.log[0] if snakemake.log else None

    if log_file:
        log = open(log_file, "w")
        old_stderr = sys.stderr
        sys.stderr = log

    try:
        print(f"Reading pileup file: {input_file}", file=sys.stderr)

        # Read the gzipped pileup file
        df = pl.read_csv(
            input_file,
            separator="\t",
            comment_prefix=None,
        )

        # Rename first column if it has # prefix
        if df.columns[0].startswith("#"):
            df = df.rename({df.columns[0]: df.columns[0].lstrip("#")})

        print(f"Input rows: {len(df)}", file=sys.stderr)

        # Filter to peaks
        result = filter_peaks(df)

        print(f"Output rows (peaks): {len(result)}", file=sys.stderr)

        # Write uncompressed output (will be piped to bgzip)
        result.write_csv(output_file, separator="\t")

        print(f"Wrote output to: {output_file}", file=sys.stderr)

    finally:
        if log_file:
            sys.stderr = old_stderr
            log.close()


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_main()
    else:
        # CLI mode for testing
        df = pl.read_csv(
            sys.stdin,
            separator="\t",
        )
        if df.columns[0].startswith("#"):
            df = df.rename({df.columns[0]: df.columns[0].lstrip("#")})

        result = filter_peaks(df)
        result.write_csv(sys.stdout, separator="\t")
