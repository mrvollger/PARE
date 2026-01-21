#!/usr/bin/env python3
"""
Add peak data to annotated REs.

Joins annotated_res.bed with peaks (keyed by re_id) to add chromatin
accessibility information to each RE.
Verifies 1:1 mapping between re_id and peaks.
"""

import sys
import subprocess
import polars as pl


def snakemake_main():
    peaks_file = snakemake.input.peaks
    annotated_res_file = snakemake.input.annotated_res
    output_bed = snakemake.output.bed
    output_tbi = snakemake.output.tbi
    log_file = snakemake.log[0] if snakemake.log else None

    if log_file:
        log = open(log_file, "w")
        old_stderr = sys.stderr
        sys.stderr = log

    try:
        print(f"Loading annotated REs: {annotated_res_file}", file=sys.stderr)
        annotations = pl.read_csv(
            annotated_res_file,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "re_id", "cluster_id", "sample", "haplotype"],
        )
        print(f"Annotations loaded: {len(annotations)} rows", file=sys.stderr)

        print(f"Loading peaks: {peaks_file}", file=sys.stderr)
        peaks = pl.read_csv(
            peaks_file,
            separator="\t",
            infer_schema_length=10000,
        )

        # Rename first column if it has # prefix
        if peaks.columns[0].startswith("#"):
            peaks = peaks.rename({peaks.columns[0]: peaks.columns[0].lstrip("#")})

        print(f"Peaks loaded: {len(peaks)} rows", file=sys.stderr)

        # Check for duplicate re_ids in peaks
        re_id_counts = peaks.group_by("name").agg(pl.len().alias("count"))
        duplicates = re_id_counts.filter(pl.col("count") > 1)
        if len(duplicates) > 0:
            print(
                f"WARNING: Found {len(duplicates)} duplicate re_ids in peaks:",
                file=sys.stderr,
            )
            print(duplicates.head(10), file=sys.stderr)

        # Rename peaks columns to avoid collision, keep 'name' as join key
        # Drop chrom, start, end from peaks since we use the annotation coords
        peak_cols = [c for c in peaks.columns if c not in ["chrom", "start", "end"]]
        peaks_to_join = peaks.select(peak_cols)

        # Join annotations with peaks on re_id
        merged = annotations.join(
            peaks_to_join,
            left_on="re_id",
            right_on="name",
            how="left",
        )

        # Check join results
        missing_peaks = merged.filter(pl.col("score").is_null())
        if len(missing_peaks) > 0:
            print(
                f"WARNING: {len(missing_peaks)} REs have no peak data",
                file=sys.stderr,
            )

        has_peaks = merged.filter(pl.col("score").is_not_null())
        print(f"REs with peak data: {len(has_peaks)}", file=sys.stderr)
        print(f"Merged result: {len(merged)} rows", file=sys.stderr)

        # Rename _H1/_H2 columns to _ref/_alt based on haplotype
        # For hap1 REs: H1 -> ref, H2 -> alt
        # For hap2 REs: H2 -> ref, H1 -> alt
        h1_cols = [c for c in merged.columns if c.endswith("_H1")]
        h2_cols = [c for c in merged.columns if c.endswith("_H2")]

        # Create renamed columns based on haplotype
        rename_exprs = []
        for col in h1_cols:
            base = col[:-3]  # Remove _H1
            ref_name = f"{base}_ref"
            alt_name = f"{base}_alt"
            # When hap1: H1 is ref, H2 is alt
            # When hap2: H2 is ref, H1 is alt
            rename_exprs.append(
                pl.when(pl.col("haplotype") == "hap1")
                .then(pl.col(col))
                .otherwise(pl.col(col.replace("_H1", "_H2")))
                .alias(ref_name)
            )
            rename_exprs.append(
                pl.when(pl.col("haplotype") == "hap1")
                .then(pl.col(col.replace("_H1", "_H2")))
                .otherwise(pl.col(col))
                .alias(alt_name)
            )

        # Select non-haplotype columns plus the new ref/alt columns
        non_hap_cols = [c for c in merged.columns if not c.endswith("_H1") and not c.endswith("_H2")]
        merged = merged.select(non_hap_cols + rename_exprs)

        print(f"Renamed _H1/_H2 columns to _ref/_alt based on haplotype", file=sys.stderr)

        # Sort by chrom and position for BED format
        merged = merged.sort(["chrom", "start"])

        # Rename first column to start with # for tabix compatibility
        merged = merged.rename({"chrom": "#chrom"})

        # Write to temp file, then bgzip and tabix
        temp_file = output_bed.replace(".gz", "")
        merged.write_csv(temp_file, separator="\t", include_header=True)

        # bgzip
        subprocess.run(["bgzip", "-f", temp_file], check=True)

        # tabix
        subprocess.run(["tabix", "-p", "bed", output_bed], check=True)

        print(f"Wrote output to: {output_bed}", file=sys.stderr)

    finally:
        if log_file:
            sys.stderr = old_stderr
            log.close()


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_main()
    else:
        print("This script is designed to be run via Snakemake", file=sys.stderr)
        sys.exit(1)
