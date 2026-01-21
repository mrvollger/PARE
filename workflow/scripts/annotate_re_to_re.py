#!/usr/bin/env python3
"""
Annotate RE-to-RE alignments with cluster IDs and calculate pairwise identity.

This script:
1. Parses RE IDs from sequence names in the PAF file
2. Joins with cluster assignments from annotated_res.bed
3. Removes self-alignments (same RE to itself)
4. Calculates sequence identity from alignment stats (both BLAST and gap-compressed)
5. Outputs annotated PAF and a simplified TSV for analysis
"""

import gzip
import sys
from pathlib import Path

import pandas as pd

# Import PAF module from same directory
sys.path.insert(0, str(Path(__file__).parent))
from paf import PAFRecord


def parse_re_sequence_name(seq_name: str) -> str:
    """
    Parse RE ID from sequence name.

    Sequence names from bedtools getfasta -name have format:
    chrom_start_end::chrom:slopped_start-slopped_end

    We want to extract the original RE ID: chrom_start_end
    """
    # Split on :: to get the name part (before ::)
    if "::" in seq_name:
        return seq_name.split("::")[0]
    return seq_name


def main():
    # Snakemake inputs
    paf_file = snakemake.input.paf
    annotated_res_file = snakemake.input.annotated_res

    # Outputs
    out_paf = snakemake.output.paf
    out_tsv = snakemake.output.tsv

    # Set up logging
    log_file = snakemake.log[0]
    with open(log_file, "w") as log:
        log.write(f"Processing RE-to-RE alignments\n")
        log.write(f"Input PAF: {paf_file}\n")
        log.write(f"Annotated REs: {annotated_res_file}\n")

        # Load cluster assignments from annotated_res.bed
        # Format: #chrom, start, end, re_id, cluster_id, sample, haplotype
        log.write("Loading cluster assignments...\n")
        re_clusters = {}
        re_info = {}

        with open(annotated_res_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 7:
                    continue
                re_id = fields[3]
                cluster_id = fields[4]
                sample = fields[5]
                haplotype = fields[6]
                re_clusters[re_id] = cluster_id
                re_info[re_id] = {
                    "cluster_id": cluster_id,
                    "sample": sample,
                    "haplotype": haplotype,
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                }

        log.write(f"Loaded {len(re_clusters)} RE cluster assignments\n")

        # Process PAF alignments
        log.write("Processing PAF alignments...\n")

        alignments = []
        stats = {
            "total": 0,
            "self_removed": 0,
            "query_not_found": 0,
            "target_not_found": 0,
            "kept": 0,
        }

        with open(paf_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue

                stats["total"] += 1
                fields = line.strip().split("\t")
                if len(fields) < 12:
                    continue

                # Parse standard PAF fields
                query_name = fields[0]
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target_name = fields[5]
                target_length = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                num_matches = int(fields[9])
                alignment_length = int(fields[10])
                mapq = int(fields[11])

                # Parse RE IDs from sequence names
                query_re_id = parse_re_sequence_name(query_name)
                target_re_id = parse_re_sequence_name(target_name)

                # Skip self-alignments
                if query_re_id == target_re_id:
                    stats["self_removed"] += 1
                    continue

                # Look up cluster IDs
                query_cluster = re_clusters.get(query_re_id)
                target_cluster = re_clusters.get(target_re_id)

                if query_cluster is None:
                    stats["query_not_found"] += 1
                    continue
                if target_cluster is None:
                    stats["target_not_found"] += 1
                    continue

                stats["kept"] += 1

                # Parse line into PAFRecord to use identity methods
                paf_record = PAFRecord.from_line(line)

                # Calculate both identity metrics using PAF module
                blast_identity = paf_record.blast_identity() / 100.0  # Convert to fraction
                gap_compressed_identity = paf_record.gap_compressed_identity()
                if gap_compressed_identity is not None:
                    gap_compressed_identity = gap_compressed_identity / 100.0  # Convert to fraction

                # Get additional info
                query_info = re_info.get(query_re_id, {})
                target_info = re_info.get(target_re_id, {})

                # Determine if same cluster
                same_cluster = query_cluster == target_cluster

                alignments.append(
                    {
                        "query_re_id": query_re_id,
                        "target_re_id": target_re_id,
                        "query_cluster": query_cluster,
                        "target_cluster": target_cluster,
                        "same_cluster": same_cluster,
                        "query_sample": query_info.get("sample", "unknown"),
                        "target_sample": target_info.get("sample", "unknown"),
                        "query_haplotype": query_info.get("haplotype", "unknown"),
                        "target_haplotype": target_info.get("haplotype", "unknown"),
                        "blast_identity": blast_identity,
                        "gap_compressed_identity": gap_compressed_identity,
                        "num_matches": num_matches,
                        "alignment_length": alignment_length,
                        "query_length": query_length,
                        "target_length": target_length,
                        "query_start": query_start,
                        "query_end": query_end,
                        "target_start": target_start,
                        "target_end": target_end,
                        "strand": strand,
                        "mapq": mapq,
                        # Store original line for PAF output
                        "_original_line": line.strip(),
                    }
                )

        log.write(f"\nAlignment statistics:\n")
        log.write(f"  Total alignments: {stats['total']}\n")
        log.write(f"  Self-alignments removed: {stats['self_removed']}\n")
        log.write(f"  Query RE not found: {stats['query_not_found']}\n")
        log.write(f"  Target RE not found: {stats['target_not_found']}\n")
        log.write(f"  Alignments kept: {stats['kept']}\n")

        # Convert to DataFrame
        df = pd.DataFrame(alignments)

        if len(df) > 0:
            # Count same-cluster alignments
            same_cluster_count = df["same_cluster"].sum()
            log.write(f"\nSame-cluster alignments: {same_cluster_count}\n")
            log.write(
                f"Cross-cluster alignments: {len(df) - same_cluster_count}\n"
            )

            # Deduplicate: keep only one alignment per RE pair (the one with highest identity)
            # Create a canonical pair key (sorted alphabetically to treat A-B same as B-A)
            log.write("\nDeduplicating RE pairs...\n")
            log.write(f"  Alignments before dedup: {len(df)}\n")

            def make_pair_key(row):
                """Create canonical pair key (sorted so A-B == B-A)."""
                pair = sorted([row["query_re_id"], row["target_re_id"]])
                return f"{pair[0]}|{pair[1]}"

            df["pair_key"] = df.apply(make_pair_key, axis=1)

            # Sort by (num_matches * blast_identity) to prefer longer, high-quality alignments
            # This balances alignment length and identity
            df["_sort_score"] = df["num_matches"] * df["blast_identity"]
            df = df.sort_values("_sort_score", ascending=False)
            df = df.drop_duplicates(subset=["pair_key"], keep="first")
            df = df.drop(columns=["pair_key", "_sort_score"])

            log.write(f"  Alignments after dedup: {len(df)}\n")

            # Verify no duplicate pairs remain
            pair_counts = df.apply(
                lambda r: tuple(sorted([r["query_re_id"], r["target_re_id"]])),
                axis=1,
            ).value_counts()
            duplicates = pair_counts[pair_counts > 1]
            if len(duplicates) > 0:
                log.write(f"  WARNING: {len(duplicates)} duplicate pairs remain!\n")
                for pair, count in duplicates.head(5).items():
                    log.write(f"    {pair}: {count}\n")
            else:
                log.write("  Verified: no duplicate RE pairs in output.\n")

            # Write annotated PAF (gzipped)
            log.write(f"\nWriting annotated PAF to {out_paf}...\n")
            with gzip.open(out_paf, "wt") as f:
                for _, row in df.iterrows():
                    # Add cluster tags to original PAF line
                    paf_line = row["_original_line"]
                    paf_line += f"\tqc:Z:{row['query_cluster']}"
                    paf_line += f"\ttc:Z:{row['target_cluster']}"
                    paf_line += f"\tsc:i:{1 if row['same_cluster'] else 0}"
                    paf_line += f"\tbi:f:{row['blast_identity']:.4f}"
                    if row['gap_compressed_identity'] is not None:
                        paf_line += f"\tgi:f:{row['gap_compressed_identity']:.4f}"
                    f.write(paf_line + "\n")

            # Write simplified TSV for R analysis (gzipped)
            log.write(f"Writing TSV to {out_tsv}...\n")
            tsv_cols = [
                "query_re_id",
                "target_re_id",
                "query_cluster",
                "target_cluster",
                "same_cluster",
                "query_sample",
                "target_sample",
                "query_haplotype",
                "target_haplotype",
                "blast_identity",
                "gap_compressed_identity",
                "num_matches",
                "alignment_length",
                "query_length",
                "target_length",
                "strand",
            ]
            df[tsv_cols].to_csv(out_tsv, sep="\t", index=False, compression="gzip")

            log.write(f"Done. Wrote {len(df)} unique RE pair alignments.\n")
        else:
            log.write("\nNo alignments to write.\n")
            # Write empty files with header
            with gzip.open(out_paf, "wt") as f:
                pass
            with gzip.open(out_tsv, "wt") as f:
                f.write("\t".join([
                    "query_re_id",
                    "target_re_id",
                    "query_cluster",
                    "target_cluster",
                    "same_cluster",
                    "query_sample",
                    "target_sample",
                    "query_haplotype",
                    "target_haplotype",
                    "blast_identity",
                    "gap_compressed_identity",
                    "num_matches",
                    "alignment_length",
                    "query_length",
                    "target_length",
                    "strand",
                ]) + "\n")


if __name__ == "__main__":
    main()
