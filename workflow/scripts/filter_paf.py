#!/usr/bin/env python3
"""
Filter PAF alignments based on identity and coverage thresholds.
Classify alignments as self, paralogous, or allelic.
Output PAF format with filtered alignments and classification tags.
"""

import sys
import logging
from paf import PAFReader


def main():
    paf_file = snakemake.input.paf
    out_paf_file = snakemake.output.paf
    min_identity = snakemake.params.min_identity
    min_coverage = snakemake.params.min_coverage

    passed = 0
    total = 0
    classification_counts = {
        "self": 0,
        "allelic": 0,
        "paralog": 0,
        "ortholog": 0,
        "unknown": 0,
    }

    # Read PAF file using PAFReader
    reader = PAFReader(paf_file)
    filtered_records = []

    for record in reader:
        total += 1

        # Calculate metrics using PAFRecord methods
        identity = record.calculate_identity()
        coverage = record.calculate_re_coverage()
        overlap = record.calculate_self_overlap()

        # Filter based on thresholds
        if identity >= min_identity and coverage >= min_coverage:
            # Classify the alignment using PAFRecord method
            classification = record.classify_alignment()
            classification_counts[classification] += 1

            # Add classification and overlap tags to the record
            record.tags["ct"] = classification
            record.tags["ov"] = overlap
            record.tags["ts"] = record.target_sample
            record.tags["qs"] = record.query_sample
            record.tags["th"] = record.target_haplotype
            record.tags["qh"] = record.query_haplotype

            filtered_records.append(record)
            passed += 1
        else:
            logging.warning(
                f"Filtered alignment: identity={identity:.2f}%, coverage={coverage:.2f}%"
            )

    # Write filtered records back to PAF
    PAFReader.write_paf(filtered_records, out_paf_file)

    # Log statistics
    with open(snakemake.log[0], "w") as log:
        log.write(f"Total alignments: {total}\n")
        log.write(f"Passed filters: {passed}\n")
        if total > 0:
            log.write(
                f"Filtered out: {total - passed} ({(total-passed)/total*100:.1f}%)\n"
            )
        log.write(f"\nClassification counts:\n")
        log.write(f"  Self alignments: {classification_counts['self']}\n")
        log.write(f"  Allelic alignments: {classification_counts['allelic']}\n")
        log.write(f"  Paralog alignments: {classification_counts['paralog']}\n")
        log.write(f"  Ortholog alignments: {classification_counts['ortholog']}\n")
        log.write(f"  Unknown: {classification_counts['unknown']}\n")


if __name__ == "__main__":
    main()
