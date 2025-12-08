#!/usr/bin/env python3
"""
Filter PAF alignments based on identity and coverage thresholds.
Classify alignments as self, paralogous, or allelic.
Output PAF format with filtered alignments and classification tags.
"""

import sys
import re
import logging


def parse_re_id(re_id):
    """
    Parse RE ID to extract chromosome and coordinates.
    Format: chrom_start_end
    Example: chr20_MATERNAL_13652400_13652567
    """
    # Match pattern: chrom_start_end (everything before last two underscores is chrom)
    match = re.match(r'^(.+)_(\d+)_(\d+)$', re_id)
    if not match:
        return None

    return {
        'chrom': match.group(1),
        'start': int(match.group(2)),
        'end': int(match.group(3)),
    }


def get_haplotype_type(chrom_name):
    """
    Determine haplotype type from chromosome name.
    Returns: 'hap1', 'hap2', or 'unknown'
    """
    chrom_lower = chrom_name.lower()
    if 'paternal' in chrom_lower or '#1#' in chrom_lower or '_h1_' in chrom_lower:
        return 'hap1'
    elif 'maternal' in chrom_lower or '#2#' in chrom_lower or '_h2_' in chrom_lower:
        return 'hap2'
    return 'unknown'


def classify_alignment(query_chrom, target_chrom, query_start, query_end, target_start, target_end, re_info):
    """
    Classify alignment as self, paralog (same haplotype), or allelic (different haplotype).

    PAF structure (WITHOUT rb invert):
      - Query = Full chromosome (after adjustment from adjust_paf_for_slop.py)
      - Target = Where the RE aligned in the target assembly
      - id:Z = Original RE location (chrom_start_end)

    The query coordinates represent positions on the full chromosome.
    The id:Z tag tells us the original RE location.

    Args:
        query_chrom: Query chromosome (full chromosome after adjustment)
        target_chrom: Target chromosome name where the alignment mapped
        query_start: Start position on query (on full chromosome)
        query_end: End position on query (on full chromosome)
        target_start: Start position of alignment on target
        target_end: End position of alignment on target
        re_info: Dictionary with RE origin information (chrom, start, end from id:Z)

    Returns:
        str: One of 'self', 'allelic', or 'paralog'
    """
    # Self alignment: same chromosome and same positions
    if query_chrom == target_chrom:
        if query_start == target_start and query_end == target_end:
            return 'self'
        else:
            # Same chromosome but different position = paralog
            return 'paralog'

    # Different chromosomes - check haplotype
    query_hap = get_haplotype_type(query_chrom)
    target_hap = get_haplotype_type(target_chrom)

    # Allelic alignment: Different haplotypes
    if query_hap != 'unknown' and target_hap != 'unknown' and query_hap != target_hap:
        return 'allelic'

    # Paralog: Different chromosomes, same haplotype or unknown haplotype
    return 'paralog'


def parse_paf_line(line):
    """Parse a PAF format line and extract all relevant information"""
    fields = line.strip().split("\t")
    if len(fields) < 12:
        return None

    # Basic PAF fields
    query_chrom = fields[0]
    query_len = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])
    target_chrom = fields[5]
    target_len = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])

    # Extract RE information from id:Z: tag
    re_id = None
    re_info = None
    for field in fields[12:]:
        if field.startswith("id:Z:"):
            re_id = field[5:]
            re_info = parse_re_id(re_id)
            break

    if re_id is None:
        raise ValueError(f"id:Z: tag not found in PAF line: {line[:100]}...")

    if re_info is None:
        raise ValueError(f"Failed to parse RE ID: {re_id}")

    re_length = re_info['end'] - re_info['start']

    aln = {
        "query_chrom": query_chrom,
        "query_start": query_start,
        "query_end": query_end,
        "target_chrom": target_chrom,
        "target_start": target_start,
        "target_end": target_end,
        "nmatch": int(fields[9]),  # Number of matching bases
        "alen": int(fields[10]),  # Alignment block length
        "re_length": re_length,  # Original RE length
        "re_info": re_info,  # Parsed RE information
    }

    return aln


def calculate_metrics(aln):
    """Calculate identity, coverage, and overlap metrics"""
    # Coverage is calculated as number of matching bases / original RE length
    coverage = aln["nmatch"] / aln["re_length"] * 100 if aln["re_length"] > 0 else 0
    identity = aln["nmatch"] / aln["alen"] * 100 if aln["alen"] > 0 else 0

    # Calculate overlap between query and target positions
    # Overlap is 0 if chromosomes are different
    if aln["query_chrom"] == aln["target_chrom"]:
        overlap_start = max(aln["query_start"], aln["target_start"])
        overlap_end = min(aln["query_end"], aln["target_end"])
        overlap = max(0, overlap_end - overlap_start)
    else:
        overlap = 0

    return identity, coverage, overlap


def main():
    paf_file = snakemake.input.paf
    out_paf_file = snakemake.output.paf
    min_identity = snakemake.params.min_identity
    min_coverage = snakemake.params.min_coverage

    passed = 0
    total = 0
    classification_counts = {
        'self': 0,
        'allelic': 0,
        'paralog': 0,
        'unknown': 0
    }

    with open(paf_file, "r") as paf_in, open(out_paf_file, "w") as paf_out:
        for line in paf_in:
            if line.startswith("#"):
                paf_out.write(line)
                continue

            total += 1
            aln = parse_paf_line(line)
            if aln is None:
                continue

            identity, coverage, overlap = calculate_metrics(aln)

            # Filter based on thresholds
            if identity >= min_identity and coverage >= min_coverage:
                # Classify the alignment
                classification = classify_alignment(
                    query_chrom=aln['query_chrom'],
                    target_chrom=aln['target_chrom'],
                    query_start=aln['query_start'],
                    query_end=aln['query_end'],
                    target_start=aln['target_start'],
                    target_end=aln['target_end'],
                    re_info=aln['re_info']
                )

                classification_counts[classification] += 1

                # Add classification and overlap tags to PAF line
                output_line = line.rstrip() + f"\tct:Z:{classification}\tov:i:{overlap}\n"
                paf_out.write(output_line)
                passed += 1
            else:
                logging.warning(
                    f"Filtered alignment: identity={identity:.2f}%, coverage={coverage:.2f}%\n{line.strip()}"
                )

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
        log.write(f"  Unknown: {classification_counts['unknown']}\n")


if __name__ == "__main__":
    main()
