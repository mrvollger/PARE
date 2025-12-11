#!/usr/bin/env python3
"""
Adjust PAF coordinates to represent positions in the original assembly.

The query sequences were extracted with slop, but we want the PAF to show
positions as if we aligned regions of the original chromosome.

Input PAF query format:
  Name: chr20_PATERNAL_13360_13478::chr20_PATERNAL:12860-13978
  qlen: 1119 (slopped length)
  qstart/qend: positions within extracted sequence (0-1119)

Output PAF query format:
  Name: chr20_PATERNAL (original chromosome)
  qlen: 67035379 (full chromosome length from FAI)
  qstart/qend: positions on original chromosome (12860 + offset)
"""

import sys
import re


def load_fai_lengths(fai_files):
    """Load chromosome lengths from FAI files"""
    lengths = {}
    for fai_file in fai_files:
        with open(fai_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) >= 2:
                    chrom = fields[0]
                    length = int(fields[1])
                    lengths[chrom] = length
    return lengths


def parse_query_name(qname):
    """
    Parse query name to extract chromosome and coordinates.
    Format: chrom_unslop_start_unslop_end::chrom:slop_start-slop_end
    Returns: (chrom, slop_start, slop_end)
    """
    match = re.match(r"(.+)_\d+_\d+::(.+):(\d+)-(\d+)", qname)
    if not match:
        return None

    chrom = match.group(2)  # Use the chrom from the coordinates part
    slop_start = int(match.group(3))
    slop_end = int(match.group(4))

    return chrom, slop_start, slop_end


def adjust_paf_line(line, chrom_lengths, log):
    """Adjust a single PAF line to reference original chromosome coordinates"""
    fields = line.strip().split("\t")
    if len(fields) < 12:
        return None

    qname = fields[0]
    qlen = int(fields[1])
    qstart = int(fields[2])
    qend = int(fields[3])

    # Parse query name
    parsed = parse_query_name(qname)
    if parsed is None:
        log.write(f"WARNING: Could not parse query name: {qname}\n")
        return line  # Return unchanged

    chrom, slop_start, slop_end = parsed

    # Get chromosome length from FAI
    if chrom not in chrom_lengths:
        log.write(f"WARNING: Chromosome {chrom} not found in FAI files\n")
        return None

    chrom_len = chrom_lengths[chrom]

    # Adjust coordinates to original chromosome positions
    # qstart/qend are positions in extracted sequence (0-based)
    # Add slop_start to get position on original chromosome
    new_qstart = slop_start + qstart
    new_qend = slop_start + qend

    # Update fields
    fields[0] = chrom  # Change query name to chromosome
    fields[1] = str(chrom_len)  # Full chromosome length
    fields[2] = str(new_qstart)  # Adjusted start
    fields[3] = str(new_qend)  # Adjusted end

    return "\t".join(fields)


def main():
    paf_in = snakemake.input.paf
    paf_out = snakemake.output.paf
    log_file = snakemake.log[0]

    # Load FAI files to get chromosome lengths
    # For now, we'll get them from the input assembly FAI files
    # The rule should provide these as input
    fai_files = []
    if hasattr(snakemake.input, "fai"):
        if isinstance(snakemake.input.fai, list):
            fai_files = snakemake.input.fai
        else:
            fai_files = [snakemake.input.fai]

    chrom_lengths = {}
    if fai_files:
        chrom_lengths = load_fai_lengths(fai_files)

    total = 0
    adjusted = 0
    filtered = 0

    with open(paf_in, "r") as infile, open(paf_out, "w") as outfile, open(
        log_file, "w"
    ) as log:
        log.write("Adjusting PAF coordinates to original assembly positions\n\n")
        log.write(f"Loaded {len(chrom_lengths)} chromosomes from FAI files\n\n")

        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            total += 1
            adjusted_line = adjust_paf_line(line, chrom_lengths, log)

            if adjusted_line is None:
                filtered += 1
            else:
                adjusted += 1
                outfile.write(adjusted_line + "\n")

        log.write(f"\nSummary:\n")
        log.write(f"  Total alignments: {total}\n")
        log.write(f"  Adjusted: {adjusted}\n")
        log.write(f"  Filtered: {filtered}\n")


if __name__ == "__main__":
    main()
