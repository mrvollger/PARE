#!/usr/bin/env python3
"""
PAF (Pairwise Alignment Format) parser class for working with alignment files.

Supports parsing standard PAF fields and custom tags (id:Z:, ct:Z:, ov:i:).
Can convert to pandas DataFrame for analysis.
"""

import re
from dataclasses import dataclass
from typing import Optional, Dict, Any, List
import pandas as pd


def parse_contig_name(contig_name: str) -> Dict[str, str]:
    """
    Parse contig name to extract sample ID and haplotype information.

    Supports various naming conventions:
    - HG00097#1#CM094075.1 -> sample: HG00097, haplotype: hap1
    - HG00097#2#CM094075.1 -> sample: HG00097, haplotype: hap2
    - chr1_PATERNAL -> sample: HG002, haplotype: hap1, chrom: chr1
    - chr1_MATERNAL -> sample: HG002, haplotype: hap2, chrom: chr1

    Args:
        contig_name: Contig/chromosome name

    Returns:
        Dictionary with 'sample', 'haplotype', 'chrom' keys
    """
    result = {
        'sample': 'unknown',
        'haplotype': 'unknown',
        'chrom': contig_name
    }

    # Pattern 1: Sample#hap#contig format (e.g., HG00097#1#CM094075.1)
    match = re.match(r'^([^#]+)#([12])#(.+)$', contig_name)
    if match:
        result['sample'] = match.group(1)
        result['haplotype'] = f"hap{match.group(2)}"
        result['chrom'] = match.group(3)
        return result

    # Pattern 2: Check for PATERNAL/MATERNAL or h1/h2 indicators
    if "_PATERNAL" in contig_name:
        result['sample'] = 'HG002'  # Default sample for this pattern
        result['haplotype'] = 'hap1'
        result['chrom'] = contig_name.replace("_PATERNAL", "")
    elif "_MATERNAL" in contig_name:
        result['sample'] = 'HG002'  # Default sample for this pattern
        result['haplotype'] = 'hap2'
        result['chrom'] = contig_name.replace("_MATERNAL", "")

    return result



@dataclass
class PAFRecord:
    """
    Represents a single PAF alignment record.

    Standard PAF fields (12 required):
        query_name, query_length, query_start, query_end, strand,
        target_name, target_length, target_start, target_end,
        num_matches, alignment_block_length, mapping_quality

    Common tags:
        id:Z: - Original RE ID (chrom_start_end)
        ct:Z: - Classification (self/allelic/paralog)
        ov:i: - Overlap between query and target positions
        cg:Z: - CIGAR string
    """
    # Standard PAF fields
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    num_matches: int
    alignment_block_length: int
    mapping_quality: int

    # Optional tags
    tags: Dict[str, Any]

    @classmethod
    def from_line(cls, line: str) -> Optional['PAFRecord']:
        """Parse a PAF line into a PAFRecord object."""
        if line.startswith('#') or not line.strip():
            return None

        fields = line.strip().split('\t')
        if len(fields) < 12:
            return None

        try:
            # Parse standard fields
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
            alignment_block_length = int(fields[10])
            mapping_quality = int(fields[11])

            # Parse tags (field 12 onwards)
            tags = {}
            for field in fields[12:]:
                match = re.match(r'([^:]+):([^:]+):(.+)', field)
                if match:
                    tag_name = match.group(1)
                    tag_type = match.group(2)
                    tag_value = match.group(3)

                    # Convert based on type
                    if tag_type == 'i':
                        tags[tag_name] = int(tag_value)
                    elif tag_type == 'f':
                        tags[tag_name] = float(tag_value)
                    else:  # Z, A, etc.
                        tags[tag_name] = tag_value

            return cls(
                query_name=query_name,
                query_length=query_length,
                query_start=query_start,
                query_end=query_end,
                strand=strand,
                target_name=target_name,
                target_length=target_length,
                target_start=target_start,
                target_end=target_end,
                num_matches=num_matches,
                alignment_block_length=alignment_block_length,
                mapping_quality=mapping_quality,
                tags=tags
            )
        except (ValueError, IndexError):
            return None

    def to_dict(self) -> Dict[str, Any]:
        """Convert PAF record to dictionary."""
        d = {
            'query_name': self.query_name,
            'query_length': self.query_length,
            'query_start': self.query_start,
            'query_end': self.query_end,
            'strand': self.strand,
            'target_name': self.target_name,
            'target_length': self.target_length,
            'target_start': self.target_start,
            'target_end': self.target_end,
            'num_matches': self.num_matches,
            'alignment_block_length': self.alignment_block_length,
            'mapping_quality': self.mapping_quality,
        }
        # Add tags as separate columns
        d.update(self.tags)

        # Add computed sample/haplotype columns
        d['query_sample'] = self.query_sample
        d['query_haplotype'] = self.query_haplotype
        d['target_sample'] = self.target_sample
        d['target_haplotype'] = self.target_haplotype

        return d

    def to_paf_line(self) -> str:
        """
        Convert PAF record back to PAF format string.

        Returns:
            PAF format line with all fields and tags
        """
        # Standard 12 fields
        fields = [
            self.query_name,
            str(self.query_length),
            str(self.query_start),
            str(self.query_end),
            self.strand,
            self.target_name,
            str(self.target_length),
            str(self.target_start),
            str(self.target_end),
            str(self.num_matches),
            str(self.alignment_block_length),
            str(self.mapping_quality),
        ]

        # Add tags
        for tag_name, tag_value in sorted(self.tags.items()):
            # Determine tag type
            if isinstance(tag_value, int):
                tag_type = 'i'
            elif isinstance(tag_value, float):
                tag_type = 'f'
            else:
                tag_type = 'Z'

            fields.append(f"{tag_name}:{tag_type}:{tag_value}")

        return '\t'.join(fields)

    @property
    def re_id(self) -> Optional[str]:
        """Get the RE ID from the id:Z: tag."""
        return self.tags.get('id')

    @property
    def classification(self) -> Optional[str]:
        """Get the classification from the ct:Z: tag."""
        return self.tags.get('ct')

    @property
    def overlap(self) -> Optional[int]:
        """Get the overlap from the ov:i: tag."""
        return self.tags.get('ov')

    @property
    def cigar(self) -> Optional[str]:
        """Get the CIGAR string from the cg:Z: tag."""
        return self.tags.get('cg')

    def calculate_identity(self) -> float:
        """Calculate alignment identity as percentage."""
        if self.alignment_block_length == 0:
            return 0.0
        return (self.num_matches / self.alignment_block_length) * 100

    def calculate_query_coverage(self) -> float:
        """Calculate query coverage as percentage."""
        if self.query_length == 0:
            return 0.0
        query_aln_len = self.query_end - self.query_start
        return (query_aln_len / self.query_length) * 100

    @property
    def query_sample(self) -> str:
        """Extract sample ID from query name."""
        return parse_contig_name(self.query_name)['sample']

    @property
    def query_haplotype(self) -> str:
        """Extract haplotype from query name."""
        return parse_contig_name(self.query_name)['haplotype']

    @property
    def query_chrom_parsed(self) -> str:
        """Extract parsed chromosome name from query name."""
        return parse_contig_name(self.query_name)['chrom']

    @property
    def target_sample(self) -> str:
        """Extract sample ID from target name."""
        return parse_contig_name(self.target_name)['sample']

    @property
    def target_haplotype(self) -> str:
        """Extract haplotype from target name."""
        return parse_contig_name(self.target_name)['haplotype']

    @property
    def target_chrom_parsed(self) -> str:
        """Extract parsed chromosome name from target name."""
        return parse_contig_name(self.target_name)['chrom']

    def shares_re_id(self, other: 'PAFRecord') -> bool:
        """
        Check if two PAF records share the same RE ID.

        Args:
            other: Another PAFRecord to compare with

        Returns:
            True if both records have the same re_id (id:Z: tag)
        """
        return self.re_id is not None and self.re_id == other.re_id

    def query_overlap(self, other: 'PAFRecord') -> int:
        """
        Calculate overlap in query coordinates with another PAF record.

        Args:
            other: Another PAFRecord to compare with

        Returns:
            Number of overlapping bases in query coordinates (0 if different chromosomes)
        """
        # Only calculate overlap if on same query chromosome
        if self.query_name != other.query_name:
            return 0

        overlap_start = max(self.query_start, other.query_start)
        overlap_end = min(self.query_end, other.query_end)
        return max(0, overlap_end - overlap_start)

    def target_overlap(self, other: 'PAFRecord') -> int:
        """
        Calculate overlap in target coordinates with another PAF record.

        Args:
            other: Another PAFRecord to compare with

        Returns:
            Number of overlapping bases in target coordinates (0 if different chromosomes)
        """
        # Only calculate overlap if on same target chromosome
        if self.target_name != other.target_name:
            return 0

        overlap_start = max(self.target_start, other.target_start)
        overlap_end = min(self.target_end, other.target_end)
        return max(0, overlap_end - overlap_start)

    def has_query_overlap(self, other: 'PAFRecord', min_overlap: int = 1) -> bool:
        """
        Check if two PAF records have overlapping query coordinates.

        Args:
            other: Another PAFRecord to compare with
            min_overlap: Minimum number of bases to consider as overlapping (default: 1)

        Returns:
            True if query regions overlap by at least min_overlap bases
        """
        return self.query_overlap(other) >= min_overlap

    def has_target_overlap(self, other: 'PAFRecord', min_overlap: int = 1) -> bool:
        """
        Check if two PAF records have overlapping target coordinates.

        Args:
            other: Another PAFRecord to compare with
            min_overlap: Minimum number of bases to consider as overlapping (default: 1)

        Returns:
            True if target regions overlap by at least min_overlap bases
        """
        return self.target_overlap(other) >= min_overlap

    def reciprocal_query_overlap(self, other: 'PAFRecord') -> float:
        """
        Calculate reciprocal overlap in query coordinates.

        Returns:
            Fraction of overlap relative to the shorter query span (0.0 to 1.0)
        """
        overlap = self.query_overlap(other)
        if overlap == 0:
            return 0.0

        len1 = self.query_end - self.query_start
        len2 = other.query_end - other.query_start
        min_len = min(len1, len2)

        if min_len == 0:
            return 0.0

        return overlap / min_len

    def reciprocal_target_overlap(self, other: 'PAFRecord') -> float:
        """
        Calculate reciprocal overlap in target coordinates.

        Returns:
            Fraction of overlap relative to the shorter target span (0.0 to 1.0)
        """
        overlap = self.target_overlap(other)
        if overlap == 0:
            return 0.0

        len1 = self.target_end - self.target_start
        len2 = other.target_end - other.target_start
        min_len = min(len1, len2)

        if min_len == 0:
            return 0.0

        return overlap / min_len


class PAFReader:
    """
    Reader for PAF files with support for iteration and DataFrame conversion.

    Usage:
        # Iterate through records
        reader = PAFReader('alignments.paf')
        for record in reader:
            print(record.query_name, record.target_name)

        # Convert to DataFrame
        df = reader.to_dataframe()
        print(df[['query_name', 'target_name', 'ct']])
    """

    def __init__(self, filename: str):
        """Initialize PAF reader with filename."""
        self.filename = filename

    def __iter__(self):
        """Iterate through PAF records."""
        with open(self.filename, 'r') as f:
            for line in f:
                record = PAFRecord.from_line(line)
                if record:
                    yield record

    def to_dataframe(self) -> pd.DataFrame:
        """
        Read entire PAF file into a pandas DataFrame.

        Returns:
            DataFrame with columns for all PAF fields and tags.
        """
        records = []
        for record in self:
            records.append(record.to_dict())

        if not records:
            return pd.DataFrame()

        return pd.DataFrame(records)

    def filter_by_classification(self, classification: str) -> List[PAFRecord]:
        """Get all records with a specific classification."""
        return [record for record in self if record.classification == classification]

    def filter_by_quality(self, min_mapq: int = 0) -> List[PAFRecord]:
        """Get all records with mapping quality >= threshold."""
        return [record for record in self if record.mapping_quality >= min_mapq]

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> List[PAFRecord]:
        """
        Convert a pandas DataFrame back to PAFRecord objects.

        Any columns not in the standard PAF fields will be added as tags.

        Args:
            df: DataFrame with PAF data (from to_dataframe())

        Returns:
            List of PAFRecord objects
        """
        # Standard PAF column names
        standard_cols = {
            'query_name', 'query_length', 'query_start', 'query_end', 'strand',
            'target_name', 'target_length', 'target_start', 'target_end',
            'num_matches', 'alignment_block_length', 'mapping_quality'
        }

        records = []
        for _, row in df.iterrows():
            # Extract standard fields
            standard_fields = {col: row[col] for col in standard_cols if col in df.columns}

            # Extract tags (any non-standard columns)
            tags = {col: row[col] for col in df.columns if col not in standard_cols}

            # Create record
            record = PAFRecord(
                query_name=standard_fields.get('query_name', ''),
                query_length=int(standard_fields.get('query_length', 0)),
                query_start=int(standard_fields.get('query_start', 0)),
                query_end=int(standard_fields.get('query_end', 0)),
                strand=standard_fields.get('strand', '+'),
                target_name=standard_fields.get('target_name', ''),
                target_length=int(standard_fields.get('target_length', 0)),
                target_start=int(standard_fields.get('target_start', 0)),
                target_end=int(standard_fields.get('target_end', 0)),
                num_matches=int(standard_fields.get('num_matches', 0)),
                alignment_block_length=int(standard_fields.get('alignment_block_length', 0)),
                mapping_quality=int(standard_fields.get('mapping_quality', 0)),
                tags=tags
            )
            records.append(record)

        return records

    @staticmethod
    def write_paf(records: List[PAFRecord], filename: str):
        """
        Write PAFRecord objects to a PAF file.

        Args:
            records: List of PAFRecord objects
            filename: Output PAF filename
        """
        with open(filename, 'w') as f:
            for record in records:
                f.write(record.to_paf_line() + '\n')

    @staticmethod
    def dataframe_to_paf(df: pd.DataFrame, filename: str):
        """
        Convert a pandas DataFrame to PAF format and write to file.

        Any columns not in the standard PAF fields will be added as tags.

        Args:
            df: DataFrame with PAF data
            filename: Output PAF filename

        Example:
            reader = PAFReader('input.paf')
            df = reader.to_dataframe()
            # Add new column
            df['new_tag'] = df['num_matches'] * 2
            # Write back to PAF with new tag
            PAFReader.dataframe_to_paf(df, 'output.paf')
        """
        records = PAFReader.from_dataframe(df)
        PAFReader.write_paf(records, filename)


# Example usage
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='PAF file parser and analysis tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic stats
  python paf.py alignments.paf

  # Test overlap functions
  python paf.py alignments.paf --test-overlap

  # Show first N records
  python paf.py alignments.paf --head 5
        """
    )
    parser.add_argument('paf_file', help='PAF file to parse')
    parser.add_argument('--test-overlap', action='store_true',
                        help='Test overlap functions on records')
    parser.add_argument('--head', type=int, metavar='N',
                        help='Show first N records in detail')
    parser.add_argument('--test-roundtrip', metavar='OUTPUT',
                        help='Test DataFrame round-trip conversion and write to OUTPUT')

    args = parser.parse_args()

    reader = PAFReader(args.paf_file)
    df = reader.to_dataframe()

    print(f"Loaded {len(df)} alignments")
    print(f"\nColumns: {', '.join(df.columns)}")

    if 'ct' in df.columns:
        print("\nClassification counts:")
        print(df['ct'].value_counts())

    if args.head:
        print(f"\nFirst {args.head} records:")
        for i in range(min(args.head, len(df))):
            print(f"\n--- Record {i+1} ---")
            print(df.iloc[i])

    # Test overlap functions if requested
    if args.test_overlap:
        print("\n" + "="*60)
        print("Testing overlap functions")
        print("="*60)

        records = list(reader)
        if len(records) >= 2:
            rec1, rec2 = records[0], records[1]

            print(f"\nRecord 1: {rec1.query_name}:{rec1.query_start}-{rec1.query_end}")
            print(f"          -> {rec1.target_name}:{rec1.target_start}-{rec1.target_end}")
            print(f"          RE ID: {rec1.re_id}, Classification: {rec1.classification}")

            print(f"\nRecord 2: {rec2.query_name}:{rec2.query_start}-{rec2.query_end}")
            print(f"          -> {rec2.target_name}:{rec2.target_start}-{rec2.target_end}")
            print(f"          RE ID: {rec2.re_id}, Classification: {rec2.classification}")

            print(f"\nComparison:")
            print(f"  shares_re_id: {rec1.shares_re_id(rec2)}")
            print(f"  query_overlap: {rec1.query_overlap(rec2)} bp")
            print(f"  target_overlap: {rec1.target_overlap(rec2)} bp")
            print(f"  has_query_overlap: {rec1.has_query_overlap(rec2)}")
            print(f"  has_target_overlap: {rec1.has_target_overlap(rec2)}")
            print(f"  reciprocal_query_overlap: {rec1.reciprocal_query_overlap(rec2):.3f}")
            print(f"  reciprocal_target_overlap: {rec1.reciprocal_target_overlap(rec2):.3f}")

            # Find records with same RE ID
            print(f"\n" + "="*60)
            print("Finding records with same RE ID")
            print("="*60)
            re_groups = {}
            for rec in records[:100]:  # Test first 100
                if rec.re_id:
                    if rec.re_id not in re_groups:
                        re_groups[rec.re_id] = []
                    re_groups[rec.re_id].append(rec)

            multi_hit_res = {k: v for k, v in re_groups.items() if len(v) > 1}
            print(f"\nFound {len(multi_hit_res)} REs with multiple alignments (from first 100)")

            if multi_hit_res:
                example_re = list(multi_hit_res.keys())[0]
                example_recs = multi_hit_res[example_re]
                print(f"\nExample: RE {example_re} has {len(example_recs)} alignments:")
                for i, rec in enumerate(example_recs, 1):
                    print(f"  {i}. {rec.target_name}:{rec.target_start}-{rec.target_end} [{rec.classification}]")

                if len(example_recs) >= 2:
                    r1, r2 = example_recs[0], example_recs[1]
                    print(f"\nComparing first two alignments:")
                    print(f"  Target overlap: {r1.target_overlap(r2)} bp")
                    print(f"  Reciprocal target overlap: {r1.reciprocal_target_overlap(r2):.3f}")

    # Test round-trip conversion if requested
    if args.test_roundtrip:
        print("\n" + "="*60)
        print(f"Testing DataFrame round-trip conversion")
        print("="*60)

        # Read to DataFrame
        df = reader.to_dataframe()
        print(f"\nOriginal DataFrame: {len(df)} records")

        # Add a new column to test tag addition
        df['test_score'] = df['num_matches'] * 2
        df['identity_pct'] = (df['num_matches'] / df['alignment_block_length'] * 100).round(2)
        print(f"Added new columns: test_score, identity_pct")

        # Convert back to PAF and write
        PAFReader.dataframe_to_paf(df, args.test_roundtrip)
        print(f"\nWrote {len(df)} records to: {args.test_roundtrip}")

        # Verify by reading back
        verify_reader = PAFReader(args.test_roundtrip)
        verify_df = verify_reader.to_dataframe()
        print(f"Verified: {len(verify_df)} records read back")

        # Check new tags are present
        if 'test_score' in verify_df.columns and 'identity_pct' in verify_df.columns:
            print("✓ New tags successfully added to PAF")
            print(f"\nSample of new tags:")
            print(verify_df[['query_name', 'test_score', 'identity_pct']].head(3))
        else:
            print("✗ New tags not found in output PAF")
