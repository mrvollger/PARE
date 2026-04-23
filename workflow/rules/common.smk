"""
Shared configuration, helpers, and sample discovery for PARE.

PARE now consumes a single augmented union peaks BED (bgzipped TSV) produced
upstream by `add_slop_to_union_bed.py`. The BED provides:
  - consensus_peak_id  (cross-sample consensus ID)
  - sample_id          (composite, e.g. HG002_PS01015_PacBio_HG002_2)
  - Individual_ID, Haplotype
  - asm_chr, asm_start, asm_end, asm_slop_start, asm_slop_end
  - slop_seq           (±slop flanking sequence, pre-extracted from each assembly)
  - is_peak, is_primary_sample, is_SD_*

No per-sample assembly FASTAs are required.

Alignment is performed once per haplotype_sample = (Individual_ID, Haplotype),
using only is_primary_sample == TRUE rows. Non-primary rows share the same
underlying sequence and are kept only as graph nodes (not alignment inputs).

RE_ID format: {sample_id}__{consensus_peak_id}   (see scripts/union_bed.py)
"""

import csv
import gzip
import sys
from pathlib import Path


configfile: "config/config.yaml"


# --- defaults -----------------------------------------------------------------

config.setdefault("output_dir", "results")
config.setdefault("temp_dir", "temp")

# Alignment
config.setdefault("alignment", {})
config["alignment"].setdefault("min_identity", 80)
config["alignment"].setdefault("min_coverage", 80)
config["alignment"].setdefault("threads", 8)
config["alignment"].setdefault("kmer_size", 15)
config["alignment"].setdefault("window_size", 10)

# Clustering
config.setdefault("clustering", {})
config["clustering"].setdefault("reciprocal_overlap", 0.5)

# Filtering
config.setdefault("filtering", {})
config["filtering"].setdefault("min_re_length", 50)
config["filtering"].setdefault("max_re_length", 10000)
config["filtering"].setdefault("drop_sd", False)

# Pileup (optional: requires a separate ft_bam table)
config.setdefault("pileup", {})

# Required
if "union_bed" not in config:
    raise ValueError(
        "config is missing required key 'union_bed' (path to the augmented union peaks BED)"
    )


# --- helpers ------------------------------------------------------------------

OUTPUT = config["output_dir"]
TEMP = config["temp_dir"]


def output_path(*args):
    return str(Path(OUTPUT) / Path(*args))


def temp_path(*args):
    return str(Path(TEMP) / Path(*args))


def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 8
    return attempt * 1024 * 16


# --- load union BED once at workflow-parse time ------------------------------

_UNION_BED = config["union_bed"]

# A small sidecar TSV listing the samples in the union BED. Required so that
# we don't have to stream the (potentially multi-GB) union BED at workflow-
# parse time just to discover sample_ids.
#
# Columns (tab-separated, header required):
#   sample_id              e.g. HG002_PS01015_PacBio_HG002_2
#   Individual_ID          e.g. HG002
#   Haplotype              1 or 2
#   is_primary_sample      TRUE or FALSE
#
# One row per unique sample_id. You can generate it from the union BED once:
#
#   ( echo -e "sample_id\tIndividual_ID\tHaplotype\tis_primary_sample"; \
#     bgzip -dc <union.bed.gz> | \
#     awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next} \
#                 {print $c["sample_id"]"\t"$c["Individual_ID"]"\t"$c["Haplotype"]"\t"$c["is_primary_sample"]}' | \
#     sort -u ) > samples.tsv
#
if "samples_tsv" not in config:
    raise ValueError(
        "config is missing required key 'samples_tsv' (small TSV of sample_id / Individual_ID / Haplotype / is_primary_sample)"
    )
_SAMPLES_TSV = config["samples_tsv"]


def _open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def _load_samples_tsv(path):
    """Read the small samples sidecar; returns (haplotype_samples, sample_ids)."""
    haplotype_samples = set()
    sample_ids = set()
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "Individual_ID", "Haplotype", "is_primary_sample"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"samples_tsv {path} missing required columns: {sorted(missing)}"
            )
        for row in reader:
            sid = row["sample_id"].strip()
            if not sid:
                continue
            sample_ids.add(sid)
            if row["is_primary_sample"].strip().upper() == "TRUE":
                haplotype_samples.add(
                    f"{row['Individual_ID'].strip()}_{row['Haplotype'].strip()}"
                )
    return sorted(haplotype_samples), sorted(sample_ids)


HAPLOTYPE_SAMPLES, SAMPLE_IDS = _load_samples_tsv(_SAMPLES_TSV)


def get_haplotype_samples():
    """Ordered list of haplotype_sample keys (e.g. HG002_1, HG002_2, ...)."""
    return HAPLOTYPE_SAMPLES


def get_all_sample_ids():
    """All unique sample_ids (primary and non-primary)."""
    return SAMPLE_IDS


# --- optional ft_bam side-table (for pileup) ---------------------------------

def _load_ft_bam_table():
    """
    Optional TSV mapping sample_id -> ft_bam path. Enables the pileup rules.

    Format:
        sample_id\tft_bam
        HG002_HG002_1\t/path/to/HG002.dsa.bam
        ...
    """
    path = config.get("ft_bam_table")
    if not path:
        return {}
    df = pl.read_csv(path, separator="\t", has_header=True)
    need = {"sample_id", "ft_bam"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"ft_bam_table at {path} missing columns: {missing}")
    return dict(zip(df.get_column("sample_id"), df.get_column("ft_bam")))


_FT_BAM = _load_ft_bam_table()


def has_ft_bam(sample_id):
    return sample_id in _FT_BAM


def get_ft_bam(wildcards):
    return _FT_BAM[wildcards.sample_id]


def get_sample_ids_with_ft_bam():
    return [sid for sid in SAMPLE_IDS if sid in _FT_BAM]


def ft_pileup_enabled():
    return len(_FT_BAM) > 0
