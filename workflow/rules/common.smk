import pandas as pd
from pathlib import Path


# Load configuration
configfile: "config/config.yaml"


# Set default configuration values if not specified in config file
config.setdefault("output_dir", "results")
config.setdefault("temp_dir", "temp")

# Alignment defaults
config.setdefault("alignment", {})
config["alignment"].setdefault("min_identity", 80)
config["alignment"].setdefault("min_coverage", 80)
config["alignment"].setdefault("tool", "minimap2")
config["alignment"].setdefault("threads", 8)
config["alignment"].setdefault("minimap_preset", "asm20")
config["alignment"].setdefault("max_secondary", 100)
config["alignment"].setdefault("secondary_threshold", 0.8)
config["alignment"].setdefault("kmer_size", 15)
config["alignment"].setdefault("window_size", 10)

# Clustering defaults
config.setdefault("clustering", {})
config["clustering"].setdefault("allelic_distance", 10000)
config["clustering"].setdefault("reciprocal_overlap", 0.5)
config["clustering"].setdefault("method", "graph")

# Filtering defaults
config.setdefault("filtering", {})
config["filtering"].setdefault("min_re_length", 50)
config["filtering"].setdefault("max_re_length", 10000)
config["filtering"].setdefault("skip_repeats", False)
config["filtering"].setdefault("slop", 2000)
# sd_bed: Optional path to segmental duplication BED file for filtering
# If provided, REs overlapping SDs will be excluded


# Helper function for memory allocation with retries
def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 8
    return attempt * 1024 * 16


# Load sample table
def load_sample_table():
    """Load and validate the sample table"""
    sample_file = config["sample_table"]
    df = pd.read_csv(sample_file, sep="\t")

    required_cols = ["sample", "haplotype", "assembly", "re_bed"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Sample table missing required column: {col}")

    # Create a unique identifier for each sample-haplotype combination
    df["sample_id"] = df["sample"] + "_" + df["haplotype"]

    return df


SAMPLES = load_sample_table()


# Helper functions to get file paths
def get_assembly(wildcards):
    """Get assembly path for a given sample_id"""
    row = SAMPLES[SAMPLES["sample_id"] == wildcards.sample_id]
    return row["assembly"].values[0]


def get_re_bed(wildcards):
    """Get RE bed file path for a given sample_id"""
    row = SAMPLES[SAMPLES["sample_id"] == wildcards.sample_id]
    return row["re_bed"].values[0]


def get_all_sample_ids():
    """Get all unique sample IDs"""
    return SAMPLES["sample_id"].tolist()


def get_all_assemblies():
    """Get all assembly paths"""
    return SAMPLES["assembly"].tolist()


def get_all_re_beds():
    """Get all RE bed file paths"""
    return SAMPLES["re_bed"].tolist()


def get_all_assembly_fais():
    """Get all assembly FAI file paths"""
    return [asm + ".fai" for asm in SAMPLES["assembly"].tolist()]


def sd_filtering_enabled():
    """Check if SD filtering is enabled (sd_bed is configured)"""
    return config["filtering"].get("sd_bed") is not None


# Output directory helpers
OUTPUT = config.get("output_dir", "results")
TEMP = config.get("temp_dir", "temp")


def output_path(*args):
    """Construct output path"""
    return str(Path(OUTPUT) / Path(*args))


def temp_path(*args):
    """Construct temp path"""
    return str(Path(TEMP) / Path(*args))
