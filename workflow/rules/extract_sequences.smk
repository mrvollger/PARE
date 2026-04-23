"""
Extract per-haplotype_sample FASTAs and BEDs from the augmented union BED.

One rule produces all per-haplotype_sample outputs in a single pass (reading
the bed once), plus the global re_index TSV used downstream for expanding
paralog edges from the haplotype_sample level back to per-sample nodes.
"""


rule extract_union_sequences:
    """
    Parse the union BED and emit slopped/unslopped FASTAs and BEDs per
    haplotype_sample (only is_primary_sample == TRUE rows), plus a single
    re_index.tsv.gz covering every (sample_id x consensus_peak_id) row.
    """
    input:
        bed=config["union_bed"],
    output:
        fastas=temp(
            expand(
                temp_path("sequences/{haplotype_sample}.slop.fa"),
                haplotype_sample=get_haplotype_samples(),
            )
        ),
        fastas_unslop=temp(
            expand(
                temp_path("sequences/{haplotype_sample}.unslopped.fa"),
                haplotype_sample=get_haplotype_samples(),
            )
        ),
        beds=temp(
            expand(
                temp_path("sequences/{haplotype_sample}.bed"),
                haplotype_sample=get_haplotype_samples(),
            )
        ),
        slop_beds=temp(
            expand(
                temp_path("sequences/{haplotype_sample}.slop.bed"),
                haplotype_sample=get_haplotype_samples(),
            )
        ),
        re_index=output_path("re_index.tsv.gz"),
    log:
        "logs/extract_union_sequences.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        runtime=60,
    params:
        haplotype_samples=get_haplotype_samples(),
        min_len=config["filtering"]["min_re_length"],
        max_len=config["filtering"]["max_re_length"],
        drop_sd=config["filtering"]["drop_sd"],
    script:
        "../scripts/extract_union_sequences.py"


rule merge_unslopped_beds:
    """
    Concatenate all per-haplotype_sample unslopped BEDs into one sorted file.
    Used by the graph rules for paf-vs-RE intersection.
    """
    input:
        beds=expand(
            temp_path("sequences/{haplotype_sample}.bed"),
            haplotype_sample=get_haplotype_samples(),
        ),
    output:
        bed=temp(temp_path("merged_re_sequences.bed")),
    log:
        "logs/merge_beds/all.log",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """
