"""
Rules for aligning RE sequences to all assemblies to find paralogs
"""


rule merge_all_re_sequences:
    """
    Merge all RE sequences from all samples into a single FASTA
    Sequence names already contain sample_id prefix from bedtools getfasta -name
    """
    input:
        fastas=expand("temp/sequences/{sample_id}.re_sequences.fa", sample_id=get_all_sample_ids()),
    output:
        fasta="temp/merged_re_sequences.fa",
    log:
        "logs/merge_sequences/all.log",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    shell:
        """
        cat {input.fastas} > {output.fasta} 2> {log}
        """


rule align_all_res_to_assembly:
    """
    Align all merged RE sequences to each target assembly using minimap2.

    Key parameters:
    -N: Report up to N secondary alignments (finds paralogs)
    -p: Secondary threshold - reports if score >= primary * p
    --eqx: Output = for match, X for mismatch (better identity calculation)
    -k/-w: k-mer and window size optimized for short sequences
    """
    input:
        query="temp/merged_re_sequences.fa",
        target=get_assembly,
    output:
        paf=temp("temp/alignments/all_REs_vs_{sample_id}.paf"),
    log:
        "logs/align/all_REs_vs_{sample_id}.log",
    conda:
        "../envs/bfx.yml"
    threads: config["alignment"]["threads"]
    resources:
        mem_mb=get_mem_mb,
        runtime=240,
    params:
        preset=config["alignment"]["minimap_preset"],
        max_sec=config["alignment"]["max_secondary"],
        sec_threshold=config["alignment"]["secondary_threshold"],
        kmer=config["alignment"]["kmer_size"],
        window=config["alignment"]["window_size"],
    shell:
        """
        minimap2 \
            -x {params.preset} \
            -c \
            --eqx \
            --secondary==yes \
            -s 20 \
            -N {params.max_sec} \
            -p {params.sec_threshold} \
            -t {threads} \
            -k {params.kmer} \
            -w {params.window} \
            {input.target} \
            {input.query} \
            > {output.paf} 2> {log}
        """


rule adjust_paf_coordinates:
    """
    Adjust PAF coordinates from extracted sequence space to original assembly space.
    Changes query name from extracted RE name to chromosome name,
    and adjusts coordinates to represent positions on the full chromosome.
    """
    input:
        paf="temp/alignments/all_REs_vs_{sample_id}.paf",
        fai="temp/merged_assemblies.fai",
    output:
        paf=temp("temp/alignments/all_REs_vs_{sample_id}.adjusted.paf"),
    log:
        "logs/adjust_paf/all_REs_vs_{sample_id}.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    script:
        "../scripts/adjust_paf_for_slop.py"


rule liftover_trim_alignments:
    """
    Use rustybam liftover to trim alignments back to original RE coordinates.
    Uses --qbed to trim based on query (RE) coordinates, not target assembly.
    This removes alignments that only match flanking sequence from the slop
    and trims both coordinates and CIGAR strings accordingly.
    """
    input:
        paf="temp/alignments/all_REs_vs_{sample_id}.adjusted.paf",
        bed="temp/merged_re_sequences.bed",
    output:
        paf=temp("temp/alignments/all_REs_vs_{sample_id}.trimmed.paf"),
    log:
        "logs/liftover/all_REs_vs_{sample_id}.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=60,
    shell:
        """
        rb liftover \
            --qbed \
            --bed {input.bed} \
            {input.paf} \
            | awk 'BEGIN{{OFS="\\t"}} {{$12=255; print}}' \
            | sort \
            | uniq \
            > {output.paf} 2> {log}
        """


rule filter_alignments:
    """
    Filter trimmed alignments based on identity and coverage thresholds.
    Classify each alignment as 'self', 'allelic', or 'paralog'.
    Output PAF format with filtered alignments and ct:Z: classification tag.
    """
    input:
        paf="temp/alignments/all_REs_vs_{sample_id}.trimmed.paf",
    output:
        paf="results/filtered_alignments/all_REs_vs_{sample_id}.paf",
    log:
        "logs/filter_alignments/all_REs_vs_{sample_id}.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    params:
        min_identity=config["alignment"]["min_identity"],
        min_coverage=config["alignment"]["min_coverage"],
    script:
        "../scripts/filter_paf.py"
