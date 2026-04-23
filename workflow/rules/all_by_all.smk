"""
Cross-haplotype_sample RE-to-RE identity matrix.

This takes the merged unslopped FASTA (one entry per primary sample_id x
consensus_peak_id) and runs all-vs-all minimap2 to get pairwise sequence
identity across haplotype_samples (the within-haplotype-sample alignments in
align.smk are for paralog discovery; this step is for cross-sample comparison
on the *called-peak* sequences, not flanking regions).
"""


rule merge_all_unslopped_fastas:
    """
    Concatenate all per-haplotype_sample unslopped FASTAs.
    """
    input:
        fastas=expand(
            temp_path("sequences/{haplotype_sample}.unslopped.fa"),
            haplotype_sample=get_haplotype_samples(),
        ),
    output:
        fasta=temp(temp_path("all_by_all/unslopped.fa")),
    log:
        "logs/all_by_all/merge_unslopped.log",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    shell:
        """
        cat {input.fastas} > {output.fasta} 2> {log}
        """


rule align_re_to_re:
    """
    All-vs-all minimap2 on the merged unslopped RE sequences.
    Small-k / dense minimizers are tuned for 50-10kb peak sequences.
    """
    input:
        fasta=rules.merge_all_unslopped_fastas.output.fasta,
    output:
        paf=temp(temp_path("all_by_all/re_alignments.paf")),
    log:
        "logs/all_by_all/align_re_to_re.log",
    conda:
        "../envs/bfx.yml"
    threads: config["alignment"]["threads"]
    resources:
        mem_mb=get_mem_mb,
        runtime=240,
    shell:
        """
        minimap2 \
            -X \
            -c \
            --eqx \
            -k 11 \
            -w 5 \
            -n 2 \
            -m 20 \
            -s 0 \
            -t {threads} \
            {input.fasta} \
            {input.fasta} \
            > {output.paf} 2> {log}
        """


rule annotate_re_to_re_alignments:
    """
    Join with cluster assignments and emit a deduped pairwise identity TSV.
    """
    input:
        paf=rules.align_re_to_re.output.paf,
        re_clusters=rules.paf_to_graphml.output.re_clusters,
        re_index=rules.extract_union_sequences.output.re_index,
    output:
        paf=output_path("all_by_all/re_alignments_with_clusters.paf.gz"),
        tsv=output_path("all_by_all/re_pairwise_identity.tsv.gz"),
    log:
        "logs/all_by_all/annotate_alignments.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=60,
    script:
        "../scripts/annotate_re_to_re.py"
