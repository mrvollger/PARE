"""
Rules for all-vs-all RE alignment to get pairwise sequence identity.

This module handles:
1. Running all-vs-all minimap2 alignment on RE sequences
2. Annotating alignments with cluster IDs and calculating identity

RE sequences (both annotated and implicit) are extracted in the graph.smk
paf_to_graphml rule, which outputs results/graphs/annotated_res.fa
"""


rule align_re_to_re:
    """
    All-vs-all alignment of RE sequences to get pairwise sequence identity.

    This aligns all RE sequences (annotated + implicit) against themselves to find
    direct RE-to-RE alignments with true sequence identity scores. Used for
    correlating sequence divergence with chromatin accessibility variance.

    Input includes both:
    - Annotated REs (from original RE BED files)
    - Implicit REs (discovered during graph construction)

    Key parameters for short sequences (50-10,000bp REs):
    - -X: All-vs-all mode (equivalent to -DP --dual=no --no-long-join)
          Skips self mappings and redundant A->B/B->A pairs
    - -k 11: Smaller k-mer for better sensitivity on short sequences
    - -w 5: Smaller window for denser minimizer sampling
    - -n 2: Allow chains with only 2 minimizers (short sequences)
    - -m 20: Lower chaining score threshold for short alignments
    - -s 0: No minimum DP score threshold to capture small REs
    - --eqx: Use =/X CIGAR operators for accurate identity calculation

    Note: -N and -p have no effect when -X is used (all chains retained)
    """
    input:
        fasta=rules.merge_annotated_res_fastas.output.fasta,
    output:
        paf=temp("temp/all_by_all/re_alignments.paf"),
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
    Annotate RE-to-RE alignments with cluster IDs and filter.

    - Parses RE IDs from sequence names (format: chrom_start_end or IMPLICIT_chrom_start_end)
    - Joins with cluster assignments from annotated_res.bed
    - Removes self-alignments (same RE to itself)
    - Calculates sequence identity from alignment stats
    - Outputs alignments with cluster info for downstream analysis
    """
    input:
        paf=rules.align_re_to_re.output.paf,
        annotated_res=rules.paf_to_graphml.output.annotated_res,
    output:
        paf="results/all_by_all/re_alignments_with_clusters.paf.gz",
        tsv="results/all_by_all/re_pairwise_identity.tsv.gz",
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
