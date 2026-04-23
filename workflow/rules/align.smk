"""
Per-haplotype_sample all-vs-all alignment of slopped RE sequences.

Each haplotype_sample = (Individual_ID, Haplotype) is aligned to itself to
discover paralog pairs at the consensus_peak_id level. Paralogs found in any
haplotype_sample are unioned downstream to produce the final consensus-peak
paralog graph.

Input queries are already in slop-local coordinates (0..slop_seq length). The
RE_ID is encoded in the sequence name (both query and target); downstream code
parses it back to sample_id / consensus_peak_id.
"""


rule align_within_haplotype_sample:
    """
    minimap2 all-vs-all on one haplotype_sample's slopped sequences.

    -X enables all-vs-all overlap mode (skips self-mappings and redundant
    A->B/B->A pairs). Parameters follow the existing RE-to-RE preset in
    all_by_all.smk (small k, dense minimizers).
    """
    input:
        fasta=temp_path("sequences/{haplotype_sample}.slop.fa"),
    output:
        paf=temp(temp_path("alignments/{haplotype_sample}.paf")),
    log:
        "logs/align/{haplotype_sample}.log",
    conda:
        "../envs/bfx.yml"
    threads: config["alignment"]["threads"]
    resources:
        mem_mb=get_mem_mb,
        runtime=240,
    params:
        kmer=config["alignment"]["kmer_size"],
        window=config["alignment"]["window_size"],
    shell:
        """
        minimap2 \
            -X \
            -c \
            --eqx \
            -k {params.kmer} \
            -w {params.window} \
            -n 2 \
            -m 20 \
            -s 0 \
            -t {threads} \
            {input.fasta} \
            {input.fasta} \
            > {output.paf} 2> {log}
        """


rule filter_alignments:
    """
    Filter by identity / coverage and add classification tags.

    All alignments are within a single haplotype_sample, so classification
    reduces to 'paralog' (different consensus_peak_id within same haplotype
    sample); self-mappings are already skipped by -X. The classifier still runs
    to be robust to future modes.
    """
    input:
        paf=rules.align_within_haplotype_sample.output.paf,
        re_index=rules.extract_union_sequences.output.re_index,
    output:
        paf=temp(temp_path("filtered/{haplotype_sample}.paf")),
    log:
        "logs/filter_alignments/{haplotype_sample}.log",
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


rule merge_filtered_pafs:
    """
    Concatenate all per-haplotype_sample filtered PAFs into one unified file.
    """
    input:
        pafs=expand(
            temp_path("filtered/{haplotype_sample}.paf"),
            haplotype_sample=get_haplotype_samples(),
        ),
    output:
        paf=output_path("filtered_alignments/all.paf"),
    log:
        "logs/merge_filtered_pafs.log",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=20,
    shell:
        """
        cat {input.pafs} > {output.paf} 2> {log}
        """
