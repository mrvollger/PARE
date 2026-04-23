"""
Optional pileup of ft_bam coverage on clustered RE regions.

Enabled when `ft_bam_table` (a TSV mapping sample_id -> ft_bam path) is set
in the config. The annotated_res.bed from graph.smk provides per-sample_id
rows via column 7 (sample_id).
"""


rule split_annotated_res_by_sample:
    """Extract this sample_id's REs from annotated_res.bed."""
    input:
        bed=rules.paf_to_graphml.output.annotated_res,
    output:
        bed=temp(temp_path("pileup/{sample_id}.res.bed")),
    log:
        "logs/pileup/split_res_{sample_id}.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    params:
        sample_id=lambda wc: wc.sample_id,
    shell:
        r"""
        awk -F'\t' -v sid={params.sample_id} 'BEGIN{{OFS="\t"}} $1 !~ /^#/ && $7==sid {{print $1,$2,$3,$4}}' \
            {input.bed} > {output.bed} 2> {log}
        """


rule ft_pileup:
    input:
        bed=rules.split_annotated_res_by_sample.output.bed,
        bam=get_ft_bam,
        bai=lambda wc: get_ft_bam(wc) + ".bai",
    output:
        pileup=temp(temp_path("pileup/{sample_id}.pileup.bed")),
    log:
        "logs/pileup/ft_pileup_{sample_id}.log",
    conda:
        "../envs/bfx.yml"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=120,
    params:
        rolling_max_window=config["pileup"].get("rolling_max_window", 200),
    shell:
        """
        ft pileup \
            --bed {input.bed} \
            -t {threads} \
            --rolling-max {params.rolling_max_window} \
            --fiber-coverage \
            --haps \
            {input.bam} \
            2> {log} \
        > {output.pileup}
        """


rule filter_pileup_peaks:
    input:
        pileup=rules.ft_pileup.output.pileup,
    output:
        peaks=pipe(temp_path("pileup/{sample_id}.peaks.bed")),
    log:
        "logs/pileup/filter_peaks_{sample_id}.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    script:
        "../scripts/filter_pileup_peaks.py"


rule compress_peaks:
    input:
        bed=rules.filter_pileup_peaks.output.peaks,
    output:
        gz=output_path("pileup/{sample_id}.peaks.bed.gz"),
        index=output_path("pileup/{sample_id}.peaks.bed.gz.gzi"),
    log:
        "logs/pileup/compress_peaks_{sample_id}.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    shell:
        """
        bgzip -i -o {output.gz} {input.bed} 2> {log}
        """


rule merge_all_peaks:
    """Concatenate per-sample peaks into one sorted bgzipped + tabix'd file."""
    input:
        peaks=expand(
            output_path("pileup/{sample_id}.peaks.bed.gz"),
            sample_id=get_sample_ids_with_ft_bam(),
        ),
    output:
        bed=output_path("pileup/all_peaks.bed.gz"),
        tbi=output_path("pileup/all_peaks.bed.gz.tbi"),
    log:
        "logs/pileup/merge_all_peaks.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    shell:
        """
        (
            gzip -dc {input.peaks[0]} | head -1 || true
            for f in {input.peaks}; do
                gzip -dc "$f" | tail -n +2
            done | sort -k1,1 -k2,2n
        ) | bgzip -c > {output.bed} 2> {log}
        tabix -p bed {output.bed} 2>> {log}
        """


rule merge_peaks_with_annotations:
    """Join peaks (keyed by re_id) with annotated_res.bed to add cluster_id."""
    input:
        peaks=rules.merge_all_peaks.output.bed,
        annotated_res=rules.paf_to_graphml.output.annotated_res,
    output:
        bed=output_path("pileup/peaks_with_clusters.bed.gz"),
        tbi=output_path("pileup/peaks_with_clusters.bed.gz.tbi"),
    log:
        "logs/pileup/merge_peaks_with_annotations.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    script:
        "../scripts/merge_peaks_with_annotations.py"
