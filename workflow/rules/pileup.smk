"""
Rules for running ft pileup on clustered RE regions to extract chromatin information.

These rules are conditional on ft_bam being present in the sample table.
"""


rule split_annotated_res_by_sample:
    """Split annotated REs BED by sample_id for per-sample pileup."""
    input:
        bed=rules.paf_to_graphml.output.annotated_res,
    output:
        bed=temp("temp/pileup/{sample_id}.res.bed"),
    log:
        "logs/pileup/split_res_{sample_id}.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    params:
        sample_id=lambda wc: wc.sample_id,
    run:
        with open(input.bed) as f_in, open(output.bed, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                chrom = fields[0]
                # Check if this RE belongs to this sample/haplotype
                if contig_matches_sample_id(chrom, params.sample_id):
                    # Output: chrom, start, end, re_id (use re_id as name for ft pileup)
                    f_out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\n")


rule ft_pileup:
    """Run ft pileup on RE regions to extract chromatin accessibility."""
    input:
        bed="temp/pileup/{sample_id}.res.bed",
        bam=get_ft_bam,
        bai=lambda wc: get_ft_bam(wc) + ".bai",
    output:
        pileup=temp("temp/pileup/{sample_id}.pileup.bed"),
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
    """Filter pileup to best peak per RE.

    For each continuous run of the same RE name:
    1. Find max score
    2. Break ties with max fire_coverage
    3. Break ties with most central position
    """
    input:
        pileup=rules.ft_pileup.output.pileup,
    output:
        peaks=pipe("temp/pileup/{sample_id}.peaks.bed"),
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
    """Compress peaks with bgzip and create index."""
    input:
        bed=rules.filter_pileup_peaks.output.peaks,
    output:
        gz="results/pileup/{sample_id}.peaks.bed.gz",
        index="results/pileup/{sample_id}.peaks.bed.gz.gzi",
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


rule index_pileup:
    """Index the pileup output."""
    input:
        pileup=rules.ft_pileup.output.pileup,
    output:
        tbi="results/pileup/{sample_id}.pileup.bed.gz.tbi",
    log:
        "logs/pileup/index_{sample_id}.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    shell:
        """
        tabix -p bed {input.pileup} 2> {log}
        """


rule merge_all_peaks:
    """Merge all peak files into a single sorted bed file.

    Each re_id is unique across all samples, so we can simply concatenate.
    Output is bgzipped and tabix indexed.
    """
    input:
        peaks=expand(
            "results/pileup/{sample_id}.peaks.bed.gz",
            sample_id=get_sample_ids_with_ft_bam(),
        ),
    output:
        bed="results/pileup/all_peaks.bed.gz",
        tbi="results/pileup/all_peaks.bed.gz.tbi",
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
            # Header from first file
            gzip -dc {input.peaks[0]} | head -1 || true

            # Data from all files (skip headers), sort by chrom and position
            for f in {input.peaks}; do
                gzip -dc "$f" | tail -n +2
            done | sort -k1,1 -k2,2n
        ) | bgzip -c > {output.bed} 2> {log}
        tabix -p bed {output.bed} 2>> {log}
        """


rule merge_peaks_with_annotations:
    """Merge peak data with annotated REs to add cluster_id.

    Joins peaks (keyed by re_id in name column) with annotated_res.bed.
    Verifies 1:1 mapping between re_id and peaks.
    """
    input:
        peaks=rules.merge_all_peaks.output.bed,
        annotated_res=rules.paf_to_graphml.output.annotated_res,
    output:
        bed="results/pileup/peaks_with_clusters.bed.gz",
        tbi="results/pileup/peaks_with_clusters.bed.gz.tbi",
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


rule paralog_accessibility_analysis:
    """Run exploratory analysis of paralog chromatin accessibility variance."""
    input:
        peaks=rules.merge_peaks_with_annotations.output.bed,
        rmd="analysis/paralog_accessibility.Rmd",
    output:
        html="results/analysis/paralog_accessibility.html",
        figures=directory("results/analysis/figures"),
    log:
        "logs/analysis/paralog_accessibility.log",
    conda:
        "../envs/r.yml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=60,
    shell:
        """
        mkdir -p results/analysis
        Rscript -e "rmarkdown::render('analysis/paralog_accessibility.Rmd', output_dir='results/analysis')" 2>&1 | tee ../{log}
        mv figures results/analysis/ 2>/dev/null || true
        """
