"""
Rules for converting PAF alignments to graph formats for visualization
"""


rule merge_and_sort_pafs:
    """Merge all filtered PAF files and sort by query coordinates for stable line numbers."""
    input:
        pafs=expand(
            rules.filter_alignments.output.paf,
            sample_id=get_all_sample_ids(),
        ),
    output:
        merged=temp("temp/graph/merged_sorted.paf"),
    log:
        "logs/graph/merge_sort_paf.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=30,
    shell:
        """
        cat {input.pafs} | sort -k1,1 -k3,3n -k4,4n > {output.merged} 2> {log}
        """


rule paf_to_bed:
    """Convert merged PAF to sorted BED with source RE ID from id:Z: tag."""
    input:
        paf=rules.merge_and_sort_pafs.output.merged,
    output:
        bed=temp("temp/graph/paf_positions.bed"),
    log:
        "logs/graph/paf_to_bed.log",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    shell:
        """
        awk 'BEGIN{{OFS="\\t"}} {{
            re_id = ".";
            for(i=13; i<=NF; i++) {{
                if($i ~ /^id:Z:/) {{
                    split($i, a, ":");
                    re_id = a[3];
                    break;
                }}
            }}
            print $1, $3, $4, "paf_"NR, re_id
        }}' {input.paf} | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """


rule create_implicit_res:
    """Find PAF positions with no annotated RE overlap, merge them, and create implicit RE IDs (IMPLICIT_chrom_start_end)."""
    input:
        paf_bed=rules.paf_to_bed.output.bed,
        re_bed=rules.merge_unslopped_beds.output.bed,
    output:
        implicit_bed=temp("temp/graph/implicit_res.bed"),
    log:
        "logs/graph/create_implicit_res.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    params:
        reciprocal_overlap=config["clustering"]["reciprocal_overlap"],
    shell:
        """
        bedtools intersect \
            -a {input.paf_bed} \
            -b {input.re_bed} \
            -v \
            -f {params.reciprocal_overlap} \
            -r \
            | cut -f1-3 \
            | bedtools sort \
            | bedtools merge \
            | awk 'BEGIN{{OFS="\\t"}} {{
                implicit_id = "IMPLICIT_" $1 "_" $2 "_" $3;
                print $1, $2, $3, implicit_id
            }}' \
            > {output.implicit_bed} 2> {log}
        """


rule intersect_paf_with_res:
    """Intersect PAF positions with annotated REs (reciprocal) and implicit REs (non-reciprocal), then combine."""
    input:
        paf_bed=rules.paf_to_bed.output.bed,
        annotated_bed=rules.merge_unslopped_beds.output.bed,
        implicit_bed=rules.create_implicit_res.output.implicit_bed,
    output:
        intersect=temp("temp/graph/paf_re_intersections.tsv"),
        annotated=temp("temp/graph/paf_re_intersections.annotated.tsv"),
        implicit=temp("temp/graph/paf_re_intersections.implicit.tsv"),
    log:
        "logs/graph/intersect.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    params:
        reciprocal_overlap=config["clustering"]["reciprocal_overlap"],
    shell:
        """
        bedtools intersect \
            -a {input.paf_bed} \
            -b {input.annotated_bed} \
            -loj \
            -f {params.reciprocal_overlap} \
            -r \
            | tee >(awk '$9 != "."' > {output.annotated}) \
            | awk '$9 == "."' \
            | cut -f1-5 \
            | bedtools intersect \
                -a stdin \
                -b {input.implicit_bed} \
                -wa -wb \
                -f {params.reciprocal_overlap} \
            > {output.implicit} 2> {log}

        cat {output.annotated} {output.implicit} \
            | sort -k1,1 -k2,2n -k3,3n \
            > {output.intersect}
        """


rule paf_to_graphml:
    """Convert PAF alignments to GraphML with individual RE nodes, direct edges, and cluster_id annotations."""
    input:
        paf=rules.merge_and_sort_pafs.output.merged,
        intersections=rules.intersect_paf_with_res.output.intersect,
        annotated_bed=rules.merge_unslopped_beds.output.bed,
        implicit_bed=rules.create_implicit_res.output.implicit_bed,
    output:
        graph="results/graphs/re_graph.graphml",
        cluster_map="results/graphs/re_clusters.tsv",
        annotated_paf="results/graphs/annotated.paf",
        annotated_res="results/graphs/annotated_res.bed",
    log:
        "logs/graph/re_graph.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=60,
    params:
        skip_self=False,
        format="graphml",
    script:
        "../scripts/paf_to_graph.py"
