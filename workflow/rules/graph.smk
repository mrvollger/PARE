"""
Build the consensus-peak paralog graph.

Because RE_IDs are directly encoded in the merged filtered PAF's query/target
names and `cp:Z:`/`tp:Z:` tags, this stage no longer needs bedtools
intersections — paf_to_graph.py does both clustering (at consensus_peak_id
level) and expansion to the full (sample_id x consensus_peak_id) node set
using the re_index.
"""


rule paf_to_graphml:
    """
    Cluster consensus_peak_ids from the merged filtered PAF via Union-Find,
    then expand clusters back to every (sample_id x consensus_peak_id) node
    using the re_index. Emits GraphML, cluster TSVs, annotated PAF, and a
    node BED.
    """
    input:
        paf=rules.merge_filtered_pafs.output.paf,
        re_index=rules.extract_union_sequences.output.re_index,
    output:
        graph=output_path("graphs/re_graph.graphml"),
        re_clusters=output_path("graphs/re_clusters.tsv"),
        consensus_peak_clusters=output_path("graphs/consensus_peak_clusters.tsv"),
        annotated_paf=output_path("graphs/annotated.paf"),
        annotated_res=output_path("graphs/annotated_res.bed"),
    log:
        "logs/graph/paf_to_graph.log",
    conda:
        "../envs/python.yml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=60,
    script:
        "../scripts/paf_to_graph.py"
