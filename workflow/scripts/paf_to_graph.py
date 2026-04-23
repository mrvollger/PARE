#!/usr/bin/env python3
"""
Build the PARE RE graph from a merged filtered PAF.

Input PAF:
  - query/target names are union-BED RE_IDs: {sample_id}__{consensus_peak_id}
  - tags cp:Z:<query consensus_peak_id> and tp:Z:<target consensus_peak_id>
  - ct:Z:<classification>

Clustering happens at the consensus_peak_id level (paralog edges). Every
(sample_id, consensus_peak_id) node inherits its cluster from its
consensus_peak_id. The re_index.tsv.gz supplies the full set of nodes
(including non-primary sample_ids).

Outputs:
  results/graphs/re_graph.graphml            full graph (sample_id x consensus_peak_id nodes)
  results/graphs/re_clusters.tsv             re_id -> cluster_id
  results/graphs/consensus_peak_clusters.tsv consensus_peak_id -> cluster_id (collapsed view)
  results/graphs/annotated.paf               input PAF with ci:Z:<cluster_id> appended
  results/graphs/annotated_res.bed           BED of every RE node with cluster + metadata
"""

from __future__ import annotations

import gzip
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from xml.dom import minidom

sys.path.insert(0, str(Path(__file__).resolve().parent))

from paf import PAFRecord  # noqa: E402


# --- Union-Find ---------------------------------------------------------------


class UnionFind:
    def __init__(self):
        self.parent: Dict[str, str] = {}
        self.rank: Dict[str, int] = {}

    def find(self, x: str) -> str:
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        # iterative path compression
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, x: str, y: str) -> None:
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

    def clusters(self) -> Dict[str, List[str]]:
        out: Dict[str, List[str]] = defaultdict(list)
        for x in list(self.parent):
            out[self.find(x)].append(x)
        return dict(out)


# --- Loading ------------------------------------------------------------------


def _tag(rec: PAFRecord, name: str, default=None):
    return rec.tags.get(name, default)


def load_paf_edges(
    paf_path: str,
) -> Tuple[List[PAFRecord], List[Tuple[str, str]]]:
    """
    Return (records, consensus_peak_edges). Each edge is (src_peak, tgt_peak).
    Self-loops at the consensus_peak level are skipped.
    """
    records: List[PAFRecord] = []
    peak_edges: List[Tuple[str, str]] = []
    with open(paf_path) as f:
        for line in f:
            rec = PAFRecord.from_line(line)
            if rec is None:
                continue
            records.append(rec)
            src = _tag(rec, "cp")
            tgt = _tag(rec, "tp")
            if src is None or tgt is None or src == tgt:
                continue
            peak_edges.append((src, tgt))
    return records, peak_edges


def cluster_consensus_peaks(
    edges: Iterable[Tuple[str, str]],
) -> Dict[str, str]:
    """
    Union-Find cluster assignment for consensus_peak_ids.
    Returns {consensus_peak_id -> cluster_id}.
    """
    uf = UnionFind()
    for a, b in edges:
        uf.union(a, b)
    clusters = uf.clusters()

    # Drop singletons — they aren't interesting as a cluster.
    peak_to_cluster: Dict[str, str] = {}
    i = 0
    for members in clusters.values():
        if len(members) < 2:
            continue
        cid = f"cluster_{i}"
        for m in members:
            peak_to_cluster[m] = cid
        i += 1
    return peak_to_cluster


# --- RE index -----------------------------------------------------------------


def load_re_index(path: str) -> List[dict]:
    """Read re_index.tsv.gz into a list of dicts (order preserved)."""
    opener = gzip.open if path.endswith(".gz") else open
    rows: List[dict] = []
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, values)))
    return rows


# --- GraphML writer -----------------------------------------------------------


CLASS_COLORS = {
    "self": "#808080",
    "paralog": "#E74C3C",
    "allelic": "#3498DB",
    "allelic-paralog": "#9B59B6",
    "ortholog": "#2ECC71",
    "ortholog-paralog": "#16A085",
    "unknown": "#95A5A6",
}


def build_graphml(
    re_index: List[dict],
    re_to_cluster: Dict[str, str],
    records: List[PAFRecord],
    paf_to_cluster: Dict[str, str],
) -> str:
    root = ET.Element(
        "graphml",
        {
            "xmlns": "http://graphml.graphdrawing.org/xmlns",
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:schemaLocation": "http://graphml.graphdrawing.org/xmlns "
            "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd",
        },
    )

    node_attrs = [
        ("re_id", "string"),
        ("consensus_peak_id", "string"),
        ("sample_id", "string"),
        ("haplotype_sample", "string"),
        ("Individual_ID", "string"),
        ("Haplotype", "string"),
        ("cluster_id", "string"),
        ("is_peak", "string"),
        ("is_primary_sample", "string"),
    ]
    for i, (name, typ) in enumerate(node_attrs):
        ET.SubElement(
            root,
            "key",
            {"id": f"n{i}", "for": "node", "attr.name": name, "attr.type": typ},
        )

    edge_attrs = [
        ("classification", "string"),
        ("identity", "double"),
        ("query_coverage", "double"),
        ("target_coverage", "double"),
        ("cluster_id", "string"),
        ("color", "string"),
    ]
    for i, (name, typ) in enumerate(edge_attrs):
        ET.SubElement(
            root,
            "key",
            {"id": f"e{i}", "for": "edge", "attr.name": name, "attr.type": typ},
        )

    graph = ET.SubElement(
        root, "graph", {"id": "RE_Graph", "edgedefault": "undirected"}
    )

    # Nodes: every sample_id × consensus_peak_id from the re_index
    for row in re_index:
        cluster = re_to_cluster.get(row["re_id"], "unclustered")
        node = ET.SubElement(graph, "node", {"id": row["re_id"]})
        values = {
            "re_id": row["re_id"],
            "consensus_peak_id": row["consensus_peak_id"],
            "sample_id": row["sample_id"],
            "haplotype_sample": row["haplotype_sample"],
            "Individual_ID": row["Individual_ID"],
            "Haplotype": row["Haplotype"],
            "cluster_id": cluster,
            "is_peak": row.get("is_peak", ""),
            "is_primary_sample": row.get("is_primary_sample", ""),
        }
        for i, (name, _) in enumerate(node_attrs):
            d = ET.SubElement(node, "data", {"key": f"n{i}"})
            d.text = values.get(name, "")

    # Edges: one per paralog PAF record, at the RE_ID level
    for idx, rec in enumerate(records):
        klass = _tag(rec, "ct") or rec.classification or "unknown"
        identity = float(_tag(rec, "bi") or rec.blast_identity())
        qcov = float(_tag(rec, "qc") or 0.0)
        tcov = float(_tag(rec, "tc") or 0.0)
        cid = paf_to_cluster.get(id(rec), "unclustered")
        edge = ET.SubElement(
            graph,
            "edge",
            {"id": f"e{idx}", "source": rec.query_name, "target": rec.target_name},
        )
        values = {
            "classification": klass,
            "identity": f"{identity:.4f}",
            "query_coverage": f"{qcov:.4f}",
            "target_coverage": f"{tcov:.4f}",
            "cluster_id": cid,
            "color": CLASS_COLORS.get(klass, CLASS_COLORS["unknown"]),
        }
        for i, (name, _) in enumerate(edge_attrs):
            d = ET.SubElement(edge, "data", {"key": f"e{i}"})
            d.text = str(values[name])

    return minidom.parseString(ET.tostring(root, encoding="unicode")).toprettyxml(
        indent="  "
    )


# --- Main ---------------------------------------------------------------------


def main() -> None:
    sm = snakemake  # noqa: F821

    paf_path = sm.input.paf
    re_index_path = sm.input.re_index

    out_graph = sm.output.graph
    out_re_clusters = sm.output.re_clusters
    out_peak_clusters = sm.output.consensus_peak_clusters
    out_annotated_paf = sm.output.annotated_paf
    out_annotated_res = sm.output.annotated_res
    log_path = sm.log[0]

    with open(log_path, "w") as log:
        def say(m: str) -> None:
            print(m, file=log, flush=True)

        say(f"loading PAF: {paf_path}")
        records, peak_edges = load_paf_edges(paf_path)
        say(f"  records: {len(records)}, consensus-peak edges: {len(peak_edges)}")

        say("clustering consensus peaks")
        peak_to_cluster = cluster_consensus_peaks(peak_edges)
        num_clusters = len({c for c in peak_to_cluster.values()})
        say(f"  clusters: {num_clusters}  peaks in clusters: {len(peak_to_cluster)}")

        say(f"loading re_index: {re_index_path}")
        re_index = load_re_index(re_index_path)
        say(f"  rows: {len(re_index)}")

        # Expand consensus_peak clusters to every (sample_id × consensus_peak_id)
        re_to_cluster = {
            row["re_id"]: peak_to_cluster.get(row["consensus_peak_id"], "unclustered")
            for row in re_index
        }

        # Map each PAF record to its cluster (by source consensus_peak_id)
        paf_to_cluster: Dict[int, str] = {}
        for rec in records:
            src = _tag(rec, "cp")
            paf_to_cluster[id(rec)] = peak_to_cluster.get(src, "unclustered")

        say(f"writing re_clusters: {out_re_clusters}")
        with open(out_re_clusters, "w") as fh:
            fh.write("re_id\tconsensus_peak_id\tsample_id\tcluster_id\n")
            for row in re_index:
                cid = re_to_cluster[row["re_id"]]
                fh.write(
                    f"{row['re_id']}\t{row['consensus_peak_id']}\t{row['sample_id']}\t{cid}\n"
                )

        say(f"writing consensus_peak_clusters: {out_peak_clusters}")
        with open(out_peak_clusters, "w") as fh:
            fh.write("consensus_peak_id\tcluster_id\n")
            for peak, cid in sorted(peak_to_cluster.items()):
                fh.write(f"{peak}\t{cid}\n")

        say(f"writing annotated PAF: {out_annotated_paf}")
        with open(out_annotated_paf, "w") as fh:
            for rec in records:
                cid = paf_to_cluster[id(rec)]
                line = rec.to_paf_line()
                fh.write(f"{line}\tci:Z:{cid}\n")

        say(f"writing annotated_res BED: {out_annotated_res}")
        with open(out_annotated_res, "w") as fh:
            fh.write(
                "#asm_chr\tasm_start\tasm_end\tre_id\tcluster_id\tconsensus_peak_id\t"
                "sample_id\thaplotype_sample\tIndividual_ID\tHaplotype\tis_peak\tis_primary_sample\n"
            )
            for row in re_index:
                if not row.get("asm_chr") or not row.get("asm_start") or not row.get("asm_end"):
                    continue
                cid = re_to_cluster[row["re_id"]]
                fh.write(
                    f"{row['asm_chr']}\t{row['asm_start']}\t{row['asm_end']}\t"
                    f"{row['re_id']}\t{cid}\t{row['consensus_peak_id']}\t"
                    f"{row['sample_id']}\t{row['haplotype_sample']}\t"
                    f"{row['Individual_ID']}\t{row['Haplotype']}\t"
                    f"{row.get('is_peak', '')}\t{row.get('is_primary_sample', '')}\n"
                )

        say(f"writing GraphML: {out_graph}")
        graphml = build_graphml(re_index, re_to_cluster, records, paf_to_cluster)
        with open(out_graph, "w") as fh:
            fh.write(graphml)

        say("done")


if __name__ == "__main__":
    main()
