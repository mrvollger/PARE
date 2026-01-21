#!/usr/bin/env python3
"""
Convert PAF alignment files to graph format for visualization.

Reads pre-computed bedtools intersections to build RE graph.

Creates an RE-centric graph where:
- Nodes: Individual RE IDs (both annotated and implicit)
- Edges: PAF alignments connecting REs (source RE -> target RE)
- Clustering: Cluster ID computed via Union-Find and stored as node attribute

Supports GraphML format for visualization in Cytoscape, Gephi, etc.
"""

import sys
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from dataclasses import dataclass
import xml.etree.ElementTree as ET
from xml.dom import minidom

try:
    from paf import PAFReader, PAFRecord, parse_contig_name
except ImportError:
    import os

    script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, script_dir)
    from paf import PAFReader, PAFRecord, parse_contig_name


class UnionFind:
    """Union-Find data structure for clustering."""

    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        root_x = self.find(x)
        root_y = self.find(y)
        if root_x == root_y:
            return
        if self.rank[root_x] < self.rank[root_y]:
            self.parent[root_x] = root_y
        elif self.rank[root_x] > self.rank[root_y]:
            self.parent[root_y] = root_x
        else:
            self.parent[root_y] = root_x
            self.rank[root_x] += 1

    def get_clusters(self) -> Dict[str, List[str]]:
        clusters = defaultdict(list)
        for x in list(self.parent.keys()):
            root = self.find(x)
            clusters[root].append(x)
        return dict(clusters)


def find_bridges_and_split(
    edges: List[Tuple[str, str, float]],
    min_cluster_size: int = 3,
    min_density_threshold: float = 0.5,
) -> Dict[str, List[str]]:
    """
    Find weakly connected components and split them at bridge edges if it improves density.

    Uses Tarjan's bridge-finding algorithm to identify edges whose removal
    disconnects the graph. Only splits if:
    1. Cluster has at least min_cluster_size nodes
    2. Original cluster density is below threshold (0.5)
    3. Both resulting sides improve density AND exceed threshold

    Args:
        edges: List of (node1, node2, weight) tuples
        min_cluster_size: Minimum cluster size to consider for splitting (default: 3)
        min_density_threshold: Density threshold - original must be below,
                              results must be at or above (default: 0.5)

    Returns:
        Dictionary mapping cluster root -> list of member node IDs
    """
    # Build adjacency list with total weight and edge counts
    adj = defaultdict(lambda: defaultdict(float))  # total weight
    edge_counts = defaultdict(lambda: defaultdict(int))  # number of edges
    nodes = set()

    for n1, n2, weight in edges:
        if n1 == n2:
            continue
        nodes.add(n1)
        nodes.add(n2)
        adj[n1][n2] += weight
        adj[n2][n1] += weight
        edge_counts[n1][n2] += 1
        edge_counts[n2][n1] += 1

    if not nodes:
        return {}

    # Tarjan's bridge-finding algorithm
    discovery = {}
    low = {}
    parent = {}
    bridges = []  # (u, v, total_weight, edge_count)
    time_counter = [0]

    def dfs(u):
        discovery[u] = low[u] = time_counter[0]
        time_counter[0] += 1

        for v in adj[u]:
            if v not in discovery:
                parent[v] = u
                dfs(v)
                low[u] = min(low[u], low[v])

                # Check if u-v is a bridge
                if low[v] > discovery[u]:
                    total_weight = adj[u][v]
                    count = edge_counts[u][v]
                    bridges.append((u, v, total_weight, count))
            elif v != parent.get(u):
                low[u] = min(low[u], discovery[v])

    # Run DFS from each unvisited node (handles disconnected components)
    for node in nodes:
        if node not in discovery:
            parent[node] = None
            dfs(node)

    print(f"  Found {len(bridges)} bridges in the graph", file=sys.stderr)

    if not bridges:
        # No bridges - use simple union-find
        uf = UnionFind()
        for n1, n2, _ in edges:
            if n1 != n2:
                uf.union(n1, n2)
            else:
                uf.find(n1)
        return uf.get_clusters()

    # For each bridge, determine if splitting improves cluster quality
    # We need to know which nodes are on each side of the bridge
    def get_component_sides(bridge_u, bridge_v):
        """Get nodes on each side of a bridge using BFS, excluding the bridge edge."""
        side_u = set()
        side_v = set()

        # BFS from u, not crossing the bridge
        queue = [bridge_u]
        side_u.add(bridge_u)
        while queue:
            node = queue.pop(0)
            for neighbor in adj[node]:
                if neighbor not in side_u:
                    # Don't cross the bridge
                    if not (node == bridge_u and neighbor == bridge_v):
                        side_u.add(neighbor)
                        queue.append(neighbor)

        # BFS from v, not crossing the bridge
        queue = [bridge_v]
        side_v.add(bridge_v)
        while queue:
            node = queue.pop(0)
            for neighbor in adj[node]:
                if neighbor not in side_v:
                    if not (node == bridge_v and neighbor == bridge_u):
                        side_v.add(neighbor)
                        queue.append(neighbor)

        return side_u, side_v

    def calculate_density(node_set):
        """
        Calculate graph density: actual_edges / possible_edges.

        A singleton (n=1) is considered fully connected (density=1.0)
        to allow splitting off loosely connected singletons.
        """
        n = len(node_set)
        if n <= 1:
            return 1.0  # Singleton is fully connected by definition
        possible_edges = n * (n - 1) / 2
        actual_edges = 0
        node_list = list(node_set)
        for i, n1 in enumerate(node_list):
            for n2 in node_list[i+1:]:
                if n2 in adj[n1]:
                    actual_edges += edge_counts[n1][n2]
        return actual_edges / possible_edges

    # Evaluate each bridge for potential splitting
    bridges_to_cut = set()
    skipped_cluster_too_small = 0
    skipped_no_density_improvement = 0

    for u, v, bridge_weight, bridge_count in bridges:
        side_u, side_v = get_component_sides(u, v)
        combined = side_u | side_v

        # Only apply to clusters with at least min_cluster_size nodes total
        if len(combined) < min_cluster_size:
            skipped_cluster_too_small += 1
            continue

        # Calculate density before and after split
        density_before = calculate_density(combined)
        density_u = calculate_density(side_u)
        density_v = calculate_density(side_v)

        # Check density improvement:
        # 1. Original cluster must have low density (below threshold)
        # 2. Both resulting sides must improve AND exceed threshold
        # This prevents splitting already well-connected clusters
        density_improves = (
            density_before < min_density_threshold
            and density_u > density_before and density_v > density_before
            and density_u >= min_density_threshold and density_v >= min_density_threshold
        )

        # Cut if density improves
        if density_improves:
            bridges_to_cut.add((min(u, v), max(u, v)))
            print(f"    Cutting bridge: sides={len(side_u)},{len(side_v)}, "
                  f"density: {density_before:.3f}->{density_u:.3f},{density_v:.3f}", file=sys.stderr)
        else:
            skipped_no_density_improvement += 1

    print(f"  Bridges skipped (cluster too small): {skipped_cluster_too_small}", file=sys.stderr)
    print(f"  Bridges where density didn't improve: {skipped_no_density_improvement}", file=sys.stderr)
    print(f"  Cutting {len(bridges_to_cut)} bridges that improve quality", file=sys.stderr)

    # Rebuild clusters excluding bridges to cut
    uf = UnionFind()
    for n1, n2, weight in edges:
        if n1 == n2:
            uf.find(n1)
            continue
        bridge_key = (min(n1, n2), max(n1, n2))
        if bridge_key not in bridges_to_cut:
            uf.union(n1, n2)
        else:
            # Ensure both nodes exist in union-find even if not merged
            uf.find(n1)
            uf.find(n2)

    # Also ensure singleton nodes are included
    for node in nodes:
        uf.find(node)

    return uf.get_clusters()


@dataclass
class RENode:
    """Represents a single RE node in the graph."""

    re_id: str
    cluster_id: str = None

    def to_dict(self) -> Dict:
        return {
            "re_id": self.re_id,
            "cluster_id": self.cluster_id or "unclustered",
        }


class REGraph:
    """Graph representation of RE alignments with clustering."""

    def __init__(self):
        self.paf_records: Dict[str, PAFRecord] = {}
        self.paf_to_res: Dict[
            str, Tuple[Optional[str], Optional[str]]
        ] = {}  # paf_id -> (source_re, target_re)
        self.nodes: Dict[str, RENode] = {}  # re_id -> RENode
        self.re_to_cluster: Dict[str, str] = {}  # re_id -> cluster_id
        self.edges: List[Dict] = []

    def load_paf_file(self, paf_file: str, skip_self: bool = True):
        print(f"Loading PAF file: {paf_file}", file=sys.stderr)
        line_num = 0
        with open(paf_file) as f:
            for line in f:
                line_num += 1
                record = PAFRecord.from_line(line)
                if not record:
                    continue
                classification = record.classification or record.classify_alignment()
                if skip_self and classification == "self":
                    continue
                paf_id = f"paf_{line_num}"
                self.paf_records[paf_id] = record
        print(f"  Loaded {len(self.paf_records)} PAF records", file=sys.stderr)

    def load_intersections(
        self,
        intersect_file: str,
        split_bridges: bool = True,
        min_cluster_size: int = 3,
    ):
        """
        Load bedtools intersection results and build RE nodes with cluster assignments.

        Format (9 columns):
        paf_chrom, paf_start, paf_end, paf_id, source_re_id,
        target_re_chrom, target_re_start, target_re_end, target_re_id

        Args:
            intersect_file: Path to bedtools intersection file
            split_bridges: If True, use bridge detection to split weakly connected clusters.
                          Default True enables quality-based splitting.
            min_cluster_size: Minimum size for clusters after splitting (default: 3)
        """
        print(f"Loading intersections: {intersect_file}", file=sys.stderr)

        # Collect all edges with weights from PAF records
        # Weight = identity * num_matches (quality-adjusted)
        all_edges = []  # List of (source_re_id, target_re_id, weight)

        with open(intersect_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                paf_id = fields[3]
                source_re_id = fields[4] if fields[4] != "." else None
                target_re_id = fields[8] if fields[8] != "." else None

                # Store RE-to-RE connection for this PAF
                self.paf_to_res[paf_id] = (source_re_id, target_re_id)

                # Calculate edge weight from PAF record if available
                weight = 1.0
                if paf_id in self.paf_records:
                    record = self.paf_records[paf_id]
                    identity = record.calculate_identity() / 100.0  # Convert to fraction
                    weight = identity * record.num_matches

                # Collect edges for clustering
                if source_re_id and target_re_id:
                    all_edges.append((source_re_id, target_re_id, weight))
                elif source_re_id:
                    # Singleton - add self-edge to ensure it's included
                    all_edges.append((source_re_id, source_re_id, 0.0))
                elif target_re_id:
                    all_edges.append((target_re_id, target_re_id, 0.0))

        print(f"  Processed {len(self.paf_to_res)} PAF-RE connections", file=sys.stderr)
        print(f"  Collected {len(all_edges)} edges for clustering", file=sys.stderr)

        # Use bridge detection if enabled, otherwise simple Union-Find
        if split_bridges:
            print(f"  Using bridge detection (min_cluster_size={min_cluster_size})", file=sys.stderr)
            uf_clusters = find_bridges_and_split(
                all_edges,
                min_cluster_size=min_cluster_size,
            )
        else:
            print("  Using simple Union-Find clustering (no bridge splitting)", file=sys.stderr)
            uf = UnionFind()
            for source, target, _ in all_edges:
                if source == target:
                    uf.find(source)  # Singleton
                else:
                    uf.union(source, target)
            uf_clusters = uf.get_clusters()

        # Filter out singleton clusters (size 1) - we're not interested in isolated REs
        singleton_count = sum(1 for re_ids in uf_clusters.values() if len(re_ids) == 1)
        uf_clusters = {k: v for k, v in uf_clusters.items() if len(v) > 1}
        print(f"  Found {len(uf_clusters)} RE clusters (filtered {singleton_count} singletons)", file=sys.stderr)

        cluster_counter = 0
        for _, re_ids in uf_clusters.items():
            cluster_id = f"cluster_{cluster_counter}"
            cluster_counter += 1

            # Create individual nodes for each RE with cluster assignment
            for re_id in re_ids:
                if re_id not in self.nodes:
                    self.nodes[re_id] = RENode(re_id=re_id, cluster_id=cluster_id)
                else:
                    self.nodes[re_id].cluster_id = cluster_id
                self.re_to_cluster[re_id] = cluster_id

        print(f"  Created {len(self.nodes)} RE nodes", file=sys.stderr)
        print(f"  Assigned to {len(uf_clusters)} clusters", file=sys.stderr)

    def build_edges(self):
        """
        Build graph edges from PAF records.

        Each PAF connects source_re_id -> target_re_id directly.
        Edges connect individual RE nodes.
        """
        print("Building graph edges...", file=sys.stderr)

        edges_created = 0
        edges_skipped_no_source = 0
        edges_skipped_no_target = 0

        for paf_id, record in self.paf_records.items():
            source_re_id, target_re_id = self.paf_to_res.get(paf_id, (None, None))

            if not source_re_id:
                edges_skipped_no_source += 1
                continue

            if not target_re_id:
                edges_skipped_no_target += 1
                continue

            # Ensure both RE nodes exist
            if source_re_id not in self.nodes:
                self.nodes[source_re_id] = RENode(re_id=source_re_id)

            if target_re_id not in self.nodes:
                self.nodes[target_re_id] = RENode(re_id=target_re_id)

            # Create edge connecting individual REs
            classification = record.classification or record.classify_alignment()

            # Weight for ForceAtlas2: positive attracts, negative repels
            # Paralogs cluster tightly, orthologs/allelic actively repel
            weight_map = {
                "self": 10,
                "paralog": 10,
                "allelic": 1,
                "ortholog": 0.1,
                "unknown": 0.1,
            }
            weight = weight_map.get(classification, 0.1)

            edge = {
                "source": source_re_id,
                "target": target_re_id,
                "paf_id": paf_id,
                "classification": classification,
                "weight": weight,
                "identity": record.calculate_identity(),
                "coverage": record.calculate_re_coverage(),
                "overlap": record.calculate_self_overlap(),
                "query_name": record.query_name,
                "query_start": record.query_start,
                "query_end": record.query_end,
                "target_name": record.target_name,
                "target_start": record.target_start,
                "target_end": record.target_end,
                "strand": record.strand,
                "num_matches": record.num_matches,
                "alignment_length": record.alignment_block_length,
                "mapq": record.mapping_quality,
                "query_sample": record.query_sample,
                "query_haplotype": record.query_haplotype,
                "target_sample": record.target_sample,
                "target_haplotype": record.target_haplotype,
            }

            for tag_name, tag_value in record.tags.items():
                if tag_name not in edge:
                    edge[tag_name] = tag_value

            self.edges.append(edge)
            edges_created += 1

        print(f"  Created {edges_created} edges", file=sys.stderr)
        print(
            f"  Skipped {edges_skipped_no_source} PAFs with no source RE",
            file=sys.stderr,
        )
        print(
            f"  Skipped {edges_skipped_no_target} PAFs with no target RE overlap",
            file=sys.stderr,
        )

    def get_stats(self) -> Dict:
        classification_counts = defaultdict(int)
        for edge in self.edges:
            classification_counts[edge["classification"]] += 1

        # Count cluster sizes
        cluster_sizes = defaultdict(int)
        num_clustered = 0
        for node in self.nodes.values():
            if node.cluster_id:
                cluster_sizes[node.cluster_id] += 1
                num_clustered += 1

        return {
            "num_nodes": len(self.nodes),
            "num_edges": len(self.edges),
            "num_paf_records": len(self.paf_records),
            "classifications": dict(classification_counts),
            "cluster_stats": {
                "num_clusters": len(cluster_sizes),
                "clustered_nodes": num_clustered,
                "unclustered_nodes": len(self.nodes) - num_clustered,
                "avg_cluster_size": sum(cluster_sizes.values()) / len(cluster_sizes)
                if cluster_sizes
                else 0,
                "max_cluster_size": max(cluster_sizes.values()) if cluster_sizes else 0,
            },
        }

    def to_graphml(self) -> str:
        root = ET.Element(
            "graphml",
            {
                "xmlns": "http://graphml.graphdrawing.org/xmlns",
                "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
                "xsi:schemaLocation": "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd",
            },
        )

        # Node attributes: re_id and cluster_id
        node_attrs = [
            ("re_id", "string"),
            ("cluster_id", "string"),
        ]
        for i, (name, attr_type) in enumerate(node_attrs):
            ET.SubElement(
                root,
                "key",
                {
                    "id": f"n{i}",
                    "for": "node",
                    "attr.name": name,
                    "attr.type": attr_type,
                },
            )

        # Edge attributes (including visual color attribute)
        edge_attrs = [
            ("paf_id", "string"),
            ("classification", "string"),
            ("weight", "int"),
            ("identity", "double"),
            ("coverage", "double"),
            ("overlap", "int"),
            ("query_name", "string"),
            ("query_start", "int"),
            ("query_end", "int"),
            ("target_name", "string"),
            ("target_start", "int"),
            ("target_end", "int"),
            ("strand", "string"),
            ("num_matches", "int"),
            ("alignment_length", "int"),
            ("mapq", "int"),
            ("query_sample", "string"),
            ("query_haplotype", "string"),
            ("target_sample", "string"),
            ("target_haplotype", "string"),
        ]
        edge_attr_map = {}
        for i, (name, attr_type) in enumerate(edge_attrs):
            ET.SubElement(
                root,
                "key",
                {
                    "id": f"e{i}",
                    "for": "edge",
                    "attr.name": name,
                    "attr.type": attr_type,
                },
            )
            edge_attr_map[name] = f"e{i}"

        # Add visual color attribute for Gephi/Cytoscape
        color_key_id = f"e{len(edge_attrs)}"
        ET.SubElement(
            root,
            "key",
            {
                "id": color_key_id,
                "for": "edge",
                "attr.name": "color",
                "attr.type": "string",
            },
        )

        graph = ET.SubElement(
            root, "graph", {"id": "RE_Graph", "edgedefault": "directed"}
        )

        # Add individual RE nodes
        for node in self.nodes.values():
            node_elem = ET.SubElement(graph, "node", {"id": node.re_id})
            node_dict = node.to_dict()
            for i, (attr_name, _) in enumerate(node_attrs):
                if attr_name in node_dict:
                    data = ET.SubElement(node_elem, "data", {"key": f"n{i}"})
                    data.text = str(node_dict[attr_name])

        # Color map for classifications (RGB hex colors)
        classification_colors = {
            "self": "#808080",  # Gray
            "paralog": "#E74C3C",  # Red
            "allelic": "#3498DB",  # Blue
            "ortholog": "#2ECC71",  # Green
            "unknown": "#95A5A6",  # Light gray
        }

        # Add edges connecting REs
        for i, edge in enumerate(self.edges):
            edge_elem = ET.SubElement(
                graph,
                "edge",
                {"id": f"e{i}", "source": edge["source"], "target": edge["target"]},
            )
            for attr_name, key_id in edge_attr_map.items():
                if attr_name in edge:
                    data = ET.SubElement(edge_elem, "data", {"key": key_id})
                    data.text = str(edge[attr_name])

            # Add color based on classification
            classification = edge.get("classification", "unknown")
            color = classification_colors.get(classification, "#95A5A6")
            color_data = ET.SubElement(edge_elem, "data", {"key": color_key_id})
            color_data.text = color

        xml_str = ET.tostring(root, encoding="unicode")
        dom = minidom.parseString(xml_str)
        return dom.toprettyxml(indent="  ")


def snakemake_main():
    paf_file = snakemake.input.paf
    intersect_file = snakemake.input.intersections
    annotated_bed = snakemake.input.get("annotated_bed")
    implicit_bed = snakemake.input.get("implicit_bed")

    output_graph = snakemake.output.graph
    output_cluster_map = snakemake.output.cluster_map
    output_annotated_paf = snakemake.output.annotated_paf
    output_annotated_res = snakemake.output.annotated_res

    log_file = snakemake.log[0] if snakemake.log else None
    skip_self = snakemake.params.get("skip_self", True)

    if log_file:
        log = open(log_file, "w")
        old_stderr = sys.stderr
        sys.stderr = log

    try:
        graph = REGraph()
        graph.load_paf_file(paf_file, skip_self=skip_self)
        graph.load_intersections(intersect_file)
        graph.build_edges()

        stats = graph.get_stats()
        print(f"\nGraph Statistics:")
        print(f"  Nodes (individual REs): {stats['num_nodes']}")
        print(f"  Edges: {stats['num_edges']}")
        print(f"  PAF records loaded: {stats['num_paf_records']}")
        print(f"  Cluster Statistics:")
        print(f"    Total clusters: {stats['cluster_stats']['num_clusters']}")
        print(f"    Clustered nodes: {stats['cluster_stats']['clustered_nodes']}")
        print(f"    Unclustered nodes: {stats['cluster_stats']['unclustered_nodes']}")
        print(f"    Avg cluster size: {stats['cluster_stats']['avg_cluster_size']:.2f}")
        print(f"    Max cluster size: {stats['cluster_stats']['max_cluster_size']}")
        print(f"  Edge Classifications:")
        for classification, count in sorted(stats["classifications"].items()):
            print(f"    {classification}: {count}")

        print(f"\nExporting to GraphML format...")
        output = graph.to_graphml()
        with open(output_graph, "w") as f:
            f.write(output)
        print(f"Wrote graph to: {output_graph}")

        print(f"\nExporting cluster map...")
        with open(output_cluster_map, "w") as f:
            f.write("re_id\tcluster_id\n")
            for re_id, node in sorted(graph.nodes.items()):
                cluster_id = node.cluster_id or "unclustered"
                f.write(f"{re_id}\t{cluster_id}\n")
        print(f"Wrote cluster map to: {output_cluster_map}")

        print(f"\nAnnotating PAF file...")
        with open(output_annotated_paf, "w") as f_out:
            for paf_id, record in sorted(graph.paf_records.items()):
                source_re_id, target_re_id = graph.paf_to_res.get(paf_id, (None, None))
                cluster_id = (
                    graph.nodes.get(source_re_id).cluster_id
                    if source_re_id in graph.nodes
                    else "."
                )

                line = record.to_paf_line()
                line += f"\tpi:Z:{paf_id}\tci:Z:{cluster_id or '.'}"
                f_out.write(line + "\n")
        print(f"Wrote annotated PAF to: {output_annotated_paf}")

        print("\nAnnotating RE BED file...")
        with open(output_annotated_res, "w") as f_out:
            for bed_file in [annotated_bed, implicit_bed]:
                if bed_file:
                    with open(bed_file) as f_in:
                        for line in f_in:
                            if line.startswith("#"):
                                continue
                            fields = line.strip().split("\t")
                            if len(fields) >= 4:
                                chrom, start, end, re_id = fields[:4]
                                cluster_id = (
                                    graph.nodes.get(re_id).cluster_id
                                    if re_id in graph.nodes
                                    else "."
                                )
                                # Parse sample and haplotype from chrom name
                                parsed = parse_contig_name(chrom)
                                sample = parsed["sample"]
                                haplotype = parsed["haplotype"]
                                f_out.write(
                                    f"{chrom}\t{start}\t{end}\t{re_id}\t{cluster_id or '.'}\t{sample}\t{haplotype}\n"
                                )
        print(f"Wrote annotated REs to: {output_annotated_res}")

    finally:
        if log_file:
            sys.stderr = old_stderr
            log.close()


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_main()
    else:
        print("This script is designed to run within Snakemake.", file=sys.stderr)
        sys.exit(1)
