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
    from paf import PAFReader, PAFRecord
except ImportError:
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, script_dir)
    from paf import PAFReader, PAFRecord


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


@dataclass
class RENode:
    """Represents a single RE node in the graph."""
    re_id: str
    cluster_id: str = None

    def to_dict(self) -> Dict:
        return {
            're_id': self.re_id,
            'cluster_id': self.cluster_id or 'unclustered',
        }


class REGraph:
    """Graph representation of RE alignments with clustering."""

    def __init__(self):
        self.paf_records: Dict[str, PAFRecord] = {}
        self.paf_to_res: Dict[str, Tuple[Optional[str], Optional[str]]] = {}  # paf_id -> (source_re, target_re)
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
                if skip_self and classification == 'self':
                    continue
                paf_id = f"paf_{line_num}"
                self.paf_records[paf_id] = record
        print(f"  Loaded {len(self.paf_records)} PAF records", file=sys.stderr)

    def load_intersections(self, intersect_file: str):
        """
        Load bedtools intersection results and build RE nodes with cluster assignments.

        Format (9 columns):
        paf_chrom, paf_start, paf_end, paf_id, source_re_id,
        target_re_chrom, target_re_start, target_re_end, target_re_id
        """
        print(f"Loading intersections: {intersect_file}", file=sys.stderr)
        uf = UnionFind()

        with open(intersect_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                paf_id = fields[3]
                source_re_id = fields[4] if fields[4] != '.' else None
                target_re_id = fields[8] if fields[8] != '.' else None

                # Store RE-to-RE connection for this PAF
                self.paf_to_res[paf_id] = (source_re_id, target_re_id)

                # Cluster REs that overlap
                if source_re_id and target_re_id:
                    # Both REs exist - union them
                    uf.union(source_re_id, target_re_id)
                elif source_re_id:
                    # Only source RE exists - ensure it's in the union-find
                    uf.find(source_re_id)
                elif target_re_id:
                    # Only target RE exists (shouldn't happen often)
                    uf.find(target_re_id)

        print(f"  Processed {len(self.paf_to_res)} PAF-RE connections", file=sys.stderr)

        # Build RE nodes with cluster assignments from Union-Find
        uf_clusters = uf.get_clusters()
        print(f"  Found {len(uf_clusters)} RE clusters", file=sys.stderr)

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

            # Weight hierarchy: self(4) > paralog(3) > allelic(2) > ortholog(1)
            weight_map = {'self': 4, 'paralog': 3, 'allelic': 2, 'ortholog': 1, 'unknown': 0}
            weight = weight_map.get(classification, 0)

            edge = {
                'source': source_re_id,
                'target': target_re_id,
                'paf_id': paf_id,
                'classification': classification,
                'weight': weight,
                'identity': record.calculate_identity(),
                'coverage': record.calculate_re_coverage(),
                'overlap': record.calculate_self_overlap(),
                'query_name': record.query_name,
                'query_start': record.query_start,
                'query_end': record.query_end,
                'target_name': record.target_name,
                'target_start': record.target_start,
                'target_end': record.target_end,
                'strand': record.strand,
                'num_matches': record.num_matches,
                'alignment_length': record.alignment_block_length,
                'mapq': record.mapping_quality,
                'query_sample': record.query_sample,
                'query_haplotype': record.query_haplotype,
                'target_sample': record.target_sample,
                'target_haplotype': record.target_haplotype,
            }

            for tag_name, tag_value in record.tags.items():
                if tag_name not in edge:
                    edge[tag_name] = tag_value

            self.edges.append(edge)
            edges_created += 1

        print(f"  Created {edges_created} edges", file=sys.stderr)
        print(f"  Skipped {edges_skipped_no_source} PAFs with no source RE", file=sys.stderr)
        print(f"  Skipped {edges_skipped_no_target} PAFs with no target RE overlap", file=sys.stderr)

    def get_stats(self) -> Dict:
        classification_counts = defaultdict(int)
        for edge in self.edges:
            classification_counts[edge['classification']] += 1

        # Count cluster sizes
        cluster_sizes = defaultdict(int)
        num_clustered = 0
        for node in self.nodes.values():
            if node.cluster_id:
                cluster_sizes[node.cluster_id] += 1
                num_clustered += 1

        return {
            'num_nodes': len(self.nodes),
            'num_edges': len(self.edges),
            'num_paf_records': len(self.paf_records),
            'classifications': dict(classification_counts),
            'cluster_stats': {
                'num_clusters': len(cluster_sizes),
                'clustered_nodes': num_clustered,
                'unclustered_nodes': len(self.nodes) - num_clustered,
                'avg_cluster_size': sum(cluster_sizes.values()) / len(cluster_sizes) if cluster_sizes else 0,
                'max_cluster_size': max(cluster_sizes.values()) if cluster_sizes else 0,
            },
        }

    def to_graphml(self) -> str:
        root = ET.Element('graphml', {
            'xmlns': 'http://graphml.graphdrawing.org/xmlns',
            'xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance',
            'xsi:schemaLocation': 'http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd'
        })

        # Node attributes: re_id and cluster_id
        node_attrs = [
            ('re_id', 'string'),
            ('cluster_id', 'string'),
        ]
        for i, (name, attr_type) in enumerate(node_attrs):
            ET.SubElement(root, 'key', {
                'id': f'n{i}', 'for': 'node',
                'attr.name': name, 'attr.type': attr_type
            })

        # Edge attributes
        edge_attrs = [
            ('paf_id', 'string'), ('classification', 'string'), ('weight', 'int'),
            ('identity', 'double'), ('coverage', 'double'), ('overlap', 'int'),
            ('query_name', 'string'), ('query_start', 'int'), ('query_end', 'int'),
            ('target_name', 'string'), ('target_start', 'int'), ('target_end', 'int'),
            ('strand', 'string'), ('num_matches', 'int'), ('alignment_length', 'int'),
            ('mapq', 'int'), ('query_sample', 'string'), ('query_haplotype', 'string'),
            ('target_sample', 'string'), ('target_haplotype', 'string'),
        ]
        edge_attr_map = {}
        for i, (name, attr_type) in enumerate(edge_attrs):
            ET.SubElement(root, 'key', {
                'id': f'e{i}', 'for': 'edge',
                'attr.name': name, 'attr.type': attr_type
            })
            edge_attr_map[name] = f'e{i}'

        graph = ET.SubElement(root, 'graph', {'id': 'RE_Graph', 'edgedefault': 'directed'})

        # Add individual RE nodes
        for node in self.nodes.values():
            node_elem = ET.SubElement(graph, 'node', {'id': node.re_id})
            node_dict = node.to_dict()
            for i, (attr_name, _) in enumerate(node_attrs):
                if attr_name in node_dict:
                    data = ET.SubElement(node_elem, 'data', {'key': f'n{i}'})
                    data.text = str(node_dict[attr_name])

        # Add edges connecting REs
        for i, edge in enumerate(self.edges):
            edge_elem = ET.SubElement(graph, 'edge', {
                'id': f'e{i}', 'source': edge['source'], 'target': edge['target']
            })
            for attr_name, key_id in edge_attr_map.items():
                if attr_name in edge:
                    data = ET.SubElement(edge_elem, 'data', {'key': key_id})
                    data.text = str(edge[attr_name])

        xml_str = ET.tostring(root, encoding='unicode')
        dom = minidom.parseString(xml_str)
        return dom.toprettyxml(indent='  ')


def snakemake_main():
    paf_file = snakemake.input.paf
    intersect_file = snakemake.input.intersections
    output_file = snakemake.output.graph
    log_file = snakemake.log[0] if snakemake.log else None
    skip_self = snakemake.params.get('skip_self', True)
    
    if log_file:
        log = open(log_file, 'w')
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
        for classification, count in sorted(stats['classifications'].items()):
            print(f"    {classification}: {count}")
        
        print(f"\nExporting to GraphML format...")
        output = graph.to_graphml()
        
        with open(output_file, 'w') as f:
            f.write(output)
        
        print(f"Wrote graph to: {output_file}")
    
    finally:
        if log_file:
            sys.stderr = old_stderr
            log.close()


if __name__ == "__main__":
    if 'snakemake' in globals():
        snakemake_main()
    else:
        print("This script is designed to run within Snakemake.", file=sys.stderr)
        sys.exit(1)

