#!/usr/bin/env python3
"""
Annotate RE-to-RE alignments (from all_by_all) with cluster IDs and emit a
deduped pairwise identity TSV.

Input PAF query/target names are RE_IDs: {sample_id}__{consensus_peak_id}.
Cluster assignments come from results/graphs/re_clusters.tsv
(re_id -> cluster_id) and additional metadata from the re_index.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from paf import PAFRecord  # noqa: E402


def load_re_clusters(path: str) -> dict:
    """re_clusters.tsv -> {re_id: cluster_id}"""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[col["re_id"]]] = f[col["cluster_id"]]
    return out


def load_re_index(path: str) -> dict:
    """re_index.tsv.gz -> {re_id: metadata}"""
    opener = gzip.open if path.endswith(".gz") else open
    out = {}
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[col["re_id"]]] = {
                "sample_id": f[col["sample_id"]],
                "consensus_peak_id": f[col["consensus_peak_id"]],
                "haplotype_sample": f[col["haplotype_sample"]],
                "Individual_ID": f[col["Individual_ID"]],
                "Haplotype": f[col["Haplotype"]],
            }
    return out


def main() -> None:
    sm = snakemake  # noqa: F821

    paf_in = sm.input.paf
    clusters_path = sm.input.re_clusters
    index_path = sm.input.re_index
    out_paf = sm.output.paf
    out_tsv = sm.output.tsv
    log_path = sm.log[0]

    clusters = load_re_clusters(clusters_path)
    idx = load_re_index(index_path)

    # Dedupe pairs: canonical key is tuple(sorted([q, t])). Keep the best by
    # (num_matches * blast_identity).
    best: dict = {}

    total = 0
    dropped_unknown = 0
    dropped_self = 0

    with open(paf_in) as fh:
        for line in fh:
            rec = PAFRecord.from_line(line)
            if rec is None:
                continue
            total += 1
            q, t = rec.query_name, rec.target_name
            if q == t:
                dropped_self += 1
                continue
            if q not in idx or t not in idx:
                dropped_unknown += 1
                continue
            score = rec.num_matches * (rec.blast_identity() / 100.0)
            key = tuple(sorted((q, t)))
            prev = best.get(key)
            if prev is None or score > prev[0]:
                best[key] = (score, line.rstrip("\n"), rec)

    # Write annotated PAF and TSV.
    with gzip.open(out_paf, "wt") as pfh, gzip.open(out_tsv, "wt") as tfh:
        tfh.write(
            "query_re_id\ttarget_re_id\tquery_cluster\ttarget_cluster\tsame_cluster\t"
            "query_sample\ttarget_sample\tquery_haplotype_sample\ttarget_haplotype_sample\t"
            "query_consensus_peak_id\ttarget_consensus_peak_id\t"
            "blast_identity\tgap_compressed_identity\tnum_matches\talignment_length\t"
            "query_length\ttarget_length\tstrand\n"
        )
        for (q, t), (_, line, rec) in best.items():
            qc = clusters.get(q, "unclustered")
            tc = clusters.get(t, "unclustered")
            same = qc == tc and qc != "unclustered"
            bi = rec.blast_identity() / 100.0
            gi = rec.gap_compressed_identity()
            gi_str = f"{gi/100.0:.4f}" if gi is not None else ""

            pfh.write(
                f"{line}\tqc:Z:{qc}\ttc:Z:{tc}\tsc:i:{1 if same else 0}\tbi:f:{bi:.4f}"
                + (f"\tgi:f:{gi/100.0:.4f}" if gi is not None else "")
                + "\n"
            )

            qi = idx[q]
            ti = idx[t]
            tfh.write(
                f"{q}\t{t}\t{qc}\t{tc}\t{str(same).lower()}\t"
                f"{qi['sample_id']}\t{ti['sample_id']}\t"
                f"{qi['haplotype_sample']}\t{ti['haplotype_sample']}\t"
                f"{qi['consensus_peak_id']}\t{ti['consensus_peak_id']}\t"
                f"{bi:.4f}\t{gi_str}\t"
                f"{rec.num_matches}\t{rec.alignment_block_length}\t"
                f"{rec.query_length}\t{rec.target_length}\t{rec.strand}\n"
            )

    with open(log_path, "w") as lg:
        lg.write(f"records read: {total}\n")
        lg.write(f"self-alignments dropped: {dropped_self}\n")
        lg.write(f"unknown RE_IDs dropped: {dropped_unknown}\n")
        lg.write(f"unique pairs kept: {len(best)}\n")


if __name__ == "__main__":
    main()
