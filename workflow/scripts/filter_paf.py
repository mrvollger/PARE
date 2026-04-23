#!/usr/bin/env python3
"""
Filter a per-haplotype_sample PAF by identity and peak-coverage thresholds,
and annotate each record with metadata tags.

Alignment input format:
  - Query and target names are both union-BED RE_IDs: {sample_id}__{consensus_peak_id}
  - Coordinates are in slop-local space (0..slop_seq length)

Peak coverage is computed against the *inner* peak (unslopped) region within
each slopped sequence, using the re_index to look up asm_start/asm_slop_start
per RE_ID. An alignment must cover at least `min_coverage` percent of both the
query peak AND the target peak.

Each passing record gets these tags:
  id:Z:  source RE_ID (query)  -- kept for backward compat with graph scripts
  cp:Z:  query consensus_peak_id
  tp:Z:  target consensus_peak_id
  qs:Z:  query sample_id
  ts:Z:  target sample_id
  qh:Z:  query haplotype_sample
  th:Z:  target haplotype_sample
  ct:Z:  classification (paralog / allelic-paralog / ortholog-paralog / self / unknown)
  bi:f:  BLAST identity
  qc:f:  query peak coverage %
  tc:f:  target peak coverage %
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from paf import PAFReader  # noqa: E402


def _to_int(v):
    if v == "" or v is None:
        return None
    return int(v)


def _to_bool(v):
    if v is None or v == "":
        return None
    return v in ("True", "true", "TRUE", "1", True)


def load_re_index(path: str) -> dict:
    """Load re_index.tsv.gz into {re_id: metadata dict}."""
    opener = gzip.open if path.endswith(".gz") else open
    idx = {}
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        required = {
            "re_id",
            "sample_id",
            "consensus_peak_id",
            "haplotype_sample",
            "Individual_ID",
            "Haplotype",
            "asm_start",
            "asm_end",
            "asm_slop_start",
            "asm_slop_end",
            "is_peak",
            "is_SD",
        }
        missing = required - set(col)
        if missing:
            raise ValueError(f"re_index missing columns: {missing}")

        for line in fh:
            f = line.rstrip("\n").split("\t")
            re_id = f[col["re_id"]]
            idx[re_id] = {
                "sample_id": f[col["sample_id"]],
                "consensus_peak_id": f[col["consensus_peak_id"]],
                "haplotype_sample": f[col["haplotype_sample"]],
                "Individual_ID": f[col["Individual_ID"]],
                "Haplotype": f[col["Haplotype"]],
                "asm_start": _to_int(f[col["asm_start"]]),
                "asm_end": _to_int(f[col["asm_end"]]),
                "asm_slop_start": _to_int(f[col["asm_slop_start"]]),
                "asm_slop_end": _to_int(f[col["asm_slop_end"]]),
                "is_peak": _to_bool(f[col["is_peak"]]),
                "is_SD": _to_bool(f[col["is_SD"]]),
            }
    return idx


def peak_coverage_pct(
    align_start: int,
    align_end: int,
    asm_start: int,
    asm_end: int,
    slop_start: int,
) -> float:
    """
    Fraction of the inner peak covered by the alignment, as a percentage.

    All positions are in slop-local coords:
      inner_off = asm_start - slop_start
      inner_len = asm_end - asm_start
      align_overlap = max(0, min(align_end, inner_off+inner_len) - max(align_start, inner_off))
    """
    inner_off = asm_start - slop_start
    inner_len = asm_end - asm_start
    if inner_len <= 0:
        return 0.0
    lo = max(align_start, inner_off)
    hi = min(align_end, inner_off + inner_len)
    overlap = max(0, hi - lo)
    return (overlap / inner_len) * 100.0


def classify(q_meta: dict, t_meta: dict) -> str:
    """
    Classify a paralog edge between two RE_IDs.

    Within-haplotype-sample alignments (the primary case, produced by per-hap-sample
    minimap2) classify as 'paralog'. The others only occur if records are unified
    across haplotype samples at a later merge step.
    """
    # Trivial self-RE_ID (shouldn't happen with -X, but be defensive)
    if (
        q_meta["sample_id"] == t_meta["sample_id"]
        and q_meta["consensus_peak_id"] == t_meta["consensus_peak_id"]
    ):
        return "self"

    same_ind = q_meta["Individual_ID"] == t_meta["Individual_ID"]
    same_hap = q_meta["Haplotype"] == t_meta["Haplotype"]
    same_peak = q_meta["consensus_peak_id"] == t_meta["consensus_peak_id"]

    if same_ind and same_hap:
        # Different consensus peak in the same haplotype sample
        return "paralog" if not same_peak else "self"
    if same_ind and not same_hap:
        return "allelic-paralog" if not same_peak else "allelic"
    if not same_ind:
        return "ortholog-paralog" if not same_peak else "ortholog"
    return "unknown"


def main() -> None:
    sm = snakemake  # noqa: F821

    paf_in = sm.input.paf
    re_index_path = sm.input.re_index
    paf_out = sm.output.paf
    min_identity = float(sm.params.min_identity)
    min_coverage = float(sm.params.min_coverage)
    log_path = sm.log[0]

    idx = load_re_index(re_index_path)

    total = 0
    kept = 0
    unknown_re = 0
    classifications: dict[str, int] = {}

    reader = PAFReader(paf_in)
    out_records = []
    for rec in reader:
        total += 1
        q_meta = idx.get(rec.query_name)
        t_meta = idx.get(rec.target_name)
        if q_meta is None or t_meta is None:
            unknown_re += 1
            continue

        identity = rec.blast_identity()
        q_cov = peak_coverage_pct(
            rec.query_start,
            rec.query_end,
            q_meta["asm_start"],
            q_meta["asm_end"],
            q_meta["asm_slop_start"],
        )
        t_cov = peak_coverage_pct(
            rec.target_start,
            rec.target_end,
            t_meta["asm_start"],
            t_meta["asm_end"],
            t_meta["asm_slop_start"],
        )

        if identity < min_identity:
            continue
        if q_cov < min_coverage or t_cov < min_coverage:
            continue

        klass = classify(q_meta, t_meta)
        classifications[klass] = classifications.get(klass, 0) + 1

        # Tags (overwrite to a consistent schema)
        rec.tags["id"] = rec.query_name  # kept for downstream scripts that read id:Z
        rec.tags["cp"] = q_meta["consensus_peak_id"]
        rec.tags["tp"] = t_meta["consensus_peak_id"]
        rec.tags["qs"] = q_meta["sample_id"]
        rec.tags["ts"] = t_meta["sample_id"]
        rec.tags["qh"] = q_meta["haplotype_sample"]
        rec.tags["th"] = t_meta["haplotype_sample"]
        rec.tags["ct"] = klass
        rec.tags["bi"] = round(identity, 4)
        rec.tags["qc"] = round(q_cov, 2)
        rec.tags["tc"] = round(t_cov, 2)

        out_records.append(rec)
        kept += 1

    PAFReader.write_paf(out_records, paf_out)

    with open(log_path, "w") as lg:
        lg.write(f"input PAF: {paf_in}\n")
        lg.write(f"re_index:  {re_index_path}\n")
        lg.write(f"thresholds: identity>={min_identity}% both_peak_cov>={min_coverage}%\n\n")
        lg.write(f"total records:   {total}\n")
        lg.write(f"kept:            {kept}\n")
        lg.write(f"dropped (filters): {total - kept - unknown_re}\n")
        lg.write(f"dropped (unknown RE_ID): {unknown_re}\n")
        lg.write("\nclassification counts:\n")
        for k in sorted(classifications):
            lg.write(f"  {k}: {classifications[k]}\n")


if __name__ == "__main__":
    main()
