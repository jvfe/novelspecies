#!/usr/bin/env python3
"""Estimate GGDC formula-2 (d4) digital DDH from BLASTn HSPs."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Hsp:
    qstart: int
    qend: int
    sstart: int
    send: int
    length: int
    nident: int
    bitscore: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genus", required=True)
    parser.add_argument("--ref-accession", required=True)
    parser.add_argument("--ref-organism", default="")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--blast-args", default="-task blastn -word_size 38 -dust no -evalue 1e-2 -max_target_seqs 1000000")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def run_blast(blastn: str, blast_args: str, query: Path, reference: Path, out_path: Path) -> None:
    cmd = [
        blastn,
        *blast_args.split(),
        "-query",
        str(query),
        "-subject",
        str(reference),
        "-outfmt",
        "6 qstart qend sstart send length nident bitscore",
        "-out",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def parse_blast(path: Path) -> list[Hsp]:
    hsps: list[Hsp] = []
    if not path.exists() or path.stat().st_size == 0:
        return hsps
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            hsps.append(
                Hsp(
                    qstart=int(parts[0]),
                    qend=int(parts[1]),
                    sstart=int(parts[2]),
                    send=int(parts[3]),
                    length=int(parts[4]),
                    nident=int(parts[5]),
                    bitscore=float(parts[6]),
                )
            )
    return hsps


def merge_intervals(hsps: list[Hsp], axis: str) -> list[Hsp]:
    if not hsps:
        return []
    key_start = (lambda h: h.qstart) if axis == "query" else (lambda h: h.sstart)
    key_end = (lambda h: h.qend) if axis == "query" else (lambda h: h.send)
    ordered = sorted(hsps, key=lambda h: (key_start(h), -h.bitscore))
    merged: list[Hsp] = []
    cur_start = key_start(ordered[0])
    cur_end = key_end(ordered[0])
    cur_ident = ordered[0].nident
    cur_len = ordered[0].length
    cur_score = ordered[0].bitscore

    for hsp in ordered[1:]:
        start = key_start(hsp)
        end = key_end(hsp)
        if start <= cur_end:
            if hsp.bitscore > cur_score:
                cur_start, cur_end = start, end
                cur_ident, cur_len, cur_score = hsp.nident, hsp.length, hsp.bitscore
            else:
                cur_end = max(cur_end, end)
        else:
            merged.append(
                Hsp(cur_start, cur_end, 0, 0, cur_len, cur_ident, cur_score)
            )
            cur_start, cur_end = start, end
            cur_ident, cur_len, cur_score = hsp.nident, hsp.length, hsp.bitscore
    merged.append(Hsp(cur_start, cur_end, 0, 0, cur_len, cur_ident, cur_score))
    return merged


def formula_d4(hsps: list[Hsp]) -> tuple[float, float, int]:
    if not hsps:
        return math.nan, math.nan, 0
    total_ident = sum(h.nident for h in hsps)
    total_len = sum(h.length for h in hsps)
    if total_len == 0:
        return math.nan, math.nan, 0
    similarity = total_ident / total_len
    distance = 1.0 - similarity
    return similarity, distance, len(hsps)


def distance_to_ddh(distance: float) -> float:
    """Approximate GGDC formula-2 conversion (log-GLM; Meier-Kolthoff et al. 2013)."""
    if math.isnan(distance) or distance <= 0:
        return math.nan
    log_distance = math.log(distance)
    ddh = 100.0 / (1.0 + math.exp(-(3.916 * log_distance + 3.1)))
    return max(0.0, min(100.0, ddh))


def main() -> None:
    args = parse_args()
    prefix = f"{args.sample}_vs_{args.ref_accession}"
    q_vs_r = Path(f"{prefix}_q_vs_r.blast")
    r_vs_q = Path(f"{prefix}_r_vs_q.blast")

    run_blast(args.blastn, args.blast_args, args.query, args.reference, q_vs_r)
    run_blast(args.blastn, args.blast_args, args.reference, args.query, r_vs_q)

    hsps = parse_blast(q_vs_r) + parse_blast(r_vs_q)
    hsps = merge_intervals(hsps, axis="query")
    similarity, distance, hsp_count = formula_d4(hsps)

    if math.isnan(similarity):
        ddh = math.nan
        coverage = 0.0
    else:
        ddh = distance_to_ddh(distance)
        coverage = similarity * 100.0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "genus",
                "ref_accession",
                "ref_organism",
                "formula",
                "hsp_count",
                "hsp_identity_fraction",
                "gbdp_distance_d4",
                "ddh_formula2",
                "hsp_coverage_percent",
                "ddh_method",
            ]
        )
        writer.writerow(
            [
                args.sample,
                args.genus,
                args.ref_accession,
                args.ref_organism,
                "GGDC_formula_2_d4",
                hsp_count,
                f"{similarity:.6f}" if not math.isnan(similarity) else "NA",
                f"{distance:.6f}" if not math.isnan(distance) else "NA",
                f"{ddh:.2f}" if not math.isnan(ddh) else "NA",
                f"{coverage:.2f}" if not math.isnan(similarity) else "NA",
                "local_blast",
            ]
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
