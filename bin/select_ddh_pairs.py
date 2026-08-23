#!/usr/bin/env python3
"""Select top-N reference genomes per query based on FastANI."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ani-tsv", required=True, type=Path)
    parser.add_argument("--references-tsv", required=True, type=Path)
    parser.add_argument("--assignments-tsv", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    args = parse_args()
    ani_rows = read_table(args.ani_tsv)
    ref_rows = read_table(args.references_tsv)
    assignments = {row["sample"]: row for row in read_table(args.assignments_tsv)}

    ref_lookup = {row["accession"]: row for row in ref_rows}
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)

    for row in ani_rows:
        sample = row.get("query", row.get("sample", "")).strip()
        ref = row.get("reference", row.get("ref_accession", "")).strip()
        identity = float(row.get("ani_percent", row.get("identity", "0")))
        by_sample[sample].append({"sample": sample, "ref_accession": ref, "ani_percent": identity})

    output_rows: list[dict[str, str]] = []
    for sample, hits in by_sample.items():
        assignment = assignments.get(sample, {})
        query_fasta = assignment.get("fasta", "")
        genus = assignment.get("genus", "")
        ranked = sorted(hits, key=lambda x: x["ani_percent"], reverse=True)[: args.top_n]
        for rank, hit in enumerate(ranked, start=1):
            ref = ref_lookup.get(hit["ref_accession"], {})
            output_rows.append(
                {
                    "sample": sample,
                    "genus": genus,
                    "query_fasta": query_fasta,
                    "ref_accession": hit["ref_accession"],
                    "ref_organism": ref.get("organism_name", ""),
                    "reference_fasta": ref.get("fasta", ""),
                    "ani_percent": f"{hit['ani_percent']:.2f}",
                    "rank": str(rank),
                }
            )

    if not output_rows:
        raise ValueError("No ANI hits available to select dDDH pairs")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
