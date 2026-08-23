#!/usr/bin/env python3
"""Normalize FastANI output to a standard TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--comparison-type", default="query_vs_reference")
    return parser.parse_args()


def genome_label(path_str: str) -> str:
    path = Path(path_str)
    stem = path.stem
    if stem.endswith("_genomic"):
        stem = stem[: -len("_genomic")]
    return stem


def main() -> None:
    args = parse_args()
    rows = []
    with args.input.open() as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            query = genome_label(parts[0])
            reference = genome_label(parts[1])
            ani = float(parts[2]) * 100.0
            rows.append(
                {
                    "query": query,
                    "reference": reference,
                    "ani_percent": f"{ani:.2f}",
                    "comparison_type": args.comparison_type,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["query", "reference", "ani_percent", "comparison_type"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
