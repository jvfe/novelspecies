#!/usr/bin/env python3
"""Create sample-named symlinks and a FastANI query list."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--genus", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    queries_dir = args.output.parent / "queries"
    queries_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with args.assignments.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("genus") == args.genus:
                rows.append(row)

    if not rows:
        raise ValueError(f"No assignments found for genus '{args.genus}'")

    with args.output.open("w") as out_handle:
        for row in rows:
            sample = row["sample"]
            fasta = Path(row["fasta"])
            link = queries_dir / f"{sample}.fna"
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to(fasta.resolve())
            out_handle.write(f"{link.absolute()}\n")


if __name__ == "__main__":
    main()
