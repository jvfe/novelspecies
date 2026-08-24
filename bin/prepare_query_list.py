#!/usr/bin/env python3
"""Stage sample-named query FASTAs and write a relative FastANI query list."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True, type=Path, help="TSV with sample and fasta_name columns")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def stage_fasta(src: Path, dst: Path) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()
    queries_dir = args.output.parent / "queries"
    queries_dir.mkdir(parents=True, exist_ok=True)

    with args.mapping.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise ValueError("No query genomes provided for this genus")

    with args.output.open("w") as out_handle:
        for row in rows:
            sample = row["sample"].strip()
            src = Path(row["fasta_name"].strip())
            if not src.is_file():
                raise FileNotFoundError(f"Staged query FASTA not found: {src}")
            dst = queries_dir / f"{sample}.fna"
            stage_fasta(src.resolve(), dst)
            out_handle.write(f"{dst.as_posix()}\n")


if __name__ == "__main__":
    main()
