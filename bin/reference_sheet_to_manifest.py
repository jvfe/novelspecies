#!/usr/bin/env python3
"""Split a user-provided reference sheet into per-genus manifests."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import shutil

from path_utils import resolve_input_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-sheet", required=True, type=Path)
    parser.add_argument("--project-dir", type=Path, default=Path("."))
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    by_genus: dict[str, list[dict[str, str]]] = defaultdict(list)

    with args.reference_sheet.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"genus", "accession", "fasta"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError("reference_sheet must contain genus, accession and fasta columns")
        for row in reader:
            genus = row["genus"].strip()
            accession = row["accession"].strip()
            fasta = resolve_input_path(
                row["fasta"].strip(),
                args.reference_sheet.parent,
                args.project_dir,
            )
            if not genus or not accession:
                raise ValueError(f"Invalid reference sheet row: {row}")
            genus_dir = args.outdir / genus
            refs_dir = genus_dir / "staged_refs"
            refs_dir.mkdir(parents=True, exist_ok=True)
            dst = refs_dir / f"{accession}.fna"
            if fasta.resolve() != dst.resolve():
                shutil.copy2(fasta, dst)
            by_genus[genus].append(
                {
                    "genus": genus,
                    "accession": accession,
                    "organism_name": row.get("organism_name", "").strip(),
                    "strain": row.get("strain", "").strip(),
                    "assembly_level": row.get("assembly_level", "").strip(),
                    "genome_size": row.get("genome_size", "").strip(),
                    "is_type_strain": row.get("is_type_strain", "true").strip(),
                    "selection_mode": "reference_sheet",
                    "fasta": str(dst.resolve()),
                }
            )

    for genus, rows in by_genus.items():
        genus_dir = args.outdir / genus
        genus_dir.mkdir(parents=True, exist_ok=True)
        tsv_path = genus_dir / "references.tsv"
        list_path = genus_dir / "ref_list.txt"
        with tsv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with list_path.open("w") as handle:
            for row in rows:
                handle.write(f"{row['fasta']}\n")


if __name__ == "__main__":
    main()
