#!/usr/bin/env python3
"""Assign a target bacterial genus to each query genome."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from path_utils import resolve_input_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--gtdbtk-summary", type=Path, default=None)
    parser.add_argument("--genus-map", type=Path, default=None)
    parser.add_argument("--project-dir", type=Path, default=Path("."))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--versions", type=Path, default=None)
    return parser.parse_args()


def read_samplesheet(path: Path, project_dir: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "sample" not in reader.fieldnames or "fasta" not in reader.fieldnames:
            raise ValueError("Samplesheet must contain 'sample' and 'fasta' columns")
        rows = []
        for row in reader:
            sample = row["sample"].strip()
            fasta_raw = row["fasta"].strip()
            if not sample or not fasta_raw:
                raise ValueError("Empty sample or fasta entry in samplesheet")
            fasta = resolve_input_path(fasta_raw, path.parent, project_dir)
            genus = row.get("genus", "").strip()
            rows.append({"sample": sample, "fasta": str(fasta), "genus": genus})
        return rows


def read_genus_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "sample" not in reader.fieldnames or "genus" not in reader.fieldnames:
            raise ValueError("Genus map must contain 'sample' and 'genus' columns")
        for row in reader:
            sample = row["sample"].strip()
            genus = row["genus"].strip()
            if sample and genus:
                mapping[sample] = genus
    return mapping


def parse_gtdb_genus(classification: str) -> str | None:
    match = re.search(r"g__([^;]+)", classification)
    if not match:
        return None
    genus = match.group(1).strip()
    if genus in {"", "Unknown", "unknown"}:
        return None
    return genus


def read_gtdbtk_summary(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("GTDB-Tk summary is empty")
        genome_col = None
        for candidate in ("user_genome", "genome", "genome_id"):
            if candidate in reader.fieldnames:
                genome_col = candidate
                break
        if genome_col is None:
            raise ValueError("GTDB-Tk summary must contain 'user_genome' or 'genome'")
        if "classification" not in reader.fieldnames:
            raise ValueError("GTDB-Tk summary must contain 'classification'")
        for row in reader:
            genome_id = row[genome_col].strip()
            genus = parse_gtdb_genus(row["classification"])
            if genus:
                mapping[genome_id] = genus
                base = genome_id.split(".")[0]
                mapping.setdefault(base, genus)
    return mapping


def normalize_sample_id(sample: str) -> str:
    return sample.split(".")[0]


def main() -> None:
    args = parse_args()
    rows = read_samplesheet(args.samplesheet, args.project_dir)
    genus_map = read_genus_map(args.genus_map) if args.genus_map else {}
    gtdb_map = read_gtdbtk_summary(args.gtdbtk_summary) if args.gtdbtk_summary else {}

    assignments: list[dict[str, str]] = []
    for row in rows:
        sample = row["sample"]
        genus = row["genus"]
        source = "samplesheet" if genus else ""

        if not genus and sample in genus_map:
            genus = genus_map[sample]
            source = "genus_map"
        if not genus and sample in gtdb_map:
            genus = gtdb_map[sample]
            source = "gtdbtk"
        if not genus:
            base = normalize_sample_id(sample)
            if base in gtdb_map:
                genus = gtdb_map[base]
                source = "gtdbtk"

        if not genus:
            raise ValueError(
                f"Could not assign genus for sample '{sample}'. "
                "Provide genus in the samplesheet, --genus_map, or --gtdbtk_summary."
            )

        assignments.append(
            {
                "sample": sample,
                "fasta": row["fasta"],
                "genus": genus,
                "genus_source": source,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "fasta", "genus", "genus_source"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(assignments)

    if args.versions:
        args.versions.write_text(f'"{Path(__file__).stem}":\n    assign_genus: 1.0.0\n')


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
