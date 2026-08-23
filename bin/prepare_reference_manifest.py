#!/usr/bin/env python3
"""Filter NCBI datasets output to a genus-specific reference manifest."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


TYPE_MATERIAL_LABELS = {
    "TYPE_MATERIAL",
    "type material",
    "assembly from type material",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-tsv", required=True, type=Path)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--genus", required=True)
    parser.add_argument("--type-material-only", action="store_true")
    parser.add_argument("--allow-non-type-fallback", action="store_true")
    parser.add_argument("--max-references", type=int, default=50)
    parser.add_argument("--references-tsv", required=True, type=Path)
    parser.add_argument("--ref-list", required=True, type=Path)
    return parser.parse_args()


def find_fasta(data_dir: Path, accession: str) -> Path | None:
    patterns = [
        f"{accession}*_genomic.fna",
        f"{accession}*.fna",
        f"{accession}*.fasta",
        f"{accession}*.fa",
    ]
    for pattern in patterns:
        matches = sorted(data_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def is_type_material(row: dict[str, str]) -> bool:
    label = row.get("type_material_label", "").strip()
    display = row.get("type_material_display", "").strip()
    if label.upper() == "TYPE_MATERIAL":
        return True
    if display.lower() in {x.lower() for x in TYPE_MATERIAL_LABELS}:
        return True
    return "type material" in display.lower()


def read_report(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for raw in reader:
            rows.append(
                {
                    "accession": raw.get("Assembly Accession", raw.get("accession", "")).strip(),
                    "organism_name": raw.get("Organism Name", raw.get("organism-name", "")).strip(),
                    "strain": raw.get("Organism Infraspecific Names Strain", raw.get("organism-infraspecific-strain", "")).strip(),
                    "assembly_level": raw.get("Assembly Level", raw.get("assminfo-level", "")).strip(),
                    "genome_size": raw.get("Assembly Stats Total Sequence Length", raw.get("assmstats-total-sequence-len", "")).strip(),
                    "source_database": raw.get("Source Database", raw.get("source_database", "")).strip(),
                    "type_material_label": raw.get("Type Material Label", raw.get("type_material-label", "")).strip(),
                    "type_material_display": raw.get("Type Material Display Text", raw.get("type_material-display_text", "")).strip(),
                }
            )
        return [row for row in rows if row["accession"]]


def rank_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    level_rank = {
        "complete genome": 0,
        "chromosome": 1,
        "scaffold": 2,
        "contig": 3,
    }

    def sort_key(row: dict[str, str]) -> tuple[int, int, str]:
        level = row.get("assembly_level", "").lower()
        size = int(row.get("genome_size") or 0)
        return (level_rank.get(level, 99), -size, row["accession"])

    return sorted(rows, key=sort_key)


def main() -> None:
    args = parse_args()
    rows = read_report(args.report_tsv)
    if not rows:
        raise ValueError(f"No assemblies found in report for genus '{args.genus}'")

    selected = rows
    selection_mode = "all_assemblies"
    if args.type_material_only:
        typed = [row for row in rows if is_type_material(row)]
        if typed:
            selected = typed
            selection_mode = "type_material"
        elif args.allow_non_type_fallback:
            selected = rows
            selection_mode = "fallback_all_assemblies"
        else:
            raise ValueError(
                f"No type-material assemblies found for genus '{args.genus}'. "
                "Use --allow-non-type-fallback or disable --type_material_only."
            )

    selected = rank_rows(selected)[: args.max_references]
    references: list[dict[str, str]] = []
    for row in selected:
        fasta = find_fasta(args.data_dir, row["accession"])
        if fasta is None:
            continue
        references.append(
            {
                "genus": args.genus,
                "accession": row["accession"],
                "organism_name": row["organism_name"],
                "strain": row["strain"],
                "assembly_level": row["assembly_level"],
                "genome_size": row["genome_size"],
                "is_type_strain": str(is_type_material(row)).lower(),
                "selection_mode": selection_mode,
                "fasta": str(fasta.resolve()),
            }
        )

    if not references:
        raise ValueError(f"No downloadable FASTA files found for genus '{args.genus}'")

    args.references_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.references_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(references[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(references)

    with args.ref_list.open("w") as handle:
        for ref in references:
            handle.write(f"{ref['fasta']}\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
