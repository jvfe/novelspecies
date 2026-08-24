#!/usr/bin/env python3
"""Filter NCBI datasets output to a genus-specific reference manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
from pathlib import Path


TYPE_MATERIAL_LABELS = {
    "TYPE_MATERIAL",
    "type material",
    "assembly from type material",
}

SKIP_FASTA_SUBSTRINGS = (
    "cds_from_genomic",
    "rna_from_genomic",
    "translated_cds",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-tsv", required=True, type=Path)
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--genus", required=True)
    parser.add_argument("--type-material-only", action="store_true")
    parser.add_argument("--allow-non-type-fallback", action="store_true")
    parser.add_argument("--max-references", type=int, default=50)
    parser.add_argument("--accessions-out", type=Path, default=None)
    parser.add_argument("--downloaded-accessions", type=Path, default=None)
    parser.add_argument("--staged-dir", type=Path, default=None)
    parser.add_argument("--references-tsv", type=Path, default=None)
    parser.add_argument("--ref-list", type=Path, default=None)
    return parser.parse_args()


def discover_data_dir(start: Path) -> Path:
    if start.is_dir():
        if any(start.rglob("*_genomic.fna*")) or any(start.rglob("*.fna*")):
            return start
        nested = start / "ncbi_dataset" / "data"
        if nested.is_dir():
            return nested
        data = start / "data"
        if data.is_dir():
            return data
    for candidate in Path(".").rglob("assembly_data_report.jsonl"):
        return candidate.parent
    for candidate in Path(".").rglob("*_genomic.fna"):
        return candidate.parent.parent if candidate.parent.name.startswith(("GCF_", "GCA_")) else candidate.parent
    return start


def is_skipped_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(token in name for token in SKIP_FASTA_SUBSTRINGS)


def find_fasta(data_dir: Path, accession: str) -> Path | None:
    accession = accession.strip()
    patterns = [
        f"{accession}*_genomic.fna",
        f"{accession}*_genomic.fna.gz",
        f"{accession}*.fna",
        f"{accession}*.fna.gz",
        f"{accession}*.fasta",
        f"{accession}*.fasta.gz",
        f"{accession}*.fa",
        f"{accession}*.fa.gz",
        "genomic.fna",
        "genomic.fna.gz",
    ]
    matches: list[Path] = []
    accession_dir = data_dir / accession
    search_roots = [accession_dir, data_dir] if accession_dir.is_dir() else [data_dir]
    for root in search_roots:
        for pattern in patterns:
            for hit in sorted(root.rglob(pattern)):
                if hit.is_file() and not is_skipped_fasta(hit):
                    matches.append(hit)
        if matches:
            break
    if not matches:
        return None
    uncompressed = [path for path in matches if not path.name.endswith(".gz")]
    return (uncompressed or matches)[0]


def stage_fasta(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.name.endswith(".gz"):
        with gzip.open(src, "rb") as incoming, dst.open("wb") as outgoing:
            shutil.copyfileobj(incoming, outgoing)
    else:
        shutil.copy2(src, dst)


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


def select_rows(rows: list[dict[str, str]], args: argparse.Namespace) -> tuple[list[dict[str, str]], str]:
    selected = rows
    selection_mode = "all_assemblies"
    if args.downloaded_accessions and args.downloaded_accessions.exists():
        wanted = {
            line.strip()
            for line in args.downloaded_accessions.read_text().splitlines()
            if line.strip()
        }
        selected = [row for row in rows if row["accession"] in wanted]
        selection_mode = "downloaded_accessions"
        if selected:
            return rank_rows(selected)[: args.max_references], selection_mode

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
    return rank_rows(selected)[: args.max_references], selection_mode


def describe_data_dir(data_dir: Path) -> str:
    files = [str(path) for path in sorted(data_dir.rglob("*")) if path.is_file()]
    preview = "\n".join(files[:40])
    extra = f"\n... ({len(files) - 40} more)" if len(files) > 40 else ""
    return f"{data_dir} ({len(files)} files)\n{preview}{extra}"


def main() -> None:
    args = parse_args()
    rows = read_report(args.report_tsv)
    if not rows:
        raise ValueError(f"No assemblies found in report for genus '{args.genus}'")

    selected, selection_mode = select_rows(rows, args)

    if args.accessions_out:
        args.accessions_out.parent.mkdir(parents=True, exist_ok=True)
        args.accessions_out.write_text("\n".join(row["accession"] for row in selected) + "\n")
        if args.data_dir is None:
            return

    if args.data_dir is None or args.references_tsv is None or args.ref_list is None:
        raise ValueError("Provide --data-dir, --references-tsv and --ref-list after genomes are downloaded")

    data_dir = discover_data_dir(args.data_dir)
    references: list[dict[str, str]] = []
    staged_dir = args.staged_dir
    if staged_dir:
        staged_dir.mkdir(parents=True, exist_ok=True)

    for row in selected:
        fasta = find_fasta(data_dir, row["accession"])
        if fasta is None:
            continue
        fasta_out = fasta
        if staged_dir:
            fasta_out = staged_dir / f"{row['accession']}.fna"
            stage_fasta(fasta, fasta_out)
            fasta_out = Path("staged_refs") / fasta_out.name
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
                "fasta": fasta_out.as_posix(),
            }
        )

    if not references:
        available = []
        for fasta in data_dir.rglob("*"):
            if fasta.is_file() and fasta.suffix in {".fna", ".fa", ".fasta", ".gz"}:
                available.append(fasta.name)
        extra = ""
        if available:
            extra = f" FASTAs present: {', '.join(available[:20])}."
        raise ValueError(
            f"No downloadable FASTA files found for genus '{args.genus}' "
            f"in {describe_data_dir(data_dir)}.{extra} "
            f"Looked for {len(selected)} selected accessions "
            f"({', '.join(row['accession'] for row in selected[:10])})."
        )

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
