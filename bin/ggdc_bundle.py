#!/usr/bin/env python3
"""Bundle query/reference FASTA pairs for optional GGDC web submission."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tarfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-tsv", required=True, type=Path)
    parser.add_argument("--email", default="")
    parser.add_argument("--bundle-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--instructions", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.bundle_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    with args.pairs_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)

    if not rows:
        raise ValueError("No GGDC comparison pairs provided")

    manifest_rows = []
    for idx, row in enumerate(rows, start=1):
        job_id = f"job_{idx:03d}_{row['sample']}"
        job_dir = args.bundle_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        query_src = Path(row["query_fasta"])
        ref_src = Path(row["reference_fasta"])
        query_dst = job_dir / f"{row['sample']}.fna"
        ref_dst = job_dir / f"{row['ref_accession']}.fna"
        shutil.copy2(query_src, query_dst)
        shutil.copy2(ref_src, ref_dst)
        manifest_rows.append(
            {
                "job_id": job_id,
                "sample": row["sample"],
                "genus": row.get("genus", ""),
                "ref_accession": row["ref_accession"],
                "ref_organism": row.get("ref_organism", ""),
                "query_fasta": str(query_dst.name),
                "reference_fasta": str(ref_dst.name),
                "ggdc_formula": "2",
                "ggdc_url": "https://ggdc.dsmz.de/ggdc.php",
            }
        )

    with args.manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    with args.archive.open("wb") as archive_handle:
        with tarfile.open(fileobj=archive_handle, mode="w:gz") as tar:
            tar.add(args.bundle_dir, arcname=args.bundle_dir.name)

    email_line = args.email if args.email else "<your-email>"
    args.instructions.write_text(
        "\n".join(
            [
                "GGDC manual submission bundle",
                "=============================",
                "",
                "GGDC does not expose a stable public REST API. This bundle prepares",
                "one query/reference FASTA pair per job for upload at:",
                "  https://ggdc.dsmz.de/ggdc.php",
                "",
                "For each job directory:",
                "  1. Upload the query FASTA as 'Query genome (FASTA)'",
                "  2. Upload the reference FASTA as 'Reference genomes (FASTA)'",
                "  3. Select Formula 2 (recommended for draft genomes)",
                "  4. Provide contact email: {}".format(email_line),
                "  5. Download the emailed results and ingest with:",
                "       --ggdc_results /path/to/ggdc_results/",
                "",
                "Manifest: {}".format(args.manifest.name),
                "Archive:  {}".format(args.archive.name),
                "",
            ]
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
