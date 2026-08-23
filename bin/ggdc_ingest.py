#!/usr/bin/env python3
"""Parse GGDC result files (email attachments or downloaded reports)."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--pairs-tsv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_ggdc_text(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    patterns = {
        "ddh_formula2": r"Formula\s*2[^\n]*?([0-9]+\.[0-9]+)\s*%",
        "distance_d4": r"distance[^\n]*?([0-9]+\.[0-9]+)",
        "probability_same_species": r"Probability that DDH > 70%[^\n]*?([0-9]+\.[0-9]+)\s*%",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            values[key] = match.group(1)

    # Fallback: look for explicit "DDH estimate" lines near formula 2 blocks
    if "ddh_formula2" not in values:
        match = re.search(r"DDH estimate \(GLM-based\):\s*([0-9]+\.[0-9]+)%", text, flags=re.IGNORECASE)
        if match:
            values["ddh_formula2"] = match.group(1)
    return values


def load_pairs(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_result_file(results_dir: Path, sample: str, accession: str) -> Path | None:
    candidates = sorted(results_dir.rglob("*"))
    sample_lower = sample.lower()
    accession_lower = accession.lower()
    for candidate in candidates:
        if not candidate.is_file():
            continue
        name = candidate.name.lower()
        if accession_lower in name or sample_lower in name:
            if candidate.suffix.lower() in {".txt", ".eml", ".html", ".tsv", ".csv"}:
                return candidate
    return None


def main() -> None:
    args = parse_args()
    pairs = load_pairs(args.pairs_tsv)
    rows: list[dict[str, str]] = []

    for pair in pairs:
        result_file = find_result_file(args.results_dir, pair["sample"], pair["ref_accession"])
        if result_file is None:
            rows.append(
                {
                    "sample": pair["sample"],
                    "genus": pair.get("genus", ""),
                    "ref_accession": pair["ref_accession"],
                    "ref_organism": pair.get("ref_organism", ""),
                    "ddh_formula2": "NA",
                    "distance_d4": "NA",
                    "probability_same_species": "NA",
                    "ddh_method": "ggdc_web",
                    "result_file": "NA",
                }
            )
            continue

        parsed = parse_ggdc_text(result_file.read_text(errors="ignore"))
        rows.append(
            {
                "sample": pair["sample"],
                "genus": pair.get("genus", ""),
                "ref_accession": pair["ref_accession"],
                "ref_organism": pair.get("ref_organism", ""),
                "ddh_formula2": parsed.get("ddh_formula2", "NA"),
                "distance_d4": parsed.get("distance_d4", "NA"),
                "probability_same_species": parsed.get("probability_same_species", "NA"),
                "ddh_method": "ggdc_web",
                "result_file": str(result_file),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
