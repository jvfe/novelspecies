#!/usr/bin/env python3
"""Merge FastANI and dDDH outputs and apply species boundary rules."""

from __future__ import annotations

import argparse
import csv
import html
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments-tsv", required=True, type=Path)
    parser.add_argument("--query-ani-tsv", required=True, type=Path)
    parser.add_argument("--ddh-tsv", default=None, type=Path)
    parser.add_argument("--pairs-tsv", default=None, type=Path)
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument("--report-html", required=True, type=Path)
    parser.add_argument("--ani-threshold", type=float, default=95.0)
    parser.add_argument("--ddh-threshold", type=float, default=70.0)
    parser.add_argument("--borderline-low", type=float, default=94.0)
    parser.add_argument("--borderline-high", type=float, default=96.0)
    parser.add_argument("--ddh-method", default="local")
    return parser.parse_args()


def read_table(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def best_ani(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row.get("query", row.get("sample", "")).strip()
        identity = float(row.get("ani_percent", row.get("identity", "0")))
        ref = row.get("reference", row.get("ref_accession", "")).strip()
        organism = row.get("ref_organism", row.get("organism_name", "")).strip()
        current = best.get(sample)
        if current is None or identity > float(current["ani_percent"]):
            best[sample] = {
                "closest_type_strain": ref,
                "ref_organism": organism,
                "ani_percent": f"{identity:.2f}",
            }
    return best


def best_ddh(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row["sample"]
        ddh_raw = row.get("ddh_formula2", "NA")
        if ddh_raw in {"", "NA"}:
            continue
        ddh = float(ddh_raw)
        current = best.get(sample)
        if current is None or ddh > float(current["ddh_formula2"]):
            best[sample] = {
                "closest_type_strain": row.get("ref_accession", ""),
                "ref_organism": row.get("ref_organism", ""),
                "ddh_formula2": f"{ddh:.2f}",
                "ddh_method_used": row.get("ddh_method", "local_blast"),
            }
    return best


def classify(ani: float | None, ddh: float | None, args: argparse.Namespace) -> tuple[str, str]:
    if ani is None:
        return "insufficient_data", "false"

    borderline = args.borderline_low <= ani < args.borderline_high
    if ddh is None:
        if ani >= args.ani_threshold:
            verdict = "probable_same_species_ani_only"
        elif borderline:
            verdict = "borderline_ani_only"
        else:
            verdict = "probable_novel_species_ani_only"
        return verdict, "true" if ani < args.ani_threshold else "false"

    same_species = ani >= args.ani_threshold and ddh >= args.ddh_threshold
    novel = ani < args.ani_threshold and ddh < args.ddh_threshold

    if same_species:
        return "same_species", "false"
    if novel:
        return "novel_species_candidate", "true"
    if borderline:
        return "borderline_manual_review", "true"
    return "ambiguous_mixed_signals", "true"


def render_html(rows: list[dict[str, str]], output: Path) -> None:
    table_rows = []
    for row in rows:
        table_rows.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row[col]))}</td>" for col in rows[0].keys())
            + "</tr>"
        )
    headers = "".join(f"<th>{html.escape(col)}</th>" for col in rows[0].keys())
    content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Species delimitation summary</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ccc; padding: 0.4rem 0.6rem; text-align: left; }}
    th {{ background: #f5f5f5; }}
    .novel {{ background: #fff4e5; }}
    .same {{ background: #eef9ee; }}
  </style>
</head>
<body>
  <h1>Species delimitation summary</h1>
  <p>Species boundary rules: ANI &ge; 95% and Formula-2 dDDH &ge; 70% indicate the same species.</p>
  <table>
    <thead><tr>{headers}</tr></thead>
    <tbody>
      {''.join(table_rows)}
    </tbody>
  </table>
</body>
</html>
"""
    output.write_text(content)


def main() -> None:
    args = parse_args()
    assignments = read_table(args.assignments_tsv)
    ani_best = best_ani(read_table(args.query_ani_tsv))
    ddh_best = best_ddh(read_table(args.ddh_tsv))

    summary_rows: list[dict[str, str]] = []
    for assignment in assignments:
        sample = assignment["sample"]
        genus = assignment["genus"]
        ani = ani_best.get(sample, {})
        ddh = ddh_best.get(sample, {})

        ani_val = float(ani["ani_percent"]) if ani else None
        ddh_val = float(ddh["ddh_formula2"]) if ddh else None
        verdict, is_novel = classify(ani_val, ddh_val, args)

        closest = ani.get("closest_type_strain") or ddh.get("closest_type_strain", "NA")
        organism = ani.get("ref_organism") or ddh.get("ref_organism", "")

        summary_rows.append(
            {
                "sample": sample,
                "predicted_genus": genus,
                "closest_type_strain": closest,
                "closest_type_strain_name": organism,
                "top_fastani_identity": f"{ani_val:.2f}" if ani_val is not None else "NA",
                "top_ddh_formula2": f"{ddh_val:.2f}" if ddh_val is not None else "NA",
                "ddh_method_used": ddh.get("ddh_method_used", args.ddh_method),
                "taxonomic_verdict": verdict,
                "is_novel_candidate": is_novel,
            }
        )

    args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.summary_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    render_html(summary_rows, args.report_html)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
