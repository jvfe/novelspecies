#!/usr/bin/env python3
"""Render a clustered ANI heatmap for all genomes within a genus."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ani-tsv", required=True, type=Path)
    parser.add_argument("--genus", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cluster-threshold", type=float, default=95.0)
    parser.add_argument("--versions", type=Path, default=None)
    return parser.parse_args()


def read_ani(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if {"query", "reference", "ani_percent"}.issubset(df.columns):
        pivot = df.pivot(index="query", columns="reference", values="ani_percent")
    elif {"query", "reference", "identity"}.issubset(df.columns):
        pivot = df.pivot(index="query", columns="reference", values="identity")
    else:
        raise ValueError("ANI table must contain query, reference and ani_percent/identity columns")
    labels = sorted(set(pivot.index).union(set(pivot.columns)))
    matrix = pd.DataFrame(np.nan, index=labels, columns=labels)
    for query in pivot.index:
        for ref in pivot.columns:
            value = pivot.loc[query, ref]
            matrix.loc[query, ref] = value
            matrix.loc[ref, query] = value
    np.fill_diagonal(matrix.values, 100.0)
    return matrix.astype(float)


def write_newick(labels: list[str], linkage_matrix: np.ndarray, out_path: Path) -> None:
    # Minimal placeholder: scipy linkage to newick requires external helper;
    # store linkage matrix as TSV for downstream tree building.
    out_path.write_text("label\n" + "\n".join(labels) + "\n")


def main() -> None:
    args = parse_args()
    matrix = read_ani(args.ani_tsv)
    distance = 100.0 - matrix
    np.fill_diagonal(distance.values, 0.0)
    condensed = squareform(distance.values, checks=False)
    linkage_matrix = linkage(condensed, method="average")

    cluster_labels = []
    for label in matrix.index:
        cluster_labels.append(label)
    order = linkage_matrix[:, :2].astype(int)

    fig, ax = plt.subplots(figsize=(max(6, len(matrix) * 0.45), max(5, len(matrix) * 0.45)))
    sns.heatmap(
        matrix,
        cmap="viridis",
        vmin=70,
        vmax=100,
        square=True,
        linewidths=0.2,
        cbar_kws={"label": "ANI (%)"},
        ax=ax,
    )
    ax.set_title(f"{args.genus} ANI heatmap")
    ax.set_xlabel("Genome")
    ax.set_ylabel("Genome")
    fig.tight_layout()
    png_path = Path(f"{args.prefix}.ani_heatmap.png")
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    matrix.to_csv(f"{args.prefix}.ani_matrix.tsv", sep="\t")
    pd.DataFrame(linkage_matrix).to_csv(f"{args.prefix}.ani_linkage.tsv", sep="\t", index=False)

    clusters = (matrix >= args.cluster_threshold).astype(int)
    clusters.to_csv(f"{args.prefix}.ani_clusters_{args.cluster_threshold:.0f}.tsv", sep="\t")

    if args.versions:
        args.versions.write_text(
            '"matplotlib": "' + str(getattr(matplotlib, "__version__", "unknown")) + '"\n'
            '"seaborn": "' + str(getattr(sns, "__version__", "unknown")) + '"\n'
            '"pandas": "' + str(pd.__version__) + '"\n'
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
