#!/usr/bin/env python3
"""Render a clustered ANI heatmap for all genomes within a genus."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


# FastANI omits pairs below ~80% identity; treat those as below the species cutoff.
DEFAULT_MISSING_ANI = 70.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ani-tsv", required=True, type=Path)
    parser.add_argument("--genus", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cluster-threshold", type=float, default=95.0)
    parser.add_argument("--missing-ani", type=float, default=DEFAULT_MISSING_ANI)
    parser.add_argument("--versions", type=Path, default=None)
    return parser.parse_args()


def read_ani(path: Path, missing_ani: float) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    value_col = "ani_percent" if "ani_percent" in df.columns else "identity"
    if not {"query", "reference", value_col}.issubset(df.columns):
        raise ValueError("ANI table must contain query, reference and ani_percent/identity columns")

    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    df = df.dropna(subset=["query", "reference", value_col])
    if df.empty:
        raise ValueError("ANI table has no numeric pairwise values")

    pivot = df.pivot_table(index="query", columns="reference", values=value_col, aggfunc="max")
    labels = sorted(set(pivot.index.astype(str)).union(set(pivot.columns.astype(str))))
    matrix = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)
    for query in pivot.index:
        for ref in pivot.columns:
            value = pivot.loc[query, ref]
            if pd.isna(value):
                continue
            matrix.loc[str(query), str(ref)] = float(value)
            matrix.loc[str(ref), str(query)] = float(value)
    np.fill_diagonal(matrix.values, 100.0)
    matrix = matrix.fillna(missing_ani).clip(lower=0.0, upper=100.0)
    return matrix


def cluster_linkage(matrix: pd.DataFrame) -> np.ndarray:
    if len(matrix) < 2:
        return np.empty((0, 4))
    distance = (100.0 - matrix).clip(lower=0.0)
    np.fill_diagonal(distance.values, 0.0)
    condensed = squareform(distance.to_numpy(dtype=float), checks=False)
    condensed = np.nan_to_num(condensed, nan=30.0, posinf=30.0, neginf=0.0)
    condensed = np.clip(condensed, 0.0, None)
    if condensed.size == 0 or not np.isfinite(condensed).all():
        return np.empty((0, 4))
    return linkage(condensed, method="average")


def main() -> None:
    args = parse_args()
    matrix = read_ani(args.ani_tsv, args.missing_ani)
    linkage_matrix = cluster_linkage(matrix)

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
    fig.savefig(Path(f"{args.prefix}.ani_heatmap.png"), dpi=180)
    plt.close(fig)

    matrix.to_csv(f"{args.prefix}.ani_matrix.tsv", sep="\t")
    pd.DataFrame(linkage_matrix).to_csv(f"{args.prefix}.ani_linkage.tsv", sep="\t", index=False)

    clusters = (matrix >= args.cluster_threshold).astype(int)
    clusters.to_csv(f"{args.prefix}.ani_clusters_{args.cluster_threshold:.0f}.tsv", sep="\t")

    if args.versions:
        args.versions.write_text(
            f'"matplotlib": "{matplotlib.__version__}"\n'
            f'"seaborn": "{sns.__version__}"\n'
            f'"pandas": "{pd.__version__}"\n'
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
