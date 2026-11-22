"""Canonical publication geometry from retained Superinformation values."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import importlib.metadata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

from .blocks import BLOCK_SIZES


@dataclass(frozen=True)
class GeometryResult:
    """Canonical feature-level matrix and derived deterministic analyses."""

    matrix: pd.DataFrame
    pca_coordinates: pd.DataFrame
    pca_variance: pd.DataFrame
    umap_coordinates: pd.DataFrame
    linkage_matrix: np.ndarray
    distance_matrix: pd.DataFrame
    cluster_assignments: pd.DataFrame
    cluster_composition: pd.DataFrame
    dendrogram_order: tuple[str, ...]


def feature_level_matrix(
    si: pd.DataFrame,
    *,
    species_column: str = "species",
    feature_column: str = "so_feature",
    block_column: str = "block_size",
    value_column: str = "superinformation",
    block_sizes: Sequence[int] = BLOCK_SIZES,
) -> pd.DataFrame:
    """Average SI across available species and retain complete feature rows."""
    required = {species_column, feature_column, block_column, value_column}
    missing = required.difference(si.columns)
    if missing:
        raise ValueError(f"SI table missing columns: {sorted(missing)}")
    sizes = tuple(block_sizes)
    if sizes != BLOCK_SIZES:
        raise ValueError(f"canonical geometry requires block sizes {BLOCK_SIZES}")
    grouped = (
        si.groupby([feature_column, block_column], sort=True)[value_column]
        .mean()
        .unstack()
        .reindex(columns=sizes)
    )
    return grouped.dropna(axis=0, how="any")


def _pca(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = PCA(n_components=2)
    coordinates = model.fit_transform(matrix.to_numpy(dtype=float))
    coords = pd.DataFrame(coordinates, index=matrix.index, columns=["PC1", "PC2"])
    coords.index.name = matrix.index.name
    variance = pd.DataFrame({
        "component": ["PC1", "PC2"],
        "explained_variance_ratio": model.explained_variance_ratio_,
        "explained_variance": model.explained_variance_,
    })
    return coords, variance


def _umap(matrix: pd.DataFrame) -> pd.DataFrame:
    import umap

    reducer = umap.UMAP(
        n_neighbors=10,
        min_dist=0.2,
        metric="euclidean",
        random_state=42,
        n_components=2,
    )
    coordinates = reducer.fit_transform(matrix.to_numpy(dtype=float))
    return pd.DataFrame(coordinates, index=matrix.index, columns=["UMAP1", "UMAP2"])


def analyze_geometry(si: pd.DataFrame) -> GeometryResult:
    """Run the approved 13-profile PCA, UMAP, and clustering analyses."""
    matrix = feature_level_matrix(si)
    pca_coordinates, pca_variance = _pca(matrix)
    umap_coordinates = _umap(matrix)
    values = matrix.to_numpy(dtype=float)
    distance = squareform(pdist(values, metric="euclidean"))
    distance_matrix = pd.DataFrame(distance, index=matrix.index, columns=matrix.index)
    distance_matrix.index.name = "so_feature"
    distance_matrix.columns.name = "so_feature"
    linkage_matrix = linkage(values, method="average", metric="euclidean")
    clusters = fcluster(linkage_matrix, t=4, criterion="maxclust")
    assignments = pd.DataFrame({"so_feature": matrix.index, "cluster": clusters})
    composition = (
        assignments.groupby("cluster", sort=True)
        .size()
        .reset_index(name="feature_count")
    )
    leaves = tuple(dendrogram(linkage_matrix, labels=matrix.index.tolist(), no_plot=True)["ivl"])
    return GeometryResult(
        matrix=matrix,
        pca_coordinates=pca_coordinates,
        pca_variance=pca_variance,
        umap_coordinates=umap_coordinates,
        linkage_matrix=linkage_matrix,
        distance_matrix=distance_matrix,
        cluster_assignments=assignments,
        cluster_composition=composition,
        dendrogram_order=leaves,
    )


def package_versions() -> dict[str, str]:
    """Return versions used for regenerated publication outputs."""
    names = {"numpy": "numpy", "pandas": "pandas", "scikit-learn": "scikit-learn",
             "scipy": "scipy", "umap-learn": "umap-learn"}
    return {key: importlib.metadata.version(dist) for key, dist in names.items()}
