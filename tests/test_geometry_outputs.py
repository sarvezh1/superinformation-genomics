"""Regression checks for regenerated canonical geometry source data."""

from pathlib import Path
import csv
import tomllib

import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]
GEOMETRY = ROOT / "results" / "geometry"
FEATURES = [
    "CDS", "exon", "five_prime_UTR", "lnc_RNA", "mRNA", "miRNA",
    "ncRNA_gene", "pseudogene", "rRNA", "snRNA", "snoRNA",
    "three_prime_UTR", "transposable_element",
]
BLOCKS = [20, 40, 60, 80, 100, 150, 200]


def test_canonical_matrix_shape_features_and_completeness():
    matrix = pd.read_csv(GEOMETRY / "feature_level_si_matrix.csv")
    assert list(matrix.columns) == ["so_feature", *[str(x) for x in BLOCKS]]
    assert matrix.shape == (13, 8)
    assert matrix.so_feature.tolist() == FEATURES
    assert matrix.iloc[:, 1:].notna().all().all()


def test_pca_regression_and_feature_identifiers():
    variance = pd.read_csv(GEOMETRY / "pca_variance.csv")
    assert variance.component.tolist() == ["PC1", "PC2"]
    np.testing.assert_allclose(
        variance.explained_variance_ratio.to_numpy(),
        [0.9334130852249626, 0.04517565313436555],
        rtol=0,
        atol=1e-12,
    )
    assert pd.read_csv(GEOMETRY / "pca_coordinates.csv").so_feature.tolist() == FEATURES


def test_umap_shape_finiteness_and_metadata():
    coords = pd.read_csv(GEOMETRY / "umap_coordinates.csv")
    assert coords.shape == (13, 3)
    assert coords.so_feature.tolist() == FEATURES
    assert np.isfinite(coords[["UMAP1", "UMAP2"]].to_numpy()).all()
    metadata = tomllib.loads((GEOMETRY / "umap_metadata.toml").read_text())
    assert {key: metadata["umap"][key] for key in
            ("n_neighbors", "min_dist", "metric", "random_state", "n_components")} == {
        "n_neighbors": 10,
        "min_dist": 0.2,
        "metric": "euclidean",
        "random_state": 42,
        "n_components": 2,
    }
    for package in ("numpy", "scipy", "scikit_learn", "numba", "llvmlite", "pynndescent", "umap_learn", "joblib"):
        assert package in metadata["umap"]


def test_umap_repeatability_recorded():
    metadata = tomllib.loads((GEOMETRY / "umap_metadata.toml").read_text())
    assert metadata["umap_repeat_max_abs_difference"] == 0.0
    assert metadata["umap_repeat_mean_abs_difference"] == 0.0


def test_clustering_parameter_metadata_and_feature_coverage():
    manifest = tomllib.loads((GEOMETRY / "GEOMETRY_MANIFEST.toml").read_text())
    assert manifest["clustering"] == {
        "metric": "euclidean",
        "linkage": "average",
        "criterion": "maxclust",
        "requested_clusters": 4,
    }
    assignments = pd.read_csv(GEOMETRY / "cluster_assignments.csv")
    assert assignments.so_feature.tolist() == FEATURES
    assert assignments.so_feature.is_unique
    assert len(assignments) == 13
