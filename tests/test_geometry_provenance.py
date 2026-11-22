"""Checks for canonical geometry source-data provenance and figure inputs."""

from hashlib import sha256
from pathlib import Path
import tomllib
import ast
import importlib.util

import pytest

import pandas as pd

ROOT = Path(__file__).parents[1]
GEOMETRY = ROOT / "results" / "geometry"
FEATURES = [
    "CDS", "exon", "five_prime_UTR", "lnc_RNA", "mRNA", "miRNA",
    "ncRNA_gene", "pseudogene", "rRNA", "snRNA", "snoRNA",
    "three_prime_UTR", "transposable_element",
]
SOURCE_FILES = [
    "feature_level_si_matrix.csv", "pca_coordinates.csv", "pca_variance.csv",
    "umap_coordinates.csv", "umap_metadata.toml", "clustering_distance.csv",
    "clustering_linkage.csv", "cluster_assignments.csv", "cluster_composition.csv",
]
FIGURE_FILES = [
    "umap_feature_labelled.png", "umap_feature_labelled.svg",
    "umap_broad_groups.png", "umap_broad_groups.svg",
    "pca_feature_labelled.png", "pca_feature_labelled.svg",
    "hierarchical_dendrogram.png", "hierarchical_dendrogram.svg",
    "cluster_composition.png", "cluster_composition.svg",
]


def test_geometry_manifest_hashes_all_source_files():
    manifest = tomllib.loads((GEOMETRY / "GEOMETRY_MANIFEST.toml").read_text())
    for name in SOURCE_FILES:
        path = GEOMETRY / name
        assert path.is_file()
        assert manifest["files"][name] == sha256(path.read_bytes()).hexdigest()


def test_geometry_identifiers_and_cluster_composition():
    umap = pd.read_csv(GEOMETRY / "umap_coordinates.csv")
    pca = pd.read_csv(GEOMETRY / "pca_coordinates.csv")
    assignments = pd.read_csv(GEOMETRY / "cluster_assignments.csv")
    composition = pd.read_csv(GEOMETRY / "cluster_composition.csv")
    assert umap.so_feature.tolist() == FEATURES
    assert pca.so_feature.tolist() == FEATURES
    assert assignments.so_feature.tolist() == FEATURES
    assert assignments.so_feature.is_unique
    expected = assignments.groupby("cluster").size().rename("feature_count").reset_index()
    pd.testing.assert_frame_equal(composition, expected, check_dtype=False)


def test_figure_source_paths_are_valid():
    assert (ROOT / "scripts" / "generate_geometry_figures.py").is_file()
    for name in SOURCE_FILES:
        assert (GEOMETRY / name).is_file()


def test_display_mapping_is_one_to_one_and_complete():
    namespace = {"__file__": str(ROOT / "scripts" / "generate_geometry_figures.py")}
    tree = ast.parse((ROOT / "scripts" / "generate_geometry_figures.py").read_text())
    exec(compile(tree, "generate_geometry_figures.py", "exec"), namespace)
    mapping = namespace["DISPLAY_LABELS"]
    assert list(mapping) == FEATURES
    assert len(mapping) == len(set(mapping.values()))


def test_rendering_script_contains_no_scientific_fit_calls():
    source = (ROOT / "scripts" / "generate_geometry_figures.py").read_text()
    for forbidden in ("PCA(", "UMAP(", "linkage(", "pdist(", "fit_transform("):
        assert forbidden not in source
    assert "cluster_assignments.csv" in source


def test_authoritative_broad_group_mapping():
    config = tomllib.loads((ROOT / "config" / "visualization_groups.toml").read_text())
    assert config["feature_groups"] == {
        "CDS": "Coding", "exon": "Coding", "mRNA": "Coding",
        "five_prime_UTR": "UTR", "three_prime_UTR": "UTR",
        "lnc_RNA": "ncRNA", "miRNA": "ncRNA", "ncRNA_gene": "ncRNA",
        "rRNA": "ncRNA", "snRNA": "ncRNA", "snoRNA": "ncRNA",
        "pseudogene": "Pseudogene", "transposable_element": "Transposon",
    }


def test_standalone_figure_metadata_and_hashes():
    metadata = tomllib.loads((GEOMETRY / "figures" / "FIGURE_RENDERING_METADATA.toml").read_text())
    assert metadata["scientific_source_unchanged"] is True
    assert metadata["source_sha256_before"] == metadata["source_sha256_after"]
    for name, digest in metadata["source_sha256_after"].items():
        path = ROOT / name if name.startswith("config/") else GEOMETRY / name
        assert sha256(path.read_bytes()).hexdigest() == digest
    for name in FIGURE_FILES:
        path = GEOMETRY / "figures" / name
        assert path.is_file()
        assert metadata["output_sha256"][name] == sha256(path.read_bytes()).hexdigest()
    assert set(metadata["display_labels"]) == set(FEATURES)
    assert len(set(metadata["display_labels"].values())) == 13
    assert set(metadata["figures"]) == {
        "umap_feature_labelled", "umap_broad_groups", "pca_feature_labelled",
        "hierarchical_dendrogram", "cluster_membership",
    }


def test_cluster_membership_matrix_matches_saved_assignments():
    namespace = {"__file__": str(ROOT / "scripts" / "generate_geometry_figures.py")}
    tree = ast.parse((ROOT / "scripts" / "generate_geometry_figures.py").read_text())
    exec(compile(tree, "generate_geometry_figures.py", "exec"), namespace)
    order, clusters, matrix, counts = namespace["_cluster_assignment_matrix"]()
    saved_order = pd.read_csv(GEOMETRY / "dendrogram_order.csv").dendrogram_order.tolist()
    assignments = pd.read_csv(GEOMETRY / "cluster_assignments.csv").set_index("so_feature").cluster
    assert order == saved_order
    assert clusters == [1, 2, 3, 4]
    assert matrix.shape == (13, 4)
    assert (matrix.sum(axis=1) == 1).all()
    assert [clusters[column] for column in matrix.argmax(axis=1)] == [assignments[feature] for feature in order]
    assert counts == assignments.value_counts().sort_index().to_dict()


def test_geometry_check_mode_is_non_writing_and_detects_mismatch(tmp_path, monkeypatch):
    script = ROOT / "scripts" / "generate_geometry.py"
    spec = importlib.util.spec_from_file_location("generate_geometry", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    output_dir = tmp_path / "geometry"
    output_dir.mkdir()
    frozen = output_dir / "example.csv"
    frozen.write_bytes(b"value\n1\n")
    manifest = output_dir / "GEOMETRY_MANIFEST.toml"
    manifest.write_text("[files]\n", encoding="utf-8")
    reference = tmp_path / "reference.csv"
    reference.write_bytes(b"value\n1\n")

    monkeypatch.setattr(module, "_validate_environment", lambda: None)
    monkeypatch.setattr(module, "analyze_geometry", lambda frame: object())
    monkeypatch.setattr(
        module,
        "_regenerated_outputs",
        lambda result: {"example.csv": b"value\n1\n"},
    )
    monkeypatch.setattr(module, "_validate_manifest_hashes", lambda *args: None)
    monkeypatch.setattr(
        module.tomllib,
        "loads",
        lambda text: {"files": {}, "geometry": {"input_sha256": sha256(reference.read_bytes()).hexdigest()}},
    )

    before = {path.name: path.read_bytes() for path in output_dir.iterdir()}
    module.run_check(reference, output_dir, manifest)
    after = {path.name: path.read_bytes() for path in output_dir.iterdir()}
    assert after == before

    monkeypatch.setattr(
        module,
        "_regenerated_outputs",
        lambda result: {"example.csv": b"value\n2\n"},
    )
    with pytest.raises(module.GeometryCheckError, match="regenerated bytes differ"):
        module.run_check(reference, output_dir, manifest)
