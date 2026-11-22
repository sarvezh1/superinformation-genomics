"""Render standalone publication figures from locked geometry source tables.

The renderer reads saved coordinates, linkage, order, and assignments. It does
not fit PCA/UMAP, calculate distances, or perform clustering.
"""

from __future__ import annotations

import csv
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import tomllib

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "results" / "geometry"
FIGURES = GEOMETRY / "figures"
GROUP_CONFIG = ROOT / "config" / "visualization_groups.toml"
RENDERER_VERSION = "3.0"

DISPLAY_LABELS = {
    "CDS": "CDS",
    "exon": "Exon",
    "five_prime_UTR": "5′ UTR",
    "lnc_RNA": "lncRNA",
    "mRNA": "mRNA",
    "miRNA": "miRNA",
    "ncRNA_gene": "ncRNA gene",
    "pseudogene": "Pseudogene",
    "rRNA": "rRNA",
    "snRNA": "snRNA",
    "snoRNA": "snoRNA",
    "three_prime_UTR": "3′ UTR",
    "transposable_element": "Transposable element",
}
FEATURES = list(DISPLAY_LABELS)
GROUP_ORDER = ["Coding", "UTR", "ncRNA", "Pseudogene", "Transposon"]
SOURCE_FILES = [
    "feature_level_si_matrix.csv",
    "pca_coordinates.csv",
    "pca_variance.csv",
    "umap_coordinates.csv",
    "clustering_linkage.csv",
    "dendrogram_order.csv",
    "cluster_assignments.csv",
    "cluster_composition.csv",
]
FIGURE_SOURCES = {
    "umap_feature_labelled": ["umap_coordinates.csv"],
    "umap_broad_groups": ["umap_coordinates.csv", "config/visualization_groups.toml"],
    "pca_feature_labelled": ["pca_coordinates.csv", "pca_variance.csv"],
    "hierarchical_dendrogram": ["clustering_linkage.csv", "dendrogram_order.csv"],
    "cluster_membership": ["cluster_assignments.csv", "cluster_composition.csv", "dendrogram_order.csv"],
}
FIGURE_OUTPUT_STEMS = {
    "umap_feature_labelled": "umap_feature_labelled",
    "umap_broad_groups": "umap_broad_groups",
    "pca_feature_labelled": "pca_feature_labelled",
    "hierarchical_dendrogram": "hierarchical_dendrogram",
    "cluster_membership": "cluster_composition",
}


def _read_csv(name: str) -> list[dict[str, str]]:
    with (GEOMETRY / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_coordinates(name: str, x: str, y: str) -> tuple[list[str], np.ndarray]:
    rows = _read_csv(name)
    labels = [row["so_feature"] for row in rows]
    values = np.asarray([[float(row[x]), float(row[y])] for row in rows])
    _validate_features(labels)
    return labels, values


def _validate_features(labels: list[str]) -> None:
    if labels != FEATURES or len(labels) != len(set(labels)):
        raise ValueError("feature identifiers do not match canonical order")


def _groups() -> tuple[dict[str, str], dict[str, str]]:
    config = tomllib.loads(GROUP_CONFIG.read_text(encoding="utf-8"))
    mapping = config["feature_groups"]
    colors = config["group_colors"]
    if set(mapping) != set(FEATURES) or set(mapping.values()) != set(GROUP_ORDER):
        raise ValueError("visualization-group configuration is incomplete")
    return mapping, colors


def _source_hashes() -> dict[str, str]:
    paths = {name: GEOMETRY / name for name in SOURCE_FILES}
    paths["config/visualization_groups.toml"] = GROUP_CONFIG
    return {name: sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}


def _style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 8,
        "axes.labelsize": 8.5,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def _clean_axes(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(width=0.8)


def _save(figure, stem: str) -> list[Path]:
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    outputs = [FIGURES / f"{stem}.png", FIGURES / f"{stem}.svg"]
    figure.savefig(outputs[0], dpi=300, bbox_inches="tight")
    figure.savefig(outputs[1], bbox_inches="tight")
    plt.close(figure)
    return outputs


def _annotate(axis, labels: list[str], values: np.ndarray, offsets: dict[str, tuple[int, int]]) -> None:
    for source, (x, y) in zip(labels, values):
        dx, dy = offsets.get(source, (5, 5))
        leader = abs(dx) > 10 or abs(dy) > 10
        axis.annotate(
            DISPLAY_LABELS[source],
            (x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.2,
            ha="left" if dx >= 0 else "right",
            va="bottom" if dy >= 0 else "top",
            arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.55} if leader else None,
        )


def _umap_limits(values: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return shared compact display limits for the two saved UMAP views."""
    x_min, y_min = values.min(axis=0)
    x_max, y_max = values.max(axis=0)
    x_pad = (x_max - x_min) * 0.10
    y_pad = (y_max - y_min) * 0.10
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)


def render_umap_feature_labelled() -> list[Path]:
    import matplotlib.pyplot as plt

    labels, values = _read_coordinates("umap_coordinates.csv", "UMAP1", "UMAP2")
    offsets = {
        "five_prime_UTR": (6, 6), "three_prime_UTR": (6, -8),
        "ncRNA_gene": (-10, 6), "transposable_element": (-10, 5),
        "lnc_RNA": (6, 3), "pseudogene": (6, 4),
    }
    figure, axis = plt.subplots(figsize=(4.8, 4.0), constrained_layout=True)
    axis.scatter(values[:, 0], values[:, 1], s=34, color="0.25", edgecolor="white", linewidth=0.5, zorder=3)
    _annotate(axis, labels, values, offsets)
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    x_limits, y_limits = _umap_limits(values)
    axis.set_xlim(*x_limits)
    axis.set_ylim(*y_limits)
    _clean_axes(axis)
    return _save(figure, "umap_feature_labelled")


def render_umap_broad_groups() -> list[Path]:
    import matplotlib.pyplot as plt

    labels, values = _read_coordinates("umap_coordinates.csv", "UMAP1", "UMAP2")
    mapping, colors = _groups()
    figure, axis = plt.subplots(figsize=(4.8, 4.0), constrained_layout=True)
    for group in GROUP_ORDER:
        indices = [i for i, label in enumerate(labels) if mapping[label] == group]
        axis.scatter(values[indices, 0], values[indices, 1], s=38, color=colors[group], label=group, edgecolor="white", linewidth=0.5)
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.legend(frameon=True, loc="upper right", borderpad=0.4, handletextpad=0.4, labelspacing=0.35)
    x_limits, y_limits = _umap_limits(values)
    axis.set_xlim(*x_limits)
    axis.set_ylim(*y_limits)
    _clean_axes(axis)
    return _save(figure, "umap_broad_groups")


def render_pca_feature_labelled() -> list[Path]:
    import matplotlib.pyplot as plt

    labels, values = _read_coordinates("pca_coordinates.csv", "PC1", "PC2")
    variance = {row["component"]: float(row["explained_variance_ratio"]) for row in _read_csv("pca_variance.csv")}
    offsets = {
        "five_prime_UTR": (12, 14), "three_prime_UTR": (-18, -20),
        "ncRNA_gene": (-6, -18), "transposable_element": (-8, 15),
        "lnc_RNA": (12, -12), "pseudogene": (-10, -16), "mRNA": (-10, 16),
        "exon": (-10, 10), "miRNA": (5, 5),
    }
    figure, axis = plt.subplots(figsize=(4.8, 4.0), constrained_layout=True)
    axis.scatter(values[:, 0], values[:, 1], s=34, color="0.25", edgecolor="white", linewidth=0.5, zorder=3)
    _annotate(axis, labels, values, offsets)
    axis.set_xlabel(f"PC1 ({variance['PC1'] * 100:.2f}%)")
    axis.set_ylabel(f"PC2 ({variance['PC2'] * 100:.2f}%)")
    axis.margins(x=0.12, y=0.12)
    _clean_axes(axis)
    return _save(figure, "pca_feature_labelled")


def render_hierarchical_dendrogram() -> list[Path]:
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram

    rows = _read_csv("clustering_linkage.csv")
    matrix = np.asarray([[float(row[key]) for key in ("left", "right", "distance", "cluster_size")] for row in rows])
    saved_order = [row["dendrogram_order"] for row in _read_csv("dendrogram_order.csv")]
    figure, axis = plt.subplots(figsize=(5.4, 4.6), constrained_layout=True)
    rendered = dendrogram(
        matrix,
        labels=[DISPLAY_LABELS[label] for label in FEATURES],
        orientation="left",
        link_color_func=lambda _: "#333333",
        above_threshold_color="#333333",
        leaf_font_size=7.5,
        ax=axis,
    )
    rendered_sources = [next(source for source, display in DISPLAY_LABELS.items() if display == label) for label in rendered["ivl"]]
    if rendered_sources != saved_order:
        raise ValueError("linkage leaf order differs from saved dendrogram order")
    axis.invert_yaxis()
    axis.set_xlabel("Euclidean distance")
    axis.set_ylabel("Feature")
    for collection in axis.collections:
        collection.set_linewidth(0.9)
    _clean_axes(axis)
    return _save(figure, "hierarchical_dendrogram")


def _cluster_assignment_matrix() -> tuple[list[str], list[int], np.ndarray, dict[int, int]]:
    """Arrange saved assignments in saved dendrogram order for rendering."""
    assignments = _read_csv("cluster_assignments.csv")
    labels = [row["so_feature"] for row in assignments]
    _validate_features(labels)
    expected_counts = {int(row["cluster"]): int(row["feature_count"]) for row in _read_csv("cluster_composition.csv")}
    assignment_by_feature = {row["so_feature"]: int(row["cluster"]) for row in assignments}
    order = [row["dendrogram_order"] for row in _read_csv("dendrogram_order.csv")]
    if set(order) != set(FEATURES) or len(order) != len(set(order)):
        raise ValueError("saved dendrogram order is not a permutation of canonical features")
    clusters = sorted(expected_counts)
    matrix = np.zeros((len(order), len(clusters)), dtype=int)
    cluster_to_column = {cluster: column for column, cluster in enumerate(clusters)}
    for row, feature in enumerate(order):
        matrix[row, cluster_to_column[assignment_by_feature[feature]]] = 1
    observed_counts = {cluster: int(matrix[:, column].sum()) for column, cluster in enumerate(clusters)}
    if observed_counts != expected_counts:
        raise ValueError("cluster assignments disagree with composition table")
    if not np.all(matrix.sum(axis=1) == 1):
        raise ValueError("each feature must have exactly one saved cluster assignment")
    return order, clusters, matrix, expected_counts


def render_cluster_membership() -> list[Path]:
    import matplotlib.pyplot as plt

    order, clusters, matrix, counts = _cluster_assignment_matrix()
    figure, axis = plt.subplots(figsize=(5.4, 4.6), constrained_layout=True)
    for row in range(len(order)):
        axis.hlines(row, -0.45, len(clusters) - 0.55, color="0.92", lw=0.45, zorder=0)
    assigned_rows, assigned_columns = np.where(matrix == 1)
    axis.scatter(
        assigned_columns,
        assigned_rows,
        s=42,
        marker="o",
        color="0.25",
        edgecolor="white",
        linewidth=0.55,
        zorder=3,
    )
    axis.set_xlim(-0.5, len(clusters) - 0.5)
    axis.set_ylim(len(order) - 0.5, -0.75)
    axis.set_xticks(
        range(len(clusters)),
        [f"Cluster {cluster}\n(n={counts[cluster]})" for cluster in clusters],
    )
    axis.xaxis.tick_top()
    axis.set_yticks(range(len(order)), [DISPLAY_LABELS[feature] for feature in order])
    axis.tick_params(axis="both", length=0, pad=5)
    for spine in axis.spines.values():
        spine.set_visible(False)
    return _save(figure, "cluster_composition")


def write_rendering_metadata(before: dict[str, str], after: dict[str, str], outputs: list[Path]) -> Path:
    if before != after:
        raise RuntimeError("scientific source data changed during rendering")
    mapping, colors = _groups()
    path = FIGURES / "FIGURE_RENDERING_METADATA.toml"
    lines = [
        f'renderer_version = "{RENDERER_VERSION}"',
        'figure_script = "scripts/generate_geometry_figures.py"',
        'style_evidence = "docs/figure_style_reference.md"',
        'scientific_source_unchanged = true',
        'recomputes_pca = false',
        'recomputes_umap = false',
        'recomputes_distances = false',
        'recomputes_clustering = false',
        f'generated_utc = {json.dumps(datetime.now(timezone.utc).isoformat())}',
        f'matplotlib = "{importlib.metadata.version("matplotlib")}"',
        f'scipy = "{importlib.metadata.version("scipy")}"',
        "", "[source_sha256_before]",
    ]
    lines.extend(f'{json.dumps(name)} = "{digest}"' for name, digest in before.items())
    lines += ["", "[source_sha256_after]"]
    lines.extend(f'{json.dumps(name)} = "{digest}"' for name, digest in after.items())
    lines += ["", "[output_sha256]"]
    lines.extend(f'{json.dumps(output.name)} = "{sha256(output.read_bytes()).hexdigest()}"' for output in outputs)
    lines += ["", "[display_labels]"]
    lines.extend(f'{json.dumps(source)} = {json.dumps(display)}' for source, display in DISPLAY_LABELS.items())
    lines += ["", "[broad_groups]"]
    lines.extend(f'{json.dumps(source)} = {json.dumps(group)}' for source, group in mapping.items())
    lines += ["", "[broad_group_colors]"]
    lines.extend(f'{json.dumps(group)} = {json.dumps(colors[group])}' for group in GROUP_ORDER)
    for figure_name, sources in FIGURE_SOURCES.items():
        stem = FIGURE_OUTPUT_STEMS[figure_name]
        lines += ["", f"[figures.{figure_name}]", f'conceptual_name = {json.dumps(figure_name.replace("_", " ").title())}', f'sources = {json.dumps(sources)}', f'outputs = ["{stem}.png", "{stem}.svg"]']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-membership-only", action="store_true")
    args = parser.parse_args()
    _style()
    before = _source_hashes()
    if args.cluster_membership_only:
        changed_outputs = render_cluster_membership()
    else:
        changed_outputs = []
        changed_outputs += render_umap_feature_labelled()
        changed_outputs += render_umap_broad_groups()
        changed_outputs += render_pca_feature_labelled()
        changed_outputs += render_hierarchical_dendrogram()
        changed_outputs += render_cluster_membership()
    outputs = [FIGURES / f"{stem}.{extension}" for stem in FIGURE_OUTPUT_STEMS.values() for extension in ("png", "svg")]
    after = _source_hashes()
    metadata = write_rendering_metadata(before, after, outputs)
    for output in [*changed_outputs, metadata]:
        print(output)


if __name__ == "__main__":
    main()
