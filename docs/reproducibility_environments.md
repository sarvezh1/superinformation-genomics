# Reproducibility environments

## Validated UMAP regeneration environment

Canonical UMAP coordinates were regenerated and repeated identically in an
isolated environment using Python 3.11.9 and the exact pins in
`requirements-publication.txt`:

- NumPy 2.0.2
- SciPy 1.16.3
- scikit-learn 1.6.1
- Numba 0.60.0
- llvmlite 0.43.0
- PyNNDescent 0.5.13
- umap-learn 0.5.11
- joblib 1.4.2

The maximum and mean absolute coordinate differences between two repeated runs
were both zero. This is a validated regeneration environment, not a claim about
the original analysis environment.

## Publication figure rendering environment

The canonical standalone geometry figures were rendered from saved geometry
source data without recomputing PCA, UMAP, distances, linkage, or assignments.
The retained rendering metadata records Matplotlib 3.10.7 and SciPy 1.18.0.
Exact NumPy and pandas rendering versions were not recorded, so a complete
figure-rendering lock cannot be reconstructed without inference.

The PNG and SVG figures are therefore retained immutable publication assets
with hashes recorded in `FIGURE_RENDERING_METADATA.toml`. The UMAP regeneration
pins and figure-rendering metadata describe different environments and must not
be presented as identical.
