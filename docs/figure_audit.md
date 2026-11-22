# Publication figure provenance

## Canonical geometry figures

The five standalone geometry figures are stored under
`results/geometry/figures/` and copied unchanged into `publication/overleaf/`.

| Figure asset | Numerical source | Status |
|---|---|---|
| `umap_feature_labelled.png/.svg` | `umap_coordinates.csv` | Regenerated canonical geometry rendering |
| `umap_broad_groups.png/.svg` | `umap_coordinates.csv`; visualization-group configuration | Regenerated canonical geometry rendering |
| `pca_feature_labelled.png/.svg` | `pca_coordinates.csv`; `pca_variance.csv` | Regenerated canonical geometry rendering |
| `hierarchical_dendrogram.png/.svg` | `clustering_linkage.csv`; `dendrogram_order.csv` | Regenerated canonical geometry rendering |
| `cluster_composition.png/.svg` | `cluster_assignments.csv`; `cluster_composition.csv`; dendrogram order | Regenerated canonical geometry rendering |

Rendering metadata verifies one-to-one feature labels, saved coordinate use,
source-data hashes, output hashes, and the absence of analysis fitting in the
rendering layer.

## Other publication figures

Figures 1-4 and Supplementary Figures S1-S7 are retained immutable publication
assets in `publication/overleaf/`.

- Figure 1 is an existing overview graphic. No complete public source/rendering
  script is available; it is intentionally retained unchanged and is not an
  exhaustive workflow diagram.
- Figures 2 and 3 are supported by the public retained SI table, but their exact
  rendering scripts are not part of the public repository.
- Figure 4 is supported by the retained SI and public null tables; its exact
  rendering script is not public.
- Supplementary Figure S1 is supported by the retained SI table.
- Supplementary Figure S2 is supported by the canonical 13 x 7 matrix and PCA
  outputs.
- Supplementary Figure S3 is supported by the fast and block-count-controlled
  null tables aggregated by block size.
- Supplementary Figures S4-S5 are supported by the retained population-variance
  tables.
- Supplementary Figure S6 is supported by the retained GC-controlled null.
- Supplementary Figure S7 is supported by the retained baseline tables.

These retained assets are publication evidence, not a claim that every raster
can be reproduced byte-for-byte from a public rendering script.
