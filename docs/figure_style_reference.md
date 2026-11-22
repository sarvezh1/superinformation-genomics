# Figure style reference

The standalone geometry figures use the visual conventions of the retained
publication assets: compact sans-serif typography, white backgrounds, light
axis treatment, restrained line weights, circular scatter markers, compact
legends, and minimal decoration. Standalone files omit panel letters and large
titles so that panels can be assembled in the publication source.

The broad visualization groups use a consistent five-color mapping recorded in
`config/visualization_groups.toml`. Color represents only the predefined
visualization group. Cluster membership uses separate cluster-specific markers
and does not reuse the broad-group encoding.

All geometry figures read saved numerical source data. Rendering does not fit
PCA or UMAP, calculate distances, recompute linkage, or derive new cluster
assignments. Source and output hashes are recorded in
`results/geometry/figures/FIGURE_RENDERING_METADATA.toml`.
