# Geometry publication notes

The canonical geometry uses 13 complete feature-level profiles from the
14-feature study. Each value is the arithmetic mean across available species at
one block size. The seven columns are 20, 40, 60, 80, 100, 150, and 200 bp. The
matrix contains no missing values and is not scaled.

The retained geometry features, in matrix order, are `CDS`, `exon`,
`five_prime_UTR`, `lnc_RNA`, `mRNA`, `miRNA`, `ncRNA_gene`, `pseudogene`,
`rRNA`, `snRNA`, `snoRNA`, `three_prime_UTR`, and `transposable_element`.
The incomplete tRNA profile is excluded.

PCA explained-variance ratios are `0.9334130852249626` for PC1 and
`0.04517565313436558` for PC2, totaling `0.9785887383593281`.

UMAP uses two components, Euclidean distance, `n_neighbors=10`,
`min_dist=0.2`, and `random_state=42`. The coordinates are canonical regenerated
publication results from the frozen matrix.

Clustering uses Euclidean distance, average linkage, and a four-cluster
`maxclust` cut. Full-precision inputs, coordinates, distances, linkage,
assignments, and manifests are under `results/geometry/`.

The standalone publication candidates are under
`results/geometry/figures/`. Copies of their PNG files under
`publication/overleaf/` are used by the publication source.
