# Scientific specification

## Scope and authority

This document records the accepted scientific definition implemented in the
publication repository for *Multiscale superinformation reveals structured
organization of genomic sequence features*. The finalized manuscript and
Supplementary Information define the study; retained reference outputs provide
the numerical regression targets.

Superinformation was introduced by Bose and Chouhan (2011). The present study
extends its application across genomic feature classes, species, and sequence
resolutions.

Provenance terms used here are:

- **FROZEN**: established by the finalized publication specification.
- **RECOVERED**: supported by authenticated original source or direct numerical
  reproduction.
- **REFERENCE-OUTPUT ONLY**: the retained numerical output is authoritative,
  but complete original generation state is unavailable.
- **UNRESOLVED**: a detail cannot be recovered and is not inferred.

This repository is a reproducible implementation of the finalized analysis
specification using retained reference outputs. It is not a claim of
byte-identical recovery of the original software or runtime state.

## Inputs and feature set

The analysis used Ensembl release 111 assemblies and annotations for human
(GRCh38), mouse (GRCm39), rat (mRatBN7.2), zebrafish (GRCz11), fly (BDGP6.46),
worm (WBcel235), and yeast (R64-1-1).

Fourteen feature types occur in the retained SI data:

`CDS`, `exon`, `five_prime_UTR`, `lnc_RNA`, `mRNA`, `miRNA`, `ncRNA_gene`,
`pseudogene`, `rRNA`, `snRNA`, `snoRNA`, `tRNA`, `three_prime_UTR`, and
`transposable_element`.

The final SI table contains 578 species-feature-block-size observations across
seven block sizes: 20, 40, 60, 80, 100, 150, and 200 bp.

### Input and extraction limitations

The intended assemblies, release, and source URLs were recovered, but the exact
downloaded FASTA/GFF3 bytes, retrieval timestamps, and input checksums are not
available. The recovered extractor used one-based inclusive GFF3 coordinates
through a `start - 1:end` sequence slice, reverse-complemented minus-strand
records, and treated qualifying GFF3 rows as independent records.

Sequence records were uppercased and excluded in full if they contained any
symbol outside A, C, G, and T. This eligibility check occurred before block
construction. Exact behavior for every malformed or out-of-range annotation
edge case, and a byte-level manifest of the final extracted records, remain
unresolved. The public `extraction.py` therefore defines an interface/schema
only and does not claim exact raw-data reconstruction.

## Blocks, entropy, and Superinformation

For each retained record and each block size independently, blocks start at
offset zero, do not overlap, do not cross record boundaries, and have exactly
the requested length. A terminal remainder shorter than the block size is
discarded; no padding or partial block is used.

For a complete A/C/G/T block, base-2 Shannon entropy is

\[
H=-\sum_{x\in\{A,C,G,T\}} p_x\log_2 p_x,
\]

with zero-probability terms contributing zero. Entropy ranges from 0 to 2 bits.

Superinformation is the Shannon entropy of the empirical distribution of block
entropy values. The distribution uses 50 equal-width bins spanning 0 to 2 bits.
The implementation maps a value using a 0.04-bit bin width and clamps the upper
endpoint, including exactly 2 bits, into the final bin. Zero-count bins do not
contribute to the histogram entropy.

SI is calculated independently at every evaluated block size. If a record has
no complete block at a requested size, the study-level observation is missing;
it is not replaced by numerical zero or interpolated.

## Profiles and aggregation

The retained SI table contains 88 represented species-feature trajectories;
86 contain at least two evaluated resolutions. Seventy-two of those 86 are
non-increasing, or 83.72%. The decline is therefore broad but not universal.

Feature-level geometry profiles are arithmetic means across available species
at each feature and block size. Only features complete at all seven sizes are
retained for geometry. The incomplete tRNA profile is excluded, producing an
unscaled 13 x 7 matrix. Missing values are not imputed or interpolated.

## Canonical geometry

PCA uses two components on the unscaled 13 x 7 matrix. Explained-variance
ratios are:

- PC1: `0.9334130852249626` (93.341%)
- PC2: `0.04517565313436558` (4.518%)
- combined: `0.9785887383593281` (97.859%)

UMAP uses the same matrix with two components, Euclidean distance,
`n_neighbors=10`, `min_dist=0.2`, and `random_state=42`. The saved coordinates
are the canonical regenerated publication coordinates. The validated
regeneration environment is recorded in `requirements-publication.txt` and
the geometry metadata.

Hierarchical clustering uses Euclidean distance and average linkage on the
same unscaled matrix. A four-cluster `maxclust` cut gives:

- Cluster 1: CDS, rRNA
- Cluster 2: snRNA, snoRNA
- Cluster 3: exon, five_prime_UTR, lnc_RNA, mRNA, ncRNA_gene, pseudogene,
  three_prime_UTR, transposable_element
- Cluster 4: miRNA

Cluster identifiers summarize quantitative SI-profile similarity and are not
biological classifications.

## Null models

The block-count-controlled random null generates equal-probability A/C/G/T
blocks for each species-feature-block-size case, matching the observed number
of contributing blocks. Cases with fewer than 50 observed blocks are excluded.

The species-level GC-controlled null generates 50,000 blocks per species and
block size using the retained species GC fraction, with equal A/T probabilities
and equal C/G probabilities. The fast null generates 50,000 uniform A/C/G/T
blocks per block size, independently of species and feature.

No explicit seed or runtime random-number-generator state was retained for the
original null runs. The public null tables are therefore immutable retained
reference outputs and are not regenerated by the repository. They support the
reported summaries: block-count-null observed-minus-null mean about 1.38 bits,
median about 1.46 bits, and maximum about 2.68 bits; GC-null mean about 0.43
bits and median about 0.54 bits; fast-null correlation about 0.999995 with the
block-size-aggregated block-count null and mean absolute difference about
0.005 bits.

## Baseline descriptors

The baseline analysis compares feature-level mean SI with mononucleotide
Shannon entropy, overlapping dinucleotide entropy, and GC content. The retained
baseline input contains 89 species-feature rows and uses up to 1,000 sequences
per row. Baseline rows are distinct from the 578 SI observations.

Recovered correlations are mononucleotide entropy `r=-0.3224`, `R^2=0.1039`;
dinucleotide entropy `r=-0.2574`, `R^2=0.0662`; and GC content `r=-0.3194`,
`R^2=0.1020`. The original operation selecting or truncating rows at the
1,000-sequence cap and its record manifest are unavailable. These tables are
retained reference outputs.

## Entropy variance

The final supplementary variance products use population variance across
contributing block-entropy values:

\[
\operatorname{Var}(H)=\frac{1}{n}\sum_i H_i^2-\left(\frac{1}{n}\sum_i H_i\right)^2.
\]

Both 578-row variance tables are retained publication outputs. The repository
does not recompute the underlying block-entropy dataset.

## CIS and IPS

The Characteristic Information Scale is the first B50 threshold satisfying

\[
SI(B_{50})=0.5\,SI(20).
\]

Crossings between adjacent available resolutions use linear interpolation. If
the threshold is not reached within 20 to 200 bp, CIS is not reached. The
retained table contains 22 finite CIS values. Behavior for an unobserved
`SI(20)=0` profile remains unresolved and is never selected implicitly.

Raw Information Persistence Score is trapezoidal integration of SI over the
available evaluated block sizes within 20 to 200 bp. Complete trajectories span
the full window; incomplete trajectories use only their available points and
are not extended or interpolated. Normalized IPS applies the same integration
after division by `SI(20)`. Raw IPS has units bit-bp and normalized IPS has
units bp. Behavior for normalized IPS when `SI(20)=0` remains unresolved.

## Publication assets and reproducibility

SI, CIS, IPS, null, baseline, and entropy-variance CSVs are retained reference
outputs with hashes recorded in public manifests. The geometry matrix, PCA,
UMAP, clustering, and standalone geometry figures are regenerated publication
outputs with their inputs, parameters, environments, and hashes recorded under
`results/geometry/`.

Figures 1-4 and the supplementary figures are retained immutable publication
assets. Complete public rendering pipelines are not available for all of them;
their status is documented in `docs/figure_audit.md`. Figure 1 is retained as
an intentionally unchanged overview asset and is not represented as an
exhaustive workflow diagram.

The regression suite validates deterministic core behavior, all retained SI,
CIS, and IPS rows, canonical PCA targets, geometry provenance, cluster coverage,
and figure hashes. The validated release check reports 33 passing tests with no
failures or skips.
