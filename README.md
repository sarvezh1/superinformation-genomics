# Superinformation genomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24437905.svg)](https://doi.org/10.5281/zenodo.24437905)
[![GitHub Release](https://img.shields.io/github/v/release/sarvezh1/superinformation-genomics)](https://github.com/sarvezh1/superinformation-genomics/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Code and retained analysis outputs accompanying *Multiscale superinformation reveals structured organization of genomic sequence features*.

This project applies the Superinformation measure introduced by Bose and Chouhan (2011) to genomic sequence features across species and block sizes. It includes the deterministic analysis implementation, retained reference outputs, and the canonical feature-profile geometry results.

An initial version of this work was presented as a poster by **Sarvesh Galgale** at the **Genomics India Conference 2025**, held at the **Indian Institute of Science (IISc), Bengaluru**.

## Repository contents

- `src/` — Superinformation analysis implementation
- `results/` — retained reference outputs and geometry results
- `scripts/` — analysis and figure-generation scripts
- `tests/` — regression and provenance tests
- `docs/` — scientific specification and reproducibility notes
- `publication/` — publication figure assets

## Install and test

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

Some geometry regeneration steps require the pinned environment documented in `docs/reproducibility_environments.md`. The retained outputs and their manifests provide the reference values used by the regression tests.

## Scope

The repository preserves the finalized analysis implementation and retained outputs; it is not a byte-for-byte reconstruction of the original raw-data extraction or Monte Carlo randomization workflow. See `docs/scientific_specification.md` and `docs/reproducibility_environments.md` for methods, provenance, and limitations.

## Citation

If you use this repository or its associated data and analyses, please cite the archived release:

**Multiscale superinformation reveals structured organization of genomic sequence features.** (2026). *Superinformation Genomics* (Version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.24437905

The citation metadata is also available in [`CITATION.cff`](CITATION.cff).