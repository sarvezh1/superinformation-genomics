# Superinformation genomics

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

## Scope and citation

The repository preserves the finalized analysis implementation and retained outputs; it is not a byte-for-byte reconstruction of the original raw-data extraction or Monte Carlo randomization workflow. See `docs/scientific_specification.md` and `docs/reproducibility_environments.md` for methods, provenance, and limitations.

If you use this software or its associated outputs, please cite the repository using the metadata in `CITATION.cff`.
