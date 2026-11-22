# Retained reference fixtures

This directory contains three authenticated retained outputs used by the
regression suite:

- `superinformation_values.csv`: 578 final species-feature-block-size SI rows;
- `CIS_values.csv`: 88 species-feature CIS/B50 rows;
- `IPS_values.csv`: 88 species-feature raw and normalized IPS rows.

Tests match rows by stable identifiers and compare numerical fields with an
absolute tolerance of `5e-7`, matching the six-decimal precision of the
retained CSV values. The canonical PCA regression uses an absolute tolerance
of `1e-12`. File hashes, schemas, row counts, and scientific roles are recorded
in `REFERENCE_MANIFEST.toml`.

These files are retained reference outputs. The repository validates them and
uses the SI table as the input to regenerated canonical geometry; it does not
claim to reproduce the unavailable original raw-data extraction run.
