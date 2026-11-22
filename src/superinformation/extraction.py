"""Schema-only extraction interface.

GFF3/FASTA extraction is intentionally not implemented: retained evidence does
not resolve all edge-case and provenance rules. This module only defines the
data contract that a future approved extractor must satisfy.
"""

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class FeatureInterval:
    """A parsed interval without prescribing extraction or aggregation policy."""

    species: str
    chrom: str
    start: int
    end: int
    strand: Literal["+", "-", ".", "?"]
    so_feature: str

    def __post_init__(self) -> None:
        if self.start < 1 or self.end < self.start:
            raise ValueError("interval must be one-based and start <= end")
        if not self.species or not self.chrom or not self.so_feature:
            raise ValueError("species, chrom, and so_feature are required")


class ExtractionNotImplementedError(NotImplementedError):
    """Raised until unresolved extraction provenance is closed."""
