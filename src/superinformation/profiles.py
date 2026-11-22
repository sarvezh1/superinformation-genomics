"""Deterministic SI profiles across the approved block sizes."""

from collections.abc import Iterable
import math

from .blocks import BLOCK_SIZES, complete_blocks, validate_sequence
from .entropy import shannon_entropy
from .superinformation import superinformation


def sequence_profile(sequence: str, block_sizes: Iterable[int] = BLOCK_SIZES) -> dict[int, float]:
    """Compute one SI value per block size for a validated sequence."""
    normalized = validate_sequence(sequence)
    result: dict[int, float] = {}
    for size in block_sizes:
        blocks = tuple(complete_blocks(normalized, size))
        if not blocks:
            result[size] = math.nan
            continue
        result[size] = superinformation(shannon_entropy(block) for block in blocks)
    return result
