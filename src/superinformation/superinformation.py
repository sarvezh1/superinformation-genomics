"""Histogram-based Superinformation calculation."""

import math
from collections.abc import Iterable

N_BINS = 50
ENTROPY_MIN = 0.0
ENTROPY_MAX = 2.0
BIN_WIDTH = (ENTROPY_MAX - ENTROPY_MIN) / N_BINS


def entropy_bin(entropy: float) -> int:
    """Map entropy to the specified 50-bin histogram with endpoint clamping."""
    index = int((float(entropy) - ENTROPY_MIN) / BIN_WIDTH)
    return min(max(index, 0), N_BINS - 1)


def entropy_histogram(entropies: Iterable[float]) -> tuple[int, ...]:
    """Return the integer histogram for an entropy iterable."""
    bins = [0] * N_BINS
    for value in entropies:
        bins[entropy_bin(value)] += 1
    return tuple(bins)


def superinformation(entropies: Iterable[float]) -> float:
    """Compute entropy of the entropy histogram; empty input returns zero."""
    bins = entropy_histogram(entropies)
    total = sum(bins)
    if total == 0:
        return 0.0
    return float(-sum((count / total) * math.log2(count / total)
                      for count in bins if count))
