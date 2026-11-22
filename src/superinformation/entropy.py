"""Shannon entropy for validated nucleotide blocks."""

import math
from collections import Counter


def shannon_entropy(sequence: str) -> float:
    """Compute base-2 Shannon entropy in bits for an A/C/G/T sequence."""
    if not sequence:
        raise ValueError("entropy is undefined for an empty sequence")
    counts = Counter(sequence)
    if set(counts) - set("ACGT"):
        raise ValueError("entropy input contains non-ACGT symbols")
    total = len(sequence)
    return float(-sum((n / total) * math.log2(n / total) for n in counts.values() if n))
