"""Complete zero-offset nucleotide blocks used by the SI calculation."""

from collections.abc import Iterable, Iterator

BLOCK_SIZES: tuple[int, ...] = (20, 40, 60, 80, 100, 150, 200)
VALID_BASES = frozenset("ACGT")


def validate_block_size(block_size: int) -> int:
    """Validate a public study block size against the canonical set."""
    if not isinstance(block_size, int) or block_size not in BLOCK_SIZES:
        raise ValueError(f"block_size must be one of {BLOCK_SIZES}")
    return block_size


def validate_sequence(sequence: str) -> str:
    """Return an uppercase sequence or raise for empty/non-ACGT input."""
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    normalized = sequence.upper()
    if not normalized:
        raise ValueError("sequence must not be empty")
    if set(normalized) - VALID_BASES:
        raise ValueError("sequence contains non-ACGT symbols")
    return normalized


def _complete_blocks(sequence: str, block_size: int) -> Iterator[str]:
    """Internal generic complete-block helper; caller validates study sizes."""
    if not isinstance(block_size, int) or block_size <= 0:
        raise ValueError("block_size must be a positive integer")
    normalized = validate_sequence(sequence)
    for start in range(0, len(normalized) - block_size + 1, block_size):
        yield normalized[start:start + block_size]


def complete_blocks(sequence: str, block_size: int) -> Iterator[str]:
    """Yield canonical study blocks from offset zero; trailing sequence is ignored."""
    yield from _complete_blocks(sequence, validate_block_size(block_size))


def block_profile(sequence: str, block_sizes: Iterable[int] = BLOCK_SIZES) -> dict[int, tuple[str, ...]]:
    """Return complete blocks independently for each requested block size."""
    sizes = tuple(validate_block_size(size) for size in block_sizes)
    return {size: tuple(complete_blocks(sequence, size)) for size in sizes}
