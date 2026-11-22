"""Validation helpers for deterministic outputs and retained references."""

from collections.abc import Mapping
import math


def assert_close_mapping(actual: Mapping[int, float], expected: Mapping[int, float], *,
                         atol: float = 1e-12, rtol: float = 1e-12) -> None:
    """Raise AssertionError if two numeric keyed mappings differ."""
    if actual.keys() != expected.keys():
        raise AssertionError(f"keys differ: {actual.keys()} != {expected.keys()}")
    for key in actual:
        if not math.isclose(actual[key], expected[key], abs_tol=atol, rel_tol=rtol):
            raise AssertionError(f"value differs at {key}: {actual[key]} != {expected[key]}")
