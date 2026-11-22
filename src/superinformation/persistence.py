"""CIS/B50 and raw/normalized IPS over evaluated block sizes."""

from collections.abc import Iterable, Sequence
import math


def _ordered(points: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    values = sorted((float(x), float(y)) for x, y in points)
    if not values:
        raise ValueError("at least one evaluated point is required")
    return values


def compute_cis(points: Iterable[tuple[float, float]], *, zero_initial: str = "raise") -> float:
    """Compute the first B50 crossing with explicit zero-initial policy.

    ``zero_initial='raise'`` is the default because the recovered sources
    disagree for this unobserved edge case. ``'nan'`` and ``'first'`` are
    explicit caller choices and are never selected implicitly.
    """
    values = _ordered(points)
    si20 = [y for x, y in values if x == 20.0]
    if not si20:
        return math.nan
    x0, y0 = 20.0, si20[0]
    if not math.isfinite(y0) or y0 <= 0:
        if zero_initial == "nan":
            return math.nan
        if zero_initial == "first":
            return x0
        raise ValueError("SI(20)<=0 policy is unresolved; choose zero_initial explicitly")
    target = 0.5 * y0
    crossing_values = [(x, y) for x, y in values if x >= 20.0]
    for (x1, y1), (x2, y2) in zip(crossing_values, crossing_values[1:]):
        if not (math.isfinite(y1) and math.isfinite(y2)):
            continue
        if y2 <= target:
            if y2 == target or y2 == y1:
                return x2
            return x1 + (target - y1) * (x2 - x1) / (y2 - y1)
    return math.nan


def _trapezoid(points: Sequence[tuple[float, float]]) -> float:
    return float(sum((y1 + y2) * (x2 - x1) / 2
                     for (x1, y1), (x2, y2) in zip(points, points[1:])))


def compute_ips(points: Iterable[tuple[float, float]], *, normalized: bool = False,
                lower: float = 20.0, upper: float = 200.0) -> float:
    """Integrate available evaluated points within the approved study window."""
    values = [(x, y) for x, y in _ordered(points) if lower <= x <= upper]
    if not values:
        raise ValueError("no points within the integration window")
    if normalized:
        si20 = [y for x, y in values if x == 20.0]
        if not si20:
            return math.nan
        first = si20[0]
        if first == 0 or not math.isfinite(first):
            raise ValueError("SI(20)=0 policy is unresolved for normalized IPS")
        values = [(x, y / first) for x, y in values]
    return _trapezoid(values)
