import math

import pytest

from superinformation.blocks import BLOCK_SIZES, complete_blocks, validate_sequence
from superinformation.entropy import shannon_entropy
from superinformation.persistence import compute_cis, compute_ips
from superinformation.profiles import sequence_profile
from superinformation.superinformation import entropy_bin, superinformation


def test_blocks_are_zero_offset_and_drop_trailing_sequence():
    assert list(complete_blocks("ACGT" * 25 + "AC", 20)) == ["ACGT" * 5] * 5


def test_public_blocks_enforce_canonical_sizes():
    assert BLOCK_SIZES == (20, 40, 60, 80, 100, 150, 200)
    with pytest.raises(ValueError):
        list(complete_blocks("A" * 40, 4))


def test_invalid_sequence_is_not_silently_accepted():
    with pytest.raises(ValueError):
        validate_sequence("ACN")


def test_entropy_known_values():
    assert shannon_entropy("A" * 20) == 0.0
    assert shannon_entropy("ACGT" * 5) == 2.0


def test_histogram_endpoint_clamping_and_empty_si():
    assert entropy_bin(0.0) == 0
    assert entropy_bin(2.0) == 49
    assert entropy_bin(-1.0) == 0
    assert entropy_bin(3.0) == 49
    assert superinformation([]) == 0.0


def test_profile_has_requested_resolutions():
    profile = sequence_profile("ACGT" * 100, (20, 40, 200))
    assert set(profile) == {20, 40, 200}
    assert all(math.isfinite(value) for value in profile.values())


def test_profile_preserves_missing_observation_as_nan():
    profile = sequence_profile("ACGT" * 5, (20, 40))
    assert math.isfinite(profile[20])
    assert math.isnan(profile[40])


def test_cis_interpolates_first_crossing():
    assert compute_cis([(20, 4), (40, 3), (60, 1)]) == pytest.approx(50.0)


def test_cis_requires_explicit_block_size_20():
    assert math.isnan(compute_cis([(40, 4), (60, 1)]))


def test_cis_non_monotonic_and_no_crossing():
    assert compute_cis([(20, 4), (40, 3), (60, 1), (80, 3)]) == pytest.approx(50.0)
    assert math.isnan(compute_cis([(20, 4), (40, 3), (60, 2.1)]))


def test_cis_requires_explicit_zero_initial_policy():
    with pytest.raises(ValueError):
        compute_cis([(20, 0), (40, 1)])
    assert math.isnan(compute_cis([(20, 0), (40, 1)], zero_initial="nan"))
    assert compute_cis([(20, 0), (40, 1)], zero_initial="first") == 20


def test_ips_available_span_and_one_point():
    assert compute_ips([(20, 2), (100, 1), (200, 0)]) == pytest.approx(170.0)
    assert compute_ips([(100, 3)]) == 0.0
    assert compute_ips([(20, 2), (200, 1)], normalized=True) == pytest.approx(135.0)


def test_normalized_ips_requires_explicit_zero_policy():
    with pytest.raises(ValueError):
        compute_ips([(20, 0), (40, 1)], normalized=True)


def test_normalized_ips_requires_block_size_20():
    assert math.isnan(compute_ips([(40, 2), (100, 1)], normalized=True))


def test_ips_filters_outside_window_without_fabrication():
    assert compute_ips([(10, 100), (20, 2), (100, 1), (200, 0), (220, 100)]) == pytest.approx(170.0)
