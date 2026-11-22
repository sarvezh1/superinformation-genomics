"""Regression checks for reference files added under tests/fixtures/reference/."""

from pathlib import Path
import hashlib
import zipfile

import pandas as pd
import pytest

from superinformation.persistence import compute_cis, compute_ips


REFERENCE = Path(__file__).parent / "fixtures" / "reference"


def _zip_csv(archive: Path, suffix: str) -> pd.DataFrame:
    with zipfile.ZipFile(archive) as bundle:
        name = next(name for name in bundle.namelist() if name.endswith(suffix))
        return pd.read_csv(bundle.open(name))


def _require_reference(name: str) -> Path:
    path = REFERENCE / name
    if not path.exists():
        pytest.skip(f"reference fixture not installed: {path}")
    return path


def _assert_rows(table: pd.DataFrame, expected: list[dict[str, object]], value_columns: list[str], atol: float = 5e-7) -> None:
    keys = [key for key in ("species", "so_feature", "block_size") if key in table.columns]
    for row in expected:
        mask = pd.Series(True, index=table.index)
        for key in keys:
            if key in row:
                mask &= table[key].eq(row[key])
        matches = table.loc[mask]
        assert len(matches) == 1, f"expected exactly one row for {row}"
        actual = matches.iloc[0]
        for column in value_columns:
            if pd.isna(row[column]):
                assert pd.isna(actual[column])
            else:
                assert actual[column] == pytest.approx(row[column], abs=atol)


def test_retained_si_schema_and_all_numerical_rows():
    table = pd.read_csv(_require_reference("superinformation_values.csv"))
    assert list(table.columns) == ["species", "so_feature", "block_size", "superinformation", "n_blocks"]
    assert len(table) == 578
    assert table[["species", "so_feature", "block_size"]].duplicated().sum() == 0
    assert table["superinformation"].notna().all()
    assert table["n_blocks"].gt(0).all()
    digest = hashlib.sha256((REFERENCE / "superinformation_values.csv").read_bytes()).hexdigest()
    assert digest == "be3defe739306a40c89799834ab6b64ff651dd09cffae7e934c703b046139347"
    _assert_rows(table, [
        {"species": "human", "so_feature": "exon", "block_size": 20, "superinformation": 3.544893},
        {"species": "human", "so_feature": "exon", "block_size": 40, "superinformation": 2.979802},
        {"species": "human", "so_feature": "exon", "block_size": 60, "superinformation": 2.695481},
    ], ["superinformation"])


def test_retained_cis_schema_and_numerical_rows():
    si = pd.read_csv(_require_reference("superinformation_values.csv"))
    table = pd.read_csv(_require_reference("CIS_values.csv"))
    assert list(table.columns) == ["species", "so_feature", "CIS"]
    assert len(table) == 88
    for row in table.itertuples(index=False):
        group = si[(si.species == row.species) & (si.so_feature == row.so_feature)]
        calculated = compute_cis(list(zip(group.block_size, group.superinformation)), zero_initial="nan")
        if pd.isna(row.CIS):
            assert pd.isna(calculated)
        else:
            assert calculated == pytest.approx(row.CIS, abs=5e-6)


def test_retained_ips_schema_and_numerical_rows():
    si = pd.read_csv(_require_reference("superinformation_values.csv"))
    table = pd.read_csv(_require_reference("IPS_values.csv"))
    assert list(table.columns) == ["species", "so_feature", "IPS_raw", "IPS_norm"]
    assert len(table) == 88
    for row in table.itertuples(index=False):
        group = si[(si.species == row.species) & (si.so_feature == row.so_feature)]
        points = list(zip(group.block_size, group.superinformation))
        assert compute_ips(points) == pytest.approx(row.IPS_raw, abs=5e-6)
        assert compute_ips(points, normalized=True) == pytest.approx(row.IPS_norm, abs=5e-6)


def test_canonical_feature_level_pca_reference():
    table = pd.read_csv(_require_reference("superinformation_values.csv"))
    sklearn = pytest.importorskip("sklearn.decomposition")
    matrix = table.groupby(["so_feature", "block_size"])["superinformation"].mean().unstack().dropna()
    ratios = sklearn.PCA(n_components=2).fit(matrix.values).explained_variance_ratio_
    assert ratios[0] == pytest.approx(0.9334130852249626, abs=1e-12)
    assert ratios[1] == pytest.approx(0.04517565313436555, abs=1e-12)
    assert ratios.sum() == pytest.approx(0.9785887383593281, abs=1e-12)


def test_reference_si_data_contract_and_geometry_shape():
    table = pd.read_csv(_require_reference("superinformation_values.csv"))
    expected_features = {
        "CDS", "exon", "five_prime_UTR", "three_prime_UTR", "lnc_RNA",
        "ncRNA_gene", "miRNA", "snRNA", "snoRNA", "rRNA", "tRNA",
        "mRNA", "pseudogene", "transposable_element",
    }
    assert set(table.so_feature) == expected_features
    assert set(table.block_size) == {20, 40, 60, 80, 100, 150, 200}
    matrix = table.groupby(["so_feature", "block_size"]).superinformation.mean().unstack()
    assert matrix.dropna().shape == (13, 7)
    assert "tRNA" not in matrix.dropna().index
