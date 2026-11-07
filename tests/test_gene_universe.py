import json
from pathlib import Path

import pandas as pd
import pytest

from src.flybase_io import find_raw, read_flybase_tsv
from src.gene_universe import run

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="module")
def outputs(tmp_path_factory):
    out = tmp_path_factory.mktemp("processed")
    summary = run(FIXTURES, out)
    return {
        "summary": summary,
        "universe": pd.read_parquet(out / "gene_universe.parquet"),
        "fbpp": pd.read_parquet(out / "fbgn_fbpp.parquet"),
        "id_map": pd.read_parquet(out / "fbgn_id_map.parquet"),
        "summary_file": json.loads((out / "gene_universe_summary.json").read_text()),
    }


def test_find_raw_resolves_suffixes():
    assert find_raw("fbgn_annotation_ID", FIXTURES).name == "fbgn_annotation_ID.tsv"
    assert find_raw("gene_map_table", FIXTURES).name == "gene_map_table.tsv.gz"
    with pytest.raises(FileNotFoundError):
        find_raw("does_not_exist", FIXTURES)


def test_reader_skips_comments_and_short_first_line():
    df = read_flybase_tsv(FIXTURES / "fbgn_annotation_ID.tsv", ["symbol", "org", "fbgn", "sec", "cg", "sec_cg"])
    assert not df["symbol"].str.startswith("#").any()
    assert len(df) == 5
    assert df.loc[df["symbol"] == "w", "sec"].isna().all()


def test_universe_membership(outputs):
    u = outputs["universe"]
    # Act5C and w are current, protein-coding, with polypeptides.
    # mir-1 is ncRNA; geneC has no FBpp; newgene is absent from fbgn_annotation_ID; Dsim row is not Dmel.
    assert list(u["fbgn"]) == ["FBgn0000042", "FBgn0003996"]
    assert u["fbgn"].is_unique


def test_universe_columns_and_isoforms(outputs):
    u = outputs["universe"].set_index("fbgn")
    assert u.loc["FBgn0000042", "symbol"] == "Act5C"
    assert u.loc["FBgn0000042", "annotation_id"] == "CG4027"
    assert u.loc["FBgn0000042", "n_transcripts"] == 2
    assert u.loc["FBgn0000042", "n_polypeptides"] == 2
    assert u.loc["FBgn0000042", "chrom_arm"] == "X"
    assert u.loc["FBgn0000042", "seq_start"] == 5900861
    assert u.loc["FBgn0000042", "strand"] == -1


def test_fbpp_map_is_deduplicated(outputs):
    fbpp = outputs["fbpp"]
    assert len(fbpp) == 3
    assert set(fbpp["fbgn"]) == {"FBgn0000042", "FBgn0003996"}
    assert not fbpp.duplicated().any()


def test_id_map(outputs):
    m = outputs["id_map"].set_index("from_fbgn")["fbgn"]
    assert m["FBgn0000042"] == "FBgn0000042"
    assert m["FBgn0000043"] == "FBgn0000042"
    # FBgn0099999 is a secondary ID of two different genes -> ambiguous, dropped.
    assert "FBgn0099999" not in m.index
    assert "FBgn0200001" not in m.index  # Dsim


def test_summary(outputs):
    s = outputs["summary"]
    assert s == outputs["summary_file"]
    assert s["universe_genes"] == 2
    assert s["protein_coding_not_in_annotation_id"] == 1
    assert s["protein_coding_without_fbpp"] == 1
