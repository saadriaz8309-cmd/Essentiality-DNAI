"""Stage 1: canonical protein-coding gene universe, keyed on primary FBgn.

Inputs (data/raw/):
  fbgn_annotation_ID       ID backbone: primary FBgn, symbol, CG annotation ID, secondary FBgns
  fbgn_fbtr_fbpp_expanded  gene_type + FBgn -> FBtr -> FBpp
  gene_map_table           (optional) cytogenetic + sequence location

Outputs (data/processed/):
  gene_universe.parquet          one row per protein-coding D. melanogaster gene
  fbgn_fbpp.parquet              FBgn, FBtr, FBpp for every polypeptide of a universe gene
  fbgn_id_map.parquet            any FBgn (primary or secondary) -> primary FBgn, unambiguous only
  gene_universe_summary.json     counts for sanity checks

A gene is in the universe iff it is a current Dmel gene in fbgn_annotation_ID, has
gene_type == "protein_coding_gene" in fbgn_fbtr_fbpp_expanded, and has at least one FBpp.

Note: n_transcripts / n_polypeptides are kept as metadata. Isoform counts partly reflect
cDNA/study depth, so treat them as potentially study-biased before using them as features.

Run:  python -m src.gene_universe
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from .flybase_io import PROCESSED_DIR, RAW_DIR, find_raw, read_flybase_tsv, write_parquet

ANNOTATION_COLS = [
    "symbol",
    "organism",
    "primary_fbgn",
    "secondary_fbgns",
    "annotation_id",
    "secondary_annotation_ids",
]
EXPANDED_COLS = [
    "organism",
    "gene_type",
    "fbgn",
    "gene_symbol",
    "gene_fullname",
    "annotation_id",
    "transcript_type",
    "fbtr",
    "transcript_symbol",
    "fbpp",
    "polypeptide_symbol",
]
GENE_MAP_COLS = [
    "organism",
    "symbol",
    "fbgn",
    "recombination_loc",
    "cytogenetic_loc",
    "sequence_loc",
]

PROTEIN_CODING = "protein_coding_gene"
FBGN_RE = re.compile(r"^FBgn\d{7}$")
FBPP_RE = re.compile(r"^FBpp\d{7}$")
# e.g. "2L:7529..9484(1)" or "X:19961297..19969323(-1)"
SEQ_LOC_RE = re.compile(r"^(?P<arm>[^:]+):(?P<start>\d+)\.\.(?P<end>\d+)(?:\((?P<strand>-?1)\))?$")


def load_annotation_ids(path: Path) -> pd.DataFrame:
    df = read_flybase_tsv(path, ANNOTATION_COLS)
    df = df[(df["organism"] == "Dmel") & df["primary_fbgn"].fillna("").str.match(FBGN_RE)]
    return df.drop_duplicates("primary_fbgn").reset_index(drop=True)


def load_expanded(path: Path) -> pd.DataFrame:
    df = read_flybase_tsv(path, EXPANDED_COLS)
    return df[(df["organism"] == "Dmel") & df["fbgn"].fillna("").str.match(FBGN_RE)].reset_index(drop=True)


def load_gene_map(path: Path) -> pd.DataFrame:
    df = read_flybase_tsv(path, GENE_MAP_COLS)
    df = df[(df["organism"] == "Dmel") & df["fbgn"].fillna("").str.match(FBGN_RE)]
    df = df.drop_duplicates("fbgn")
    loc = df["sequence_loc"].fillna("").str.extract(SEQ_LOC_RE)
    out = pd.DataFrame(
        {
            "fbgn": df["fbgn"].values,
            "cytogenetic_loc": df["cytogenetic_loc"].values,
            "chrom_arm": loc["arm"].values,
            "seq_start": pd.to_numeric(loc["start"], errors="coerce").astype("Int64").values,
            "seq_end": pd.to_numeric(loc["end"], errors="coerce").astype("Int64").values,
            "strand": pd.to_numeric(loc["strand"], errors="coerce").astype("Int8").values,
        }
    )
    return out


def build_id_map(ann: pd.DataFrame) -> pd.DataFrame:
    """Map every primary and secondary FBgn to its current primary FBgn.

    Secondary IDs that point to more than one primary (merged/split genes) are dropped
    as ambiguous. A primary ID always maps to itself.
    """
    primary = pd.DataFrame({"from_fbgn": ann["primary_fbgn"], "fbgn": ann["primary_fbgn"], "is_primary": True})
    sec = ann[["primary_fbgn", "secondary_fbgns"]].dropna()
    sec = sec.assign(from_fbgn=sec["secondary_fbgns"].str.split(",")).explode("from_fbgn")
    sec["from_fbgn"] = sec["from_fbgn"].str.strip()
    sec = sec[sec["from_fbgn"].str.match(FBGN_RE)]
    sec = sec[~sec["from_fbgn"].isin(primary["from_fbgn"])]  # never remap a current primary ID
    sec = sec.drop_duplicates(["from_fbgn", "primary_fbgn"])
    ambiguous = sec["from_fbgn"].duplicated(keep=False)
    sec = sec[~ambiguous]
    secondary = pd.DataFrame({"from_fbgn": sec["from_fbgn"], "fbgn": sec["primary_fbgn"], "is_primary": False})
    return pd.concat([primary, secondary], ignore_index=True)


def build_gene_universe(
    ann: pd.DataFrame, expanded: pd.DataFrame, gene_map: pd.DataFrame | None = None
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    current = set(ann["primary_fbgn"])

    pc = expanded[expanded["gene_type"] == PROTEIN_CODING]
    pc_with_pp = pc[pc["fbpp"].fillna("").str.match(FBPP_RE)]

    pc_genes = set(pc["fbgn"])
    keep = current & set(pc_with_pp["fbgn"])

    fbpp = (
        pc_with_pp[pc_with_pp["fbgn"].isin(keep)][["fbgn", "fbtr", "fbpp"]]
        .drop_duplicates()
        .sort_values(["fbgn", "fbpp"])
        .reset_index(drop=True)
    )
    iso = fbpp.groupby("fbgn").agg(n_transcripts=("fbtr", "nunique"), n_polypeptides=("fbpp", "nunique"))

    universe = (
        ann[ann["primary_fbgn"].isin(keep)][["primary_fbgn", "symbol", "annotation_id"]]
        .rename(columns={"primary_fbgn": "fbgn"})
        .merge(iso, left_on="fbgn", right_index=True, how="left")
    )
    if gene_map is not None:
        universe = universe.merge(gene_map, on="fbgn", how="left")
    universe = universe.sort_values("fbgn").reset_index(drop=True)

    if universe["fbgn"].duplicated().any():
        raise AssertionError("gene universe has duplicate FBgn keys")

    summary = {
        "annotation_id_current_dmel_genes": len(current),
        "expanded_protein_coding_genes": len(pc_genes),
        "protein_coding_not_in_annotation_id": len(pc_genes - current),
        "protein_coding_without_fbpp": len(pc_genes - set(pc_with_pp["fbgn"])),
        "universe_genes": len(universe),
        "universe_polypeptides": int(fbpp["fbpp"].nunique()),
        "universe_transcripts": int(fbpp["fbtr"].nunique()),
        "universe_missing_cg_id": int(universe["annotation_id"].isna().sum()),
    }
    if gene_map is not None:
        summary["universe_missing_sequence_loc"] = int(universe["chrom_arm"].isna().sum())
        summary["universe_by_chrom_arm"] = {
            str(k): int(v) for k, v in universe["chrom_arm"].value_counts(dropna=False).sort_index().items()
        }
    return universe, fbpp, summary


def run(raw_dir: Path = RAW_DIR, out_dir: Path = PROCESSED_DIR) -> dict:
    ann = load_annotation_ids(find_raw("fbgn_annotation_ID", raw_dir))
    expanded = load_expanded(find_raw("fbgn_fbtr_fbpp_expanded", raw_dir))
    try:
        gene_map = load_gene_map(find_raw("gene_map_table", raw_dir))
    except FileNotFoundError:
        gene_map = None

    universe, fbpp, summary = build_gene_universe(ann, expanded, gene_map)
    id_map = build_id_map(ann)
    summary["id_map_rows"] = len(id_map)
    summary["id_map_secondary_ids"] = int((~id_map["is_primary"]).sum())

    write_parquet(universe, "gene_universe", out_dir)
    write_parquet(fbpp, "fbgn_fbpp", out_dir)
    write_parquet(id_map, "fbgn_id_map", out_dir)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "gene_universe_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--raw-dir", type=lambda s: Path(s).expanduser(), default=RAW_DIR)
    parser.add_argument("--out-dir", type=Path, default=PROCESSED_DIR)
    args = parser.parse_args()
    summary = run(args.raw_dir, args.out_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
