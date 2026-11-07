"""Shared helpers for locating and reading FlyBase precomputed files."""
from __future__ import annotations

import csv
import gzip
import io
import os
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = Path(os.environ.get("DNAI_RAW_DIR", REPO_ROOT / "data" / "raw"))
PROCESSED_DIR = Path(os.environ.get("DNAI_PROCESSED_DIR", REPO_ROOT / "data" / "processed"))

_SUFFIXES = ("", ".tsv", ".tsv.gz", ".txt", ".txt.gz", ".gz")


def find_raw(name: str, raw_dir: Path | None = None) -> Path:
    """Resolve a clean file name (e.g. 'fbgn_annotation_ID') to a file in raw_dir.

    Accepts the bare name or the name with a common tabular suffix.
    """
    raw_dir = Path(raw_dir) if raw_dir is not None else RAW_DIR
    for suffix in _SUFFIXES:
        candidate = raw_dir / f"{name}{suffix}"
        if candidate.is_file():
            return candidate
    tried = ", ".join(s or "<none>" for s in _SUFFIXES)
    raise FileNotFoundError(f"Could not find '{name}' in {raw_dir} (tried suffixes: {tried})")


def open_text(path: Path):
    """Open a plain or gzip-compressed text file (detected by magic bytes)."""
    path = Path(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


def read_flybase_tsv(path: Path, columns: list[str]) -> pd.DataFrame:
    """Read a FlyBase precomputed TSV into a string DataFrame with the given column names.

    FlyBase files mix '#'/'##' comment lines (including a commented header row) with
    tab-separated data, and the first line is often a short comment. Comment and blank
    lines are dropped before parsing so pandas never infers the width from them.
    Extra trailing columns are dropped; missing ones are filled with NA.
    """
    buf = io.StringIO()
    width = len(columns)
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            line = line.rstrip("\r\n")
            width = max(width, line.count("\t") + 1)
            buf.write(line)
            buf.write("\n")
    buf.seek(0)

    df = pd.read_csv(
        buf,
        sep="\t",
        header=None,
        names=list(range(width)),
        dtype=str,
        keep_default_na=False,
        na_values=[""],
        quoting=csv.QUOTE_NONE,  # FlyBase symbols can contain quote characters
    )
    df = df.iloc[:, : len(columns)]
    df.columns = columns
    for col in columns:
        df[col] = df[col].str.strip().replace("", pd.NA)
    return df


def write_parquet(df: pd.DataFrame, name: str, out_dir: Path | None = None) -> Path:
    out_dir = Path(out_dir) if out_dir is not None else PROCESSED_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.parquet"
    df.to_parquet(path, index=False)
    return path
