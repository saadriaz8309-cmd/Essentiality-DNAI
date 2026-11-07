# DNAI

A literature-unbiased gene-essentiality discovery pipeline for *Drosophila melanogaster*.
It ranks **candidate essential genes**: genes the model confidently predicts essential but
that FlyBase marks viable or has never annotated. It uses only annotation-independent
features and runs on CPU only.

See [`docs/DESIGN.md`](docs/DESIGN.md) for the locked design decisions and raw-file parsing notes.

## Layout

```
data/raw/        downloaded FlyBase / STRING / eLife files (clean names, not committed)
data/processed/  parquet outputs per stage (not committed)
src/             pipeline stages
tests/           unit tests on small FlyBase-format fixtures
results/         evaluation outputs
app/             Streamlit candidate browser (Stage 7)
versions.txt     FB2026_03 / r6.69 / STRING v12.0
```

## Setup

CPU-only. On Windows, run everything in WSL2.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Put the raw files in `data/raw/` under the clean names in `docs/DESIGN.md`.
Suffixes `.tsv`, `.tsv.gz`, `.txt`, `.gz` or none are all accepted.

## Stages

| Stage | Module | Status |
|---|---|---|
| 1. Gene universe | `src/gene_universe.py` | done |
| 2. Labels (PU) | `src/labels.py` | todo |
| 3. Feature matrix | `src/features.py` | todo |
| 4. Research-intensity covariate | `src/bias.py` | todo |
| 5. PU LightGBM + grouped CV | `src/model.py` | todo |
| 6. Evaluation | `src/evaluate.py` | todo |
| 7. Candidates + app | `src/candidates.py`, `app/` | todo |

### Stage 1: gene universe

```bash
python -m src.gene_universe
```

Writes `data/processed/gene_universe.parquet`, `fbgn_fbpp.parquet`, `fbgn_id_map.parquet`
and `gene_universe_summary.json`. Check the summary: `universe_genes` should be roughly
13.9k protein-coding genes for r6.69.

## Tests

```bash
python -m pytest -q
```
