# ppi-interface-transfer

## Project question

How much accuracy do structure-based GNN interface predictors (PeSTo) lose
when run on AlphaFold-predicted structures instead of experimentally
determined ones, and which structural features explain that gap?

See `docs/Project_Proposal.pdf` for the full proposal.

## Layout

- `src/` — importable Python package (analysis code, data pipeline, models)
- `data/raw/` — untouched downloaded data
- `data/interim/` — intermediate, partially processed data
- `data/processed/` — final datasets ready for modeling/analysis
- `results/` — figures, tables, metrics
- `notebooks/` — exploratory notebooks
- `tests/` — pytest test suite
- `logs/` — run logs
- `external/` — third-party code/tools (e.g. PeSTo checkout)

## Setup

A Python 3.11 virtualenv already exists at `.venv`. Install/update
dependencies with:

```bash
.venv/bin/pip install -r requirements.txt
```

Run tests with:

```bash
.venv/bin/python -m pytest
```

All paths and tunable thresholds live in `src/config.py` — see
`CLAUDE.md` for conventions.
