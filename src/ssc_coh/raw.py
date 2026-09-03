"""Raw loaders: read the 11 source CSVs without modifying them.

Every fix happens downstream in clean.py; these functions only parse.
"""
from __future__ import annotations

import pandas as pd

from .config import RAW_DIR, TABLES


def load_raw() -> dict[str, pd.DataFrame]:
    """Load all source tables keyed by filename stem. Raw bytes untouched."""
    out: dict[str, pd.DataFrame] = {}
    for stem in TABLES:
        path = RAW_DIR / f"{stem}.csv"
        out[stem] = pd.read_csv(path, low_memory=False, keep_default_na=True)
    return out
