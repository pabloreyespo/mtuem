"""Datasets shipped with mtuem, converted from the R package (data/*.rda) to parquet."""

from importlib.resources import files

import pandas as pd

DATASETS = ("maed", "enut.i", "enut.i.raw", "enut.ii", "enut.ii.raw", "enut.bigcities", "enut.bigcities.raw")


def load_dataset(name: str) -> pd.DataFrame:
    """Load a bundled dataset by its R name, e.g. ``load_dataset("enut.ii")``."""
    if name not in DATASETS:
        raise ValueError(f"Unknown dataset '{name}'. Available: {', '.join(DATASETS)}")
    path = files("mtuem") / "data" / f"{name.replace('.', '_')}.parquet"
    with path.open("rb") as fh:
        return pd.read_parquet(fh)


def maed() -> pd.DataFrame:
    """Subset of the Mexican Activity and Expenditure Dataset (MAED), as ``mtuem::maed``."""
    return load_dataset("maed")
