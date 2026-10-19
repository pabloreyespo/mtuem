"""Convert the R package datasets (data/*.rda) to parquet files bundled with mtuem-py.

Run from mtuem-py/ with:  uv run python data-raw/convert_data.py
"""

from pathlib import Path

import pyreadr

SRC = Path(__file__).resolve().parents[2] / "data"
DST = Path(__file__).resolve().parents[1] / "src" / "mtuem" / "data"

NAMES = ["maed", "enut.i", "enut.i.raw", "enut.ii", "enut.ii.raw", "enut.bigcities", "enut.bigcities.raw"]


def main() -> None:
    DST.mkdir(parents=True, exist_ok=True)
    for name in NAMES:
        frames = pyreadr.read_r(str(SRC / f"{name}.rda"))
        df = frames[name]
        out = DST / f"{name.replace('.', '_')}.parquet"
        df.to_parquet(out, index=False, compression="gzip")
        print(f"{name}: {df.shape} -> {out.name}")


if __name__ == "__main__":
    main()
