"""Replication of the activity participation article: a full fit per article model and backend.

Each fit starts from the R start values of the fixture (tests/reference/article), runs EM as R
does (at most 100 iterations, tolerance 1e-5) and then BFGS (at most 500 iterations), and must
land on the R optimum: the estimates of the free parameters and their standard errors and robust
standard errors are compared with the R values. The fits take hours, so the tests are slow
(``--runslow``); run one at a time with ``-k "2015_p5 and jax"``.

Every fit appends a row to output/replication/summary.csv and saves its result in
output/replication/<run>_<backend>/ before any assertion, so an interrupted session keeps the
fits already done.
"""

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mtuem

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
ARTICLE = Path(__file__).parent / "reference" / "article"
OUTPUT = Path(__file__).resolve().parents[1] / "output" / "replication"
sys.path.insert(0, str(EXAMPLES))

import activity_participation as ap  # noqa: E402

RUNS = ["2015_p5", "2015-new_p5", "2022-bigcities_p4", "2022_p5", "2022-rest_p5"]
SE_RTOL = {"numpy": 2e-3, "jax": 5e-3}
EST_RTOL = EST_ATOL = 1e-4
COLUMNS = [
    "target", "backend", "ll_python", "ll_r", "max_abs_estimate_diff_free", "max_rel_se_diff",
    "max_rel_robse_diff", "converged", "notes", "minutes",
]


def append_summary(row: dict) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / "summary.csv"
    new = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        if new:
            writer.writeheader()
        writer.writerow(row)


def compare(result: mtuem.Result, ref: dict) -> pd.DataFrame:
    """Python against R for the free parameters of the R fit."""
    free = [n for n in ref["r_estimate"] if n not in ref["fixed"]]
    table = pd.DataFrame({
        "py": result.estimates[free], "r": pd.Series(ref["r_estimate"])[free],
        "py_se": result.se[free], "r_se": pd.Series(ref["r_se"], dtype=float)[free],
        "py_robse": result.robse[free], "r_robse": pd.Series(ref["r_robse"], dtype=float)[free],
    })
    table["diff"] = (table["py"] - table["r"]).abs()
    table["rel_se"] = (table["py_se"] - table["r_se"]).abs() / table["r_se"].abs()
    table["rel_robse"] = (table["py_robse"] - table["r_robse"]).abs() / table["r_robse"].abs()
    return table


@pytest.mark.slow
@pytest.mark.parametrize("backend", ["jax", "numpy"])
@pytest.mark.parametrize("run", RUNS)
def test_article_replication(run, backend):
    ref = json.loads((ARTICLE / f"{run}.json").read_text(encoding="utf-8"))
    model, beta, data, control = ap.from_fixture(ARTICLE / f"{run}.json")
    folder = OUTPUT / f"{run}_{backend}"
    control = control.model_copy(update={"output_dir": folder})
    settings = mtuem.EstimationSettings(
        method="em", fallback="bfgs", max_iterations=500, em_max_iterations=100, em_tolerance=1e-5,
        backend=backend,
    )

    start = time.perf_counter()
    result = model.fit(data, beta, control, settings)
    minutes = (time.perf_counter() - start) / 60

    table = compare(result, ref)
    notes = list(result.notes)
    try:
        result.save(folder)
    except Exception as exc:  # the comparison still gets recorded
        notes.append(f"save failed: {exc!r}")
    append_summary({
        "target": run, "backend": backend, "ll_python": f"{result.loglik:.6f}",
        "ll_r": f"{ref['r_loglik']:.6f}", "max_abs_estimate_diff_free": f"{table['diff'].max():.3e}",
        "max_rel_se_diff": f"{table['rel_se'].max():.3e}", "max_rel_robse_diff": f"{table['rel_robse'].max():.3e}",
        "converged": result.converged, "notes": "; ".join(notes), "minutes": f"{minutes:.1f}",
    })

    worst = table.sort_values("diff", ascending=False).head(10)
    report = (
        f"{run} {backend}: LL python {result.loglik:.6f}, LL R {ref['r_loglik']:.6f}, "
        f"difference {result.loglik - ref['r_loglik']:.6f}\n"
        f"ten largest estimate differences:\n{worst[['py', 'r', 'diff']].to_string()}"
    )
    assert result.converged, report
    np.testing.assert_allclose(table["py"], table["r"], rtol=EST_RTOL, atol=EST_ATOL, err_msg=report)
    np.testing.assert_allclose(table["py_se"], table["r_se"], rtol=SE_RTOL[backend], err_msg=report)
    np.testing.assert_allclose(table["py_robse"], table["r_robse"], rtol=SE_RTOL[backend], err_msg=report)
