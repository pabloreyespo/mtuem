"""Replication of the activity participation article: a full fit per article model and backend.

Each fit starts from the R start values of the fixture (tests/reference/article), runs EM as R
does (at most 100 iterations, tolerance 1e-5) and then BFGS (at most 500 iterations), and must
land on the polished R optimum: the estimates of the free parameters and their standard errors
and robust standard errors are compared with ``r_polished_*``. The article's own fits stopped
before the optimum (Apollo's bfgs stops once an iteration gains about 1e-3), so
tests/reference/make_article_polish.R restarted each of them from its estimate with a tight
stopping rule. The published values (``r_estimate``, ``r_se``) are only reported: the summary
gives the largest distance of the Python estimates from them in published standard errors. The
fits take hours, so the tests are slow (``--runslow``); run one at a time with
``-k "2015_p5 and jax"``.

Every fit appends a row to output/replication/summary.csv and saves its result in
output/replication/<run>_<backend>/ before any assertion, so an interrupted session keeps the
fits already done.
"""

import csv
import json
import re
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

RUNS = ["2015_p5", "2015-new_p5", "2022-bigcities_p4", "2022_p5"]
SE_RTOL = {"numpy": 2e-3, "jax": 5e-3}
EST_RTOL = EST_ATOL = 1e-4
COLUMNS = [
    "target", "backend", "ll_python", "ll_r_polished", "ll_r_published", "max_abs_estimate_diff_free",
    "max_rel_se_diff", "max_rel_robse_diff", "max_se_units_vs_published", "converged", "notes", "minutes",
]
CHOLESKY = re.compile(r"^chk(\d)(\d)_(\d+)$")


def normalize_cholesky_signs(values: pd.Series) -> pd.Series:
    """Cholesky terms with a positive diagonal in every class.

    A Cholesky factor is identified only up to the sign of each column: negating column j of
    class s (every ``chk<i><j>_<s>``) leaves ``L_s L_s'`` and the standard errors unchanged.
    Two fits of the same model can land on either sign, so both sides are compared after
    negating every column whose diagonal ``chk<j><j>_<s>`` is negative.
    """
    out = values.copy()
    terms = {name: m.groups() for name in values.index if (m := CHOLESKY.match(name))}
    for name, (i, j, s) in terms.items():
        if i == j and values[name] < 0:
            column = [n for n, (_, jj, ss) in terms.items() if jj == j and ss == s]
            out[column] = -values[column]
    return out


def append_summary(row: dict) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / "summary.csv"
    if path.exists():
        with path.open(encoding="utf-8") as f:
            if next(csv.reader(f), None) != COLUMNS:
                path.replace(path.with_name("summary_old_columns.csv"))
    new = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        if new:
            writer.writeheader()
        writer.writerow(row)


def compare(result: mtuem.Result, ref: dict) -> pd.DataFrame:
    """Python against the polished R optimum for the free parameters, Cholesky signs normalized.

    ``se_units_published`` is the distance from the published R estimate in published R
    standard errors.
    """
    free = [n for n in ref["r_estimate"] if n not in ref["fixed"]]

    def ref_series(key: str) -> pd.Series:
        return pd.Series(ref[key], dtype=float)

    py = normalize_cholesky_signs(result.estimates)
    table = pd.DataFrame({
        "py": py[free], "r": normalize_cholesky_signs(ref_series("r_polished_estimate"))[free],
        "py_se": result.se[free], "r_se": ref_series("r_polished_se")[free],
        "py_robse": result.robse[free], "r_robse": ref_series("r_polished_robse")[free],
        "r_published": normalize_cholesky_signs(ref_series("r_estimate"))[free],
        "r_published_se": ref_series("r_se")[free],
    })
    table["diff"] = (table["py"] - table["r"]).abs()
    table["rel_se"] = (table["py_se"] - table["r_se"]).abs() / table["r_se"].abs()
    table["rel_robse"] = (table["py_robse"] - table["r_robse"]).abs() / table["r_robse"].abs()
    table["se_units_published"] = (table["py"] - table["r_published"]).abs() / table["r_published_se"]
    return table


def summary_row(run: str, backend: str, result: mtuem.Result, ref: dict, table: pd.DataFrame) -> dict:
    return {
        "target": run, "backend": backend, "ll_python": f"{result.loglik:.6f}",
        "ll_r_polished": f"{ref['r_polished_loglik']:.6f}", "ll_r_published": f"{ref['r_loglik']:.6f}",
        "max_abs_estimate_diff_free": f"{table['diff'].max():.3e}",
        "max_rel_se_diff": f"{table['rel_se'].max():.3e}", "max_rel_robse_diff": f"{table['rel_robse'].max():.3e}",
        "max_se_units_vs_published": f"{table['se_units_published'].max():.3f}",
        "converged": result.converged,
    }


def check(table: pd.DataFrame, backend: str, report: str) -> None:
    """The replication tolerances against the polished R optimum."""
    np.testing.assert_allclose(table["py"], table["r"], rtol=EST_RTOL, atol=EST_ATOL, err_msg=report)
    np.testing.assert_allclose(table["py_se"], table["r_se"], rtol=SE_RTOL[backend], err_msg=report)
    np.testing.assert_allclose(table["py_robse"], table["r_robse"], rtol=SE_RTOL[backend], err_msg=report)


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
    append_summary(summary_row(run, backend, result, ref, table) | {"notes": "; ".join(notes), "minutes": f"{minutes:.1f}"})

    worst = table.sort_values("diff", ascending=False).head(10)
    report = (
        f"{run} {backend}: LL python {result.loglik:.6f}, LL polished R {ref['r_polished_loglik']:.6f}, "
        f"difference {result.loglik - ref['r_polished_loglik']:.6f}\n"
        f"ten largest estimate differences:\n{worst[['py', 'r', 'diff']].to_string()}"
    )
    assert result.converged, report
    check(table, backend, report)


def test_cholesky_sign_normalisation_keeps_the_covariance():
    values = pd.Series({"chk11_1": -2.0, "chk21_1": 0.5, "chk22_1": -1.5, "chk11_2": 3.0, "chk21_2": -0.7,
                        "chk22_2": 1.2, "th1_1": -0.3})
    out = normalize_cholesky_signs(values)

    def sigma(v, s):
        L = np.array([[v[f"chk11_{s}"], 0.0], [v[f"chk21_{s}"], v[f"chk22_{s}"]]])
        return L @ L.T

    assert all(out[f"chk{j}{j}_{s}"] > 0 for j in (1, 2) for s in (1, 2))
    assert out["th1_1"] == -0.3
    for s in (1, 2):
        np.testing.assert_allclose(sigma(out, s), sigma(values, s), rtol=1e-15)
    pd.testing.assert_series_equal(out[["chk11_2", "chk21_2", "chk22_2"]], values[["chk11_2", "chk21_2", "chk22_2"]])
