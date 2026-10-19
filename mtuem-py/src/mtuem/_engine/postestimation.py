"""Post-estimation: delta method standard errors and predictions."""

from collections.abc import Callable, Mapping

import numpy as np
import pandas as pd

from .. import numderiv
from .estimate import Estimate
from .inputs import Inputs
from .latent import ProbFn


def delta_method(
    est: Estimate, expressions: Mapping[str, Callable], varcov: str = "robust", silent: bool = False
) -> pd.DataFrame:
    """Delta method standard errors for functions of the estimates (follows ``apollo_deltaMethod``).

    ``expressions`` maps a label to a callable that receives the full parameter vector as a dict
    (fixed parameters included) and returns a number. ``varcov`` is ``"robust"`` (default) or
    ``"classical"``. Fixed parameters carry no variance.
    """
    vc = est.robvarcov if varcov == "robust" else est.varcov
    free = list(vc.index)
    b0 = est.estimate.copy()
    x0 = b0[free].to_numpy(dtype=float)
    rows = []
    for label, fn in expressions.items():
        def f(z, fn=fn):
            b = b0.copy()
            b[free] = z
            return float(fn(b.to_dict()))

        value = f(x0)
        g = numderiv.gradient(f, x0)
        se = float(np.sqrt(g @ vc.to_numpy() @ g))
        rows.append({"Expression": label, "Value": round(value, 4), "s.e.": round(se, 4),
                     "t-ratio (0)": round(value / se, 2) if se > 0 else float("nan")})
    out = pd.DataFrame(rows)
    if not silent:
        print(f"\nRunning Delta method computation for user-defined function using {varcov} standard errors\n")
        print(out.to_string(index=False))
    return out


def prediction(est: Estimate, prob_fn: ProbFn, inputs: Inputs, silent: bool = False) -> pd.DataFrame:
    """Model predictions at the estimates."""
    P = prob_fn(est.estimate.to_dict(), inputs, "prediction")
    pred = pd.DataFrame(P["model"]).reset_index(drop=True)
    if not silent:
        print("Running predictions from model using estimated parameters...\n")
        print("Prediction at user provided parameters")
        summ = pd.DataFrame([pred.sum(), pred.mean()], index=["Aggregate", "Average"]).round(2)
        print(summ.to_string())
        print()
    obs = inputs.database.groupby(inputs.indiv_index).cumcount().to_numpy() + 1
    pred.insert(0, "Observation", obs)
    pred.insert(0, "ID", inputs.ids)
    return pred
