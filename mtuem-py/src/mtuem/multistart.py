"""Estimation runs of ``fit``: one run (method and fallback) and the multistart.

The multistart follows ``mtuem::customMultiStart`` (R): uniform draws around the start values,
an optional short EM warm-up, the Theta-Phi and alpha-beta normalizations, and the best
converged non-saddle optimum.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from . import _engine as engine
from ._engine.latent import _loglik_per_indiv
from ._engine.loglike import _model_component
from .results import Result

if TYPE_CHECKING:
    from .beta import Beta
    from .config import Control, EstimationSettings, MultiStart

# Errors that mean an estimation failed numerically (as R's tryCatch in customMultiStart).
_FAILURES = (FloatingPointError, ValueError, np.linalg.LinAlgError, ZeroDivisionError)
_NORMALIZATION_KEYS = ("work_elasticities", "times_elasticities", "goods_elasticities", "sig", "rho", "cholesky")


def estimate_once(
    start: dict, beta: Beta, prob_fn, inputs, settings: EstimationSettings, silent: bool
) -> tuple[engine.Estimate, list[str]]:
    """One estimation from ``start`` by ``settings.method``; returns the estimate and its notes.

    ``beta`` gives the fixed parameters and the bounds. ``method="em"`` runs EM then BFGS. When
    EM raises a numerical error or does not converge, ``fallback="bfgs"`` estimates by BFGS from
    ``start`` and records why in a note; without a fallback the error propagates and a
    non-converged result is returned with a note.
    """
    fixed, bounds = beta.fixed_names, beta.bounds()
    kwargs = dict(
        max_iterations=settings.max_iterations,
        max_function_evals=settings.max_function_evals,
        hessian=settings.resolved_hessian(),
        saddle_escapes=settings.saddle_escapes,
        silent=silent,
    )
    if settings.method == "bfgs":
        return engine.estimate(start, fixed, prob_fn, inputs, bounds=bounds, **kwargs), []
    try:
        est = engine.lc_em(
            start, fixed, prob_fn, inputs, max_iterations=settings.em_max_iterations,
            tolerance=settings.em_tolerance, post_em=2, estimate_kwargs=kwargs, bounds=bounds,
        )
    except (FloatingPointError, ValueError, np.linalg.LinAlgError) as err:
        if settings.fallback is None:
            raise
        reason = str(err) or type(err).__name__
    else:
        if est.successfulEstimation:
            return est, []
        if settings.fallback is None:
            return est, [f"EM did not converge ({est.message}); no fallback"]
        reason = f"not converged: {est.message}"
    if not silent:
        print(f"EM failed ({reason}); estimating by BFGS from the start values.")
    est = engine.estimate(start, fixed, prob_fn, inputs, bounds=bounds, **kwargs)
    return est, [f"EM failed ({reason}); estimated by BFGS from the start values"]


# ----------------------------------------------------------------------------------------------
# Multistart
# ----------------------------------------------------------------------------------------------
def _cholesky_split(cholesky: list[str]) -> tuple[list[str], list[str]]:
    """Diagonal and off-diagonal names of a row-wise lower-triangular Cholesky vector."""
    diag_idx, k, idx = [], 1, 1
    while idx <= len(cholesky):
        diag_idx.append(idx - 1)
        k += 1
        idx += k
    diag = [cholesky[i] for i in diag_idx]
    off = [c for i, c in enumerate(cholesky) if i not in diag_idx]
    return diag, off


def _normalize(beta_matrix: pd.DataFrame, rng, variables: list[str], norm: str, work, times, goods, sig, rho, chol):
    """Push sampled candidates into the region where the MTUEM is well defined."""
    work, times, goods = list(work or []), list(times or []), list(goods or [])
    sig, rho, chol = list(sig or []), list(rho or []), list(chol or [])
    chol_diag, chol_off = _cholesky_split(chol) if chol else ([], [])
    n = len(beta_matrix)

    if norm == "alpha_beta":
        cols = [v for v in variables if v in work]
        sub = beta_matrix[cols]
        beta_matrix[cols] = sub.mask(sub >= 0.5, 0.49)
        cols = [v for v in variables if v in times + goods + sig + chol_diag]
        sub = beta_matrix[cols]
        beta_matrix[cols] = sub.mask(sub <= 0, 0.001)
        al = [x for x in work if x.startswith("alpha")]
        be = [x for x in work if x.startswith("beta")]
        cap_times = (1 - 2 * beta_matrix[be[0]]).to_numpy() if be else np.ones(n)
        cap_goods = (1 - 2 * beta_matrix[al[0]]).to_numpy() if al else np.ones(n)
    else:
        cols = [v for v in variables if v in work + times + goods + sig + chol_diag and not v.startswith("thw")]
        sub = beta_matrix[cols]
        beta_matrix[cols] = sub.mask(sub <= 0, 0.001)
        th = [x for x in work if x.startswith("Theta")]
        ph = [x for x in work if x.startswith("Phi")]
        cap_times = beta_matrix[th[0]].to_numpy() if len(th) == 1 and th[0] in variables else np.ones(n)
        cap_goods = beta_matrix[ph[0]].to_numpy() if len(ph) == 1 and ph[0] in variables else np.ones(n)

    for group, cap in ((times, cap_times), (goods, cap_goods)):
        group = [g for g in group if g in variables]
        if not group:
            continue
        rsum = beta_matrix[group].sum(axis=1).to_numpy()
        mask = rsum > cap
        rand = rng.uniform(0, 1, mask.sum())
        factor = cap[mask] / (rsum[mask] + rand)
        beta_matrix.loc[mask, group] = beta_matrix.loc[mask, group].to_numpy() * factor[:, None]

    zero = [v for v in rho + chol_off if v in variables]
    if zero:
        beta_matrix[zero] = 0.0
    return beta_matrix


def _box(beta: Beta, ms: MultiStart) -> tuple[pd.Series, pd.Series]:
    """Lower and upper ends of the draws of the free parameters: the start value minus and plus 1
    unless ``ms.lower`` or ``ms.upper`` names the parameter, kept inside the bounds of ``beta``."""
    free = beta.free_names
    values = pd.Series(beta.values(), dtype=float)[free]
    ends = {"lower": values - 1, "upper": values + 1}
    for side, given in (("lower", ms.lower), ("upper", ms.upper)):
        given = dict(given or {})
        unknown = [n for n in given if n not in beta]
        if unknown:
            raise ValueError(f"MultiStart.{side} names unknown parameters: {', '.join(unknown)}")
        fixed = [n for n in given if n not in free]
        if fixed:
            raise ValueError(
                f"MultiStart.{side} names fixed parameters: {', '.join(fixed)}; only free ones are drawn"
            )
        for n, v in given.items():
            ends[side][n] = float(v)
    lower, upper = ends["lower"], ends["upper"]
    for n in free:
        lo, hi = beta[n].lower, beta[n].upper
        if lo is not None:
            lower[n] = max(lower[n], lo)
        if hi is not None:
            upper[n] = min(upper[n], hi)
        if lower[n] > upper[n]:
            raise ValueError(f"the multistart range of {n} is empty: [{lower[n]}, {upper[n]}], bounds included")
    return lower, upper


def _normalization_sets(model) -> list[dict[str, list[str]]]:
    """The groups of ``model.normalization_groups()`` that ``_normalize`` takes together.

    One set for a single-class model. A latent class model gives one set per class (the names
    with the suffix ``_s``), then one of the shared names, as ``customMultiStart`` does with
    ``nClass``.
    """
    groups = {k: list(model.normalization_groups().get(k, [])) for k in _NORMALIZATION_KEYS}
    n_classes = int(getattr(model, "n_classes", 0) or 0)
    if n_classes < 2:
        return [groups]
    suffixes = tuple(f"_{s}" for s in range(1, n_classes + 1))
    sets = [{k: [n for n in v if n.endswith(sfx)] for k, v in groups.items()} for sfx in suffixes]
    sets.append({k: [n for n in v if not n.endswith(suffixes)] for k, v in groups.items()})
    return sets


def _candidates(model, beta: Beta, ms: MultiStart) -> pd.DataFrame:
    """The ``ms.n`` starting points, one row each over every parameter of ``beta``.

    Only free parameters are drawn (uniform over the range of ``_box``) and normalized; fixed
    and locked ones keep their values. The draws use ``ms.seed``.
    """
    rng = np.random.default_rng(ms.seed)
    values = pd.Series(beta.values(), dtype=float)
    free = beta.free_names
    lower, upper = _box(beta, ms)
    draws = pd.DataFrame(np.tile(values.to_numpy(), (ms.n, 1)), columns=values.index)
    for v in free:
        draws[v] = rng.uniform(lower[v], upper[v], ms.n)
    if ms.normalization is not None:
        for groups in _normalization_sets(model):
            draws = _normalize(draws, rng, free, ms.normalization, *(groups[k] for k in _NORMALIZATION_KEYS))
    return draws


def _usable(start: dict, beta: Beta, prob_fn, inputs) -> bool:
    """True when ``start`` is strictly inside the bounds and its log-likelihood is finite."""
    for n in beta.free_names:
        lo, hi = beta[n].lower, beta[n].upper
        if (lo is not None and not start[n] > lo) or (hi is not None and not start[n] < hi):
            return False
    try:
        with np.errstate(all="ignore"):
            P = prob_fn(start, inputs, "estimate")
            ll = float(np.sum(_loglik_per_indiv(_model_component(P), inputs)))
    except _FAILURES:
        return False
    return bool(np.isfinite(ll))


def run_multistart(model, beta: Beta, control: Control, settings: EstimationSettings, inputs, prob_fn) -> Result:
    """Estimate from ``settings.multistart.n`` random starts and keep the best.

    Every candidate with a finite starting log-likelihood is estimated as ``fit`` estimates
    (``settings.method`` and ``fallback``), after ``em_iter_max`` EM iterations when
    ``first_em``. The best converged result, at a maximum when ``non_saddle`` (largest Hessian
    eigenvalue below 0), is returned with ``starts``: the ``Result`` of every candidate, None
    for those skipped or failed. When no candidate qualifies, raises ``RuntimeError`` with the
    attempts as ``.attempts``.
    """
    ms = settings.multistart
    verbose = settings.verbose
    if ms.non_saddle and settings.resolved_hessian() == "none":
        raise ValueError("MultiStart.non_saddle needs the Hessian: set hessian, or non_saddle=False")
    draws = _candidates(model, beta, ms)
    fixed = beta.fixed_names
    if verbose:
        print("Testing candidates...", end="")
    usable = [_usable(draws.iloc[i].to_dict(), beta, prob_fn, inputs) for i in range(ms.n)]
    if verbose:
        print(f" Done. {sum(usable)} candidates available")

    starts: list[Result | None] = [None] * ms.n
    best, best_ll, last_err = None, -np.inf, None
    for i in [i for i in range(ms.n) if usable[i]]:
        start = draws.iloc[i].to_dict()
        start_beta = beta.model_copy(deep=True).set(**{n: start[n] for n in beta.free_names})
        try:
            with np.errstate(all="ignore"):
                if ms.first_em:
                    em = engine.lc_em(
                        start, fixed, prob_fn, inputs, max_iterations=ms.em_iter_max,
                        tolerance=settings.em_tolerance, post_em=0, estimate_kwargs={"silent": True},
                        bounds=beta.bounds(),
                    )
                    start = em.estimate.to_dict()
                est, notes = estimate_once(start, beta, prob_fn, inputs, settings, silent=True)
        except _FAILURES as err:
            last_err = err
            if verbose:
                print(f"Candidate {i + 1}: failed ({err})")
            continue
        starts[i] = Result(model, start_beta, control, settings, est, inputs, prob_fn, notes=notes)
        eig = est.eigValue if np.isfinite(est.eigValue) else np.inf
        saddle = ms.non_saddle and not eig < 0
        if est.code == 0 and not saddle and est.maximum > best_ll:
            best, best_ll = i, est.maximum
        if verbose:
            status = "best so far" if best == i else "not converged" if est.code != 0 else ""
            status = status or ("saddle point" if saddle else "")
            print(f"Candidate {i + 1}: LL {est.maximum:.4f}" + (f" ({status})" if status else ""))

    if best is None:
        failed = sum(r is None for r in starts)
        err = RuntimeError(
            f"no multistart candidate converged to a maximum ({ms.n} candidates, {sum(usable)} with a finite "
            f"starting log-likelihood, {failed} skipped or failed); the attempts are in .attempts"
            + (f"; last failure: {type(last_err).__name__}: {last_err}" if last_err is not None else "")
        )
        err.attempts = starts
        raise err from last_err
    chosen = starts[best]
    notes = chosen.notes + [f"Best of {ms.n} multistart candidates: candidate {best + 1}"]
    return Result(model, chosen.beta, control, settings, chosen._estimate, inputs, prob_fn, starts=starts, notes=notes)
