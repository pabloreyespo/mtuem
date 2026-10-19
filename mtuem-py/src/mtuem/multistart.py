"""Multi-start estimation, a port of ``mtuem::customMultiStart``."""

import contextlib
import io
from collections.abc import Mapping

import numpy as np
import pandas as pd

from . import _engine as engine
from ._engine.latent import _loglik_per_indiv
from ._engine.loglike import _model_component

DEFAULT_ESTIMATE_SETTINGS = {
    "writeIter": False,
    "maxIterations": 500,
    "hessianRoutine": "numDeriv",
    "scaleHessian": False,
    "scaleAfterConvergence": False,
    "bgw_settings": {"maxFunctionEvals": 1000},
}


def _estimate_kwargs(est_set: Mapping, silent: bool | None = None) -> dict:
    """Keyword arguments of ``engine.estimate`` from the settings dict of ``customMultiStart``."""
    kwargs = {
        "max_iterations": int(est_set["maxIterations"]),
        "max_function_evals": int(dict(est_set.get("bgw_settings") or {}).get("maxFunctionEvals", 100000)),
        "hessian": "none" if est_set.get("hessianRoutine") == "none" else None,
        "silent": bool(est_set.get("silent", False)) if silent is None else silent,
    }
    if "saddleEscapes" in est_set:
        kwargs["saddle_escapes"] = int(est_set["saddleEscapes"])
    return kwargs


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


def customMultiStart(
    apollo_beta: Mapping[str, float],
    apollo_fixed: list[str],
    apollo_probabilities,
    apollo_inputs,
    customMultistart_settings: Mapping | None = None,
    estimate_settings: Mapping | None = None,
    estimation_settings: Mapping | None = None,
    first_em: bool = False,
    em_iter_max: int = 5,
    non_saddle: bool = True,
    verbose: bool = False,
    normalization: Mapping | None = None,
    nClass: int = 1,
    seed: int | None = None,
) -> dict:
    """Estimate from several random starting points and keep the best non-saddle optimum.

    Returns ``best_model`` (None when no candidate qualifies; R returns NA), ``best_ll`` and
    ``models`` (a dict from candidate number to model, None for failed candidates).
    """
    beta = pd.Series(apollo_beta, dtype=float)
    fixed = list(apollo_fixed or [])
    ms = {"apolloBetaMax": beta + 1, "apolloBetaMin": beta - 1, "nCandidates": 100}
    ms.update(customMultistart_settings or {})
    if estimate_settings is None and estimation_settings is not None:
        estimate_settings = estimation_settings  # deprecated alias, kept as in R
    est_set = dict(DEFAULT_ESTIMATE_SETTINGS)
    est_set.update(estimate_settings or {})
    est_set["scaleAfterConvergence"] = False
    rng = np.random.default_rng(seed)

    n_cand = int(ms["nCandidates"])
    beta_max = pd.Series(ms["apolloBetaMax"], dtype=float)
    beta_min = pd.Series(ms["apolloBetaMin"], dtype=float)
    variables = [v for v in beta.index if v not in fixed]
    beta_matrix = pd.DataFrame(np.tile(beta.to_numpy(), (n_cand, 1)), columns=beta.index)
    for v in variables:
        beta_matrix[v] = rng.uniform(beta_min[v], beta_max[v], n_cand)

    if normalization is not None and normalization.get("normalization") in ("alpha_beta", "theta_phi"):
        norm = normalization["normalization"]
        keys = ("work_elasticities", "times_elasticities", "goods_elasticities", "sig", "rho", "cholesky")
        params = {k: list(normalization.get(k) or []) for k in keys}
        is_lc = apollo_inputs.lc_pars is not None
        if not is_lc and nClass <= 1:
            beta_matrix = _normalize(beta_matrix, rng, variables, norm, *params.values())
        else:
            for k in range(1, nClass + 1):
                suffixed = [[f"{p}_{k}" for p in v] for v in params.values()]
                beta_matrix = _normalize(beta_matrix, rng, variables, norm, *suffixed)
            # Parameters shared across classes are passed without a suffix.
            beta_matrix = _normalize(beta_matrix, rng, variables, norm, *params.values())

    print("Testing candidates...", end="")
    works = []
    for i in range(n_cand):
        with np.errstate(all="ignore"):
            P = apollo_probabilities(beta_matrix.iloc[i].to_dict(), apollo_inputs, "estimate")
            ll = float(np.sum(_loglik_per_indiv(_model_component(P), apollo_inputs)))
        works.append(not np.isnan(ll))
    print("Done...", sum(works), "candidates available")

    models: dict[int, engine.Estimate | None] = {}
    best_model, best_ll = None, -np.inf
    kwargs = _estimate_kwargs(est_set)
    quiet = _estimate_kwargs(est_set, silent=True)
    for i in [i for i in range(n_cand) if works[i]]:
        cand = i + 1
        print(f"Candidate {cand} ...", end="" if not verbose else "\n")
        start = beta_matrix.iloc[i].to_dict()
        try:
            sink = contextlib.nullcontext() if verbose else contextlib.redirect_stdout(io.StringIO())
            if first_em:
                print("Starting EM...", end="" if not verbose else "\n")
                with sink:
                    em = engine.lc_em(start, fixed, apollo_probabilities, apollo_inputs,
                                      max_iterations=em_iter_max, post_em=0, estimate_kwargs=quiet)
                start = em.estimate.to_dict()
                print(f"LL: {-round(em.maximum, 2)} ...", end="")
            print("Starting solver...", end="" if not verbose else "\n")
            sink = contextlib.nullcontext() if verbose else contextlib.redirect_stdout(io.StringIO())
            with sink, np.errstate(all="ignore"):
                model = engine.estimate(start, fixed, apollo_probabilities, apollo_inputs, **(kwargs if verbose else quiet))
            models[cand] = model
            print(f"LL: {-round(model.maximum, 2)} ...Done...", end="")
        except (FloatingPointError, ValueError, np.linalg.LinAlgError, ZeroDivisionError) as err:
            print(f"Failed ({err})")
            models[cand] = None
            continue

        converged = model.code == 0
        eig = model.eigValue if np.isfinite(model.eigValue) else np.inf
        ok = converged and model.maximum > best_ll and (eig < 0 if non_saddle else True)
        if ok:
            best_ll, best_model = model.maximum, model
            print("Best found", end="")
        print()

    return {"best_model": best_model, "best_ll": best_ll, "models": models}
