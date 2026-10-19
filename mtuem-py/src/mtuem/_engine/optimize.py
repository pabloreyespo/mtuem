"""Quasi-Newton optimiser for sums of per-individual log-likelihoods."""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

MESSAGES = {
    0: "successful convergence",
    8: "false convergence",
    9: "function evaluation limit reached without convergence",
    10: "iteration limit reached without convergence",
}


@dataclass
class _OptResult:
    x: np.ndarray
    f: float
    code: int
    message: str
    iterations: int
    function_evals: int


def _maximise(
    fun_vec: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    max_iterations: int = 500,
    max_function_evals: int = 100000,
    rel_f_tol: float = 1e-10,
    rel_grad_tol: float = 1e-6,
    silent: bool = True,
    exact_grad_and_scores: Callable | None = None,
) -> _OptResult:
    """Maximise a sum of per-individual log-likelihoods, ``sum(fun_vec(x))``.

    Quasi-Newton (BFGS) with a backtracking line search. The curvature model is seeded with the
    BHHH matrix (outer product of the individual scores), which the central differences used for
    the gradient deliver at no extra cost. Non-finite values are treated as an over-long step, so
    the search backs away from regions where the likelihood is zero. Code 0 is convergence (the
    gain predicted by the quadratic model is below ``rel_f_tol`` times |f|); 8, 9 and 10 report
    false convergence and the function evaluation and iteration limits.
    """
    evals = 0  # function evaluations outside the gradient

    def fvec(x):
        v = np.asarray(fun_vec(x), dtype=float)
        return v if np.all(np.isfinite(v)) else None

    def f(x):
        nonlocal evals
        evals += 1
        v = fvec(x)
        return float(v.sum()) if v is not None else -np.inf

    def grad_and_scores(x):
        n = x.size
        S = None
        for i in range(n):
            h = 1e-5 * max(1.0, abs(x[i]))
            e = np.zeros(n)
            e[i] = h
            up, down = fvec(x + e), fvec(x - e)
            if up is None or down is None:
                return None, None
            col = (up - down) / (2 * h)
            if S is None:
                S = np.zeros((col.size, n))
            S[:, i] = col
        return S.sum(axis=0), S

    if exact_grad_and_scores is not None:
        grad_and_scores = exact_grad_and_scores  # noqa: F811  (autodiff replaces finite differences)

    def bhhh_inverse(S):
        B = S.T @ S
        B += 1e-8 * np.trace(B) / B.shape[0] * np.eye(B.shape[0])
        try:
            return np.linalg.inv(B)
        except np.linalg.LinAlgError:
            return np.eye(B.shape[0])

    x = np.asarray(x0, dtype=float).copy()
    fx = f(x)
    if not np.isfinite(fx):
        raise FloatingPointError("Log-likelihood is not finite at the starting values")
    g, S = grad_and_scores(x)
    if g is None:
        raise FloatingPointError("Log-likelihood is not finite next to the starting values")
    n = x.size
    Hinv = bhhh_inverse(S)
    restarted = True
    code = 10
    it = 0
    for it in range(1, max_iterations + 1):
        direction = Hinv @ g
        if direction @ g <= 0:
            Hinv = bhhh_inverse(S)
            restarted = True
            direction = Hinv @ g
        slope = direction @ g
        step = 1.0
        accepted = False
        while step > 1e-14:
            xn = x + step * direction
            fn = f(xn)
            if np.isfinite(fn) and fn >= fx + 1e-4 * step * slope:
                accepted = True
                break
            step *= 0.5 if np.isfinite(fn) else 0.1
        gn = Sn = None
        if accepted:
            gn, Sn = grad_and_scores(xn)
            accepted = gn is not None
        if not accepted:
            if not restarted:
                Hinv = bhhh_inverse(S)
                restarted = True
                continue
            # No step improves the function: converged if the predicted gain is negligible
            # (numerical noise), otherwise false convergence.
            code = 0 if 0.5 * slope / max(abs(fx), 1.0) < 1e-7 else 8
            break
        s = xn - x
        y = -(gn - g)  # change in the gradient of -f, the function BFGS minimises
        sy = s @ y
        if sy > 1e-12 * np.linalg.norm(s) * np.linalg.norm(y):
            rho = 1.0 / sy
            V = np.eye(n) - rho * np.outer(s, y)
            Hinv = V @ Hinv @ V.T + rho * np.outer(s, s)
        restarted = False
        rel_change = abs(fn - fx) / max(abs(fn), 1.0)
        x, fx, g, S = xn, fn, gn, Sn
        pred = 0.5 * float(g @ Hinv @ g) / max(abs(fx), 1.0)
        rel_grad = float(np.max(np.abs(g) * np.maximum(np.abs(x), 1.0))) / max(abs(fx), 1.0)
        if not silent:
            print(f"  it {it:4d}   LL {fx: .6f}   predicted rel. gain {pred:.2e}   rel. gradient {rel_grad:.2e}")
        if pred < rel_f_tol and rel_change < 1e3 * rel_f_tol:
            if rel_grad < 1e2 * rel_grad_tol:
                code = 0
                break
            # The quasi-Newton matrix has lost curvature along directions where the gradient
            # is still large: re-seed it with BHHH.
            Hinv = bhhh_inverse(S)
            restarted = True
        if evals >= max_function_evals:
            code = 9
            break
    return _OptResult(x, fx, code, MESSAGES[code], it, evals)
