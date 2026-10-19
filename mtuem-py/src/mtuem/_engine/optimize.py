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


def _central_scores(fvec: Callable[[np.ndarray], np.ndarray | None], x: np.ndarray, steps: np.ndarray):
    """Gradient and individual scores by central differences; (None, None) if a value is not finite."""
    n = x.size
    S = None
    for i in range(n):
        h = steps[i]
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
        return _central_scores(fvec, x, 1e-5 * np.maximum(1.0, np.abs(x)))

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


def _inside(x: np.ndarray, lower: np.ndarray, upper: np.ndarray, rel: float) -> np.ndarray:
    """``x`` moved at least ``rel * max(1, |bound|)`` (at most a quarter of the width) inside its bounds."""
    width = upper - lower
    with np.errstate(invalid="ignore"):
        d_lo = np.minimum(rel * np.maximum(1.0, np.abs(lower)), 0.25 * width)
        d_hi = np.minimum(rel * np.maximum(1.0, np.abs(upper)), 0.25 * width)
    x = np.where(np.isfinite(lower), np.maximum(x, lower + np.where(np.isfinite(lower), d_lo, 0.0)), x)
    return np.where(np.isfinite(upper), np.minimum(x, upper - np.where(np.isfinite(upper), d_hi, 0.0)), x)


def maximise_bounded(
    fun_vec: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    exact_grad_and_scores: Callable | None = None,
    mu: float = 1e-4,
    outer_eps: float = 1e-5,
    outer_iterations: int = 100,
    **maximise_kwargs,
) -> _OptResult:
    """Maximise ``sum(fun_vec(x))`` subject to ``lower < x < upper`` (infinite entries are no bound).

    The adaptive logarithmic barrier of R's ``constrOptim``, which maxLik uses for Apollo's
    inequality constraints. With ``g(x)`` the stacked distances ``x - lower`` and ``upper - x``
    for the finite bounds and ``A = dg/dx``, each outer iteration maximises

        R(x) = f(x) + mu * sum(g(x_k) * log g(x) - A x)

    from ``x_k`` with ``_maximise``. The barrier gradient ``mu * A'(g(x_k) / g(x) - 1)`` vanishes
    at ``x_k``, so the penalty fades as the iterates settle and the limit is the constrained
    maximum, also when it lies on a bound. The loop stops when an outer iteration improves R
    by less than ``tol = outer_eps * (|R| + outer_eps)``, or when the gain still available by
    closing the gaps to the bounds, ``mu * sum|g(x_k) - g(x_{k+1})|``, is below ``tol``. The
    returned ``f`` is the unpenalised objective, and ``x`` is kept at least 1e-12 (relative)
    inside the bounds.

    The barrier is used only when the unconstrained maximum (from ``x0``) is not strictly
    inside the bounds; otherwise that maximum is returned as is. The barrier enters
    ``_maximise`` as one extra individual whose value and gradient are exact; its curvature
    ``A' diag(mu g(x_k) / g^2) A`` is added to the BHHH seed. The scores of ``fun_vec`` are
    central differences (with steps kept inside the bounds) unless ``exact_grad_and_scores`` is
    given. ``max_iterations`` is shared by the unconstrained run and the inner barrier runs.
    The code is that of the last inner run: an inner run that does not converge (false
    convergence, or the iteration or evaluation limit) ends the loop with its code, and a budget
    used up before a barrier run gives code 10.
    """
    x = np.asarray(x0, dtype=float).copy()
    lower = np.broadcast_to(np.asarray(lower, dtype=float), x.shape)
    upper = np.broadcast_to(np.asarray(upper, dtype=float), x.shape)
    lo, hi = np.isfinite(lower), np.isfinite(upper)
    if not (np.all(x > lower) and np.all(x < upper)):
        raise ValueError("The starting values must lie strictly inside the bounds")
    # An unconstrained maximum strictly inside the bounds is the constrained maximum: bounds that
    # do not bind leave the optimiser's path, and so its result, exactly as without bounds.
    max_iterations = int(maximise_kwargs.pop("max_iterations", 500))
    free = _maximise(fun_vec, x, max_iterations=max_iterations, exact_grad_and_scores=exact_grad_and_scores,
                     **maximise_kwargs)
    if not (lo.any() or hi.any()) or (free.code == 0 and np.all(free.x > lower) and np.all(free.x < upper)):
        return free
    eye = np.eye(x.size)
    A = np.vstack([eye[lo], -eye[hi]])
    b = np.concatenate([-lower[lo], upper[hi]])

    def g(z):
        return A @ z + b

    def fvec(z):
        v = np.asarray(fun_vec(z), dtype=float)
        return v if np.all(np.isfinite(v)) else None

    def f_grad_and_scores(z):
        if exact_grad_and_scores is not None:
            return exact_grad_and_scores(z)
        # Central differences, with each step kept within half the distance to the bounds.
        dist = np.minimum(np.where(lo, z - lower, np.inf), np.where(hi, upper - z, np.inf))
        return _central_scores(fvec, z, np.minimum(1e-5 * np.maximum(1.0, np.abs(z)), 0.5 * dist))

    def penalised(g_old):
        def barrier(z):
            gz = g(z)
            if np.any(gz <= 0):
                return -np.inf
            return mu * float(np.sum(g_old * np.log(gz)) - np.sum(A @ z))

        def r_vec(z):
            return np.append(np.asarray(fun_vec(z), dtype=float), barrier(z))

        def r_grad_and_scores(z):
            gz = g(z)
            if np.any(gz <= 0):
                return None, None
            grad, S = f_grad_and_scores(z)
            if grad is None:
                return None, None
            db = mu * (A.T @ (g_old / gz - 1.0))
            # The rows only seed the curvature model (BHHH); the barrier's own rows make it carry
            # the barrier Hessian A' diag(mu g_old / g^2) A, which dominates next to a bound.
            return grad + db, np.vstack([S, (np.sqrt(mu * g_old) / gz)[:, None] * A])

        return barrier, r_vec, r_grad_and_scores

    used, evals = free.iterations, free.function_evals
    # Start at least a small distance inside each bound. From a point next to a bound (e.g. the
    # previous EM step) the first barrier would be centred within rounding error of the bound,
    # where the inner BFGS cannot converge.
    if used < max_iterations:
        x = _inside(x, lower, upper, 1e-4)
    code = 10
    for _ in range(int(outer_iterations)):
        if used >= max_iterations:
            code = 10
            break
        g_old = g(x)
        barrier, r_vec, r_grad = penalised(g_old)
        r_old = float(np.sum(fun_vec(x))) + barrier(x)
        res = _maximise(r_vec, x, max_iterations=max_iterations - used, exact_grad_and_scores=r_grad,
                        **maximise_kwargs)
        used += res.iterations
        evals += res.function_evals
        x = res.x
        code = res.code
        if code != 0:
            break
        tol = outer_eps * (abs(r_old) + outer_eps)
        # At the inner maximum the multiplier of bound j is mu * (g_old_j / g_j - 1), so closing
        # the remaining gaps g_j gains about mu * sum |g_old - g|. Once that is below the
        # tolerance a further outer iteration cannot improve R by more than the tolerance either,
        # and it would only push x into the rounding error next to an active bound.
        if abs(res.f - r_old) < tol or mu * float(np.sum(np.abs(g_old - g(x)))) < tol:
            break
    # Keep x strictly inside by a margin far below any tolerance, so that rescaling it (the next
    # EM step, the final estimation) cannot round it onto a bound.
    x = _inside(x, lower, upper, 1e-12)
    return _OptResult(x, float(np.sum(fun_vec(x))), code, MESSAGES[code], used, evals)
