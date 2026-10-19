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
    accept_stall: bool = False,
) -> _OptResult:
    """Maximise a sum of per-individual log-likelihoods, ``sum(fun_vec(x))``.

    Quasi-Newton (BFGS) with a backtracking line search. The curvature model is seeded with the
    BHHH matrix (outer product of the individual scores), which the central differences used for
    the gradient deliver at no extra cost. Non-finite values are treated as an over-long step, so
    the search backs away from regions where the likelihood is zero. Code 0 is convergence (the
    gain predicted by the quadratic model is below ``rel_f_tol`` times |f|); 8, 9 and 10 report
    false convergence and the function evaluation and iteration limits.

    ``accept_stall`` (used for the barrier subproblems of ``maximise_bounded``): when the
    quasi-Newton model and, right after it, a fresh BHHH model both predict a negligible gain
    while the gradient test still fails, the run counts as converged. Next to an active bound
    the barrier's curvature makes the gradient test unattainable there although no further gain
    is available.
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
    stalls = 0  # consecutive iterations with a negligible predicted gain but a large gradient
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
            stalls += 1
            if accept_stall and stalls >= 2:
                code = 0
                break
            # The quasi-Newton matrix has lost curvature along directions where the gradient
            # is still large: re-seed it with BHHH.
            Hinv = bhhh_inverse(S)
            restarted = True
        else:
            stalls = 0
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
    x_tol: float = 1e-8,
    **maximise_kwargs,
) -> _OptResult:
    """Maximise ``sum(fun_vec(x))`` subject to ``lower < x < upper`` (infinite entries are no bound).

    The adaptive logarithmic barrier of R's ``constrOptim``, which maxLik uses for Apollo's
    inequality constraints. With ``g(x)`` the stacked distances ``x - lower`` and ``upper - x``
    for the finite bounds and ``A = dg/dx``, each outer iteration maximises

        R(x) = f(x) + mu * sum(g(x_k) * log g(x) - A x)

    from ``x_k`` with ``_maximise``. The barrier gradient ``mu * A'(g(x_k) / g(x) - 1)`` vanishes
    at ``x_k``, so the penalty fades as the iterates settle and the limit is the constrained
    maximum, also when it lies on a bound. The barrier is used only when the unconstrained
    maximum (from ``x0``) is not strictly inside the bounds; otherwise that maximum is returned
    as is.

    Stopping: an outer iteration ends the loop when it improves R by less than
    ``outer_eps * (|R(x_k)| + outer_eps)`` (constrOptim's rule), or when every active bound (a
    parameter with one bound whose gap the iteration at least halved) is within
    ``x_tol * max(1, |bound|)``. The second
    rule is scale-free in f; without it, f near zero runs further outer iterations that only
    push the gap into the rounding error of x, where the inner BFGS cannot converge.

    Active bounds and the barrier's stiffness: next to an active bound the barrier curvature
    ``mu g(x_k) / g^2`` dwarfs that of f, which BFGS learns slowly and the relative tolerances
    of ``_maximise`` do not resolve. Since f is near-linear over such a small gap, the
    equilibrium of R along an active parameter alone is ``g = mu g_k / (mu - f')`` (lower bound;
    ``mu + f'`` for an upper bound). Each later outer iteration moves the active parameters
    there before its inner run, and each inner run is followed by the same step for its own R;
    both steps are kept only if they raise R. The barrier enters ``_maximise`` as one extra
    individual with an exact gradient; the BHHH seed gets the barrier Hessian
    ``A' diag(mu g(x_k) / g^2) A`` as extra rows, and the inner runs accept a stall (see
    ``_maximise``), which that curvature otherwise turns into endless BHHH re-seeds. The scores of ``fun_vec`` are central differences (with steps kept inside the
    bounds where the full step gives a non-finite value) unless ``exact_grad_and_scores`` is
    given.

    ``max_iterations`` is shared by the unconstrained run and the inner barrier runs. The code
    is that of the last inner run: an inner run that does not converge (false convergence, or
    the iteration or evaluation limit) ends the loop with its code, and a budget used up before
    the barrier can run gives code 10. The returned ``f`` is the unpenalised objective, and
    ``x`` is kept at least 1e-12 (relative) inside the bounds.
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
    used, evals = free.iterations, free.function_evals
    if used >= max_iterations:
        # No budget left for the barrier: the free run's point, brought inside the bounds.
        xf = _inside(free.x, lower, upper, 1e-12)
        v = np.asarray(fun_vec(xf), dtype=float)
        if not np.all(np.isfinite(v)):
            xf, v = x, np.asarray(fun_vec(x), dtype=float)
        return _OptResult(xf, float(np.sum(v)), 10, MESSAGES[10], used, evals)

    eye = np.eye(x.size)
    A = np.vstack([eye[lo], -eye[hi]])
    b = np.concatenate([-lower[lo], upper[hi]])
    only_lo, only_hi = lo & ~hi, hi & ~lo
    tol_x = x_tol * np.maximum(1.0, np.abs(np.where(lo, lower, np.where(hi, upper, 0.0))))

    def g(z):
        return A @ z + b

    def gap(z):
        """Distance of each parameter to its nearer finite bound (inf if unbounded)."""
        return np.minimum(np.where(lo, z - lower, np.inf), np.where(hi, upper - z, np.inf))

    def fvec(z):
        v = np.asarray(fun_vec(z), dtype=float)
        return v if np.all(np.isfinite(v)) else None

    def f_grad_and_scores(z):
        if exact_grad_and_scores is not None:
            return exact_grad_and_scores(z)
        steps = 1e-5 * np.maximum(1.0, np.abs(z))
        grad, S = _central_scores(fvec, z, steps)
        if grad is None:
            # The function ends at a bound: keep each step within half the distance to it.
            grad, S = _central_scores(fvec, z, np.minimum(steps, 0.5 * gap(z)))
        return grad, S

    def penalised(g_old):
        def barrier(z):
            gz = g(z)
            if np.any(gz <= 0):
                return -np.inf
            return mu * float(np.sum(g_old * np.log(gz)) - np.sum(A @ z))

        def r_vec(z):
            return np.append(np.asarray(fun_vec(z), dtype=float), barrier(z))

        def r_grad_and_scores(z):
            if np.any(g(z) <= 0):
                return None, None
            grad, S = f_grad_and_scores(z)
            if grad is None:
                return None, None
            db = mu * (A.T @ (g_old / g(z) - 1.0))
            return grad + db, np.vstack([S, (np.sqrt(mu * g_old) / g(z))[:, None] * A])

        return r_vec, r_grad_and_scores

    def settle(z, r, xc, active, r_vec):
        """Active parameters moved to the equilibrium of the barrier centred at ``xc``, if R rises."""
        grad, _ = f_grad_and_scores(z)
        if grad is None:
            return z, r
        glo_c, ghi_c = xc - lower, upper - xc
        with np.errstate(divide="ignore", invalid="ignore"):
            g_lo = np.where(mu - grad > 0, mu * glo_c / (mu - grad), np.inf)
            g_hi = np.where(mu + grad > 0, mu * ghi_c / (mu + grad), np.inf)
            zn = np.where(active & only_lo & (g_lo < z - lower), lower + g_lo, z)
            zn = np.where(active & only_hi & (g_hi < upper - z), upper - g_hi, zn)
        if np.array_equal(zn, z) or not (np.all(zn > lower) and np.all(zn < upper)):
            return z, r
        v = np.asarray(r_vec(zn), dtype=float)
        rn = float(np.sum(v))
        return (zn, rn) if np.all(np.isfinite(v)) and rn >= r else (z, r)

    code = 10
    active = np.zeros(x.size, dtype=bool)
    for _ in range(int(outer_iterations)):
        if used >= max_iterations:
            code = 10
            break
        r_vec, r_grad = penalised(g(x))
        r_old = float(np.sum(r_vec(x)))
        start, r_start = settle(x, r_old, x, active, r_vec) if active.any() else (x, r_old)
        res = _maximise(r_vec, start, max_iterations=max_iterations - used, exact_grad_and_scores=r_grad,
                        accept_stall=True, **maximise_kwargs)
        used += res.iterations
        evals += res.function_evals
        code = res.code
        if code != 0:
            x = res.x
            break
        active = (gap(res.x) < 0.5 * gap(x)) & (only_lo | only_hi)
        xn, r_new = settle(res.x, res.f, x, active, r_vec)
        x = xn
        if abs(r_new - r_old) < outer_eps * (abs(r_old) + outer_eps):
            break
        if active.any() and np.all(gap(x)[active] <= tol_x[active]):
            break
    # Keep x strictly inside by a margin far below any tolerance, so that rescaling it (the next
    # EM step, the final estimation) cannot round it onto a bound.
    x = _inside(x, lower, upper, 1e-12)
    return _OptResult(x, float(np.sum(fun_vec(x))), code, MESSAGES[code], used, evals)
