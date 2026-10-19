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


def _sided_column(fn: Callable[[np.ndarray], np.ndarray | None], x: np.ndarray, i: int, h: float,
                  room_lo: float, room_hi: float, floor: float, f0: np.ndarray | None = None):
    """Derivative of the vector ``fn`` in coordinate ``i`` by a difference of step ``h``.

    Central with step ``h`` when both bounds are at least ``10 h`` away. Nearer a bound the step
    is a tenth of the distance to it, so that f singular at the bound (a log, say) still gets
    its local derivative; below ``floor``, the step at which roundoff would show, it is the
    second-order one-sided difference away from the bound with step ``floor``. Both are exact
    for a quadratic, so an optimum closer to the bound than the step keeps the right sign. None
    if a value is not finite.
    """
    e = np.zeros(x.size)
    near, far = min(room_lo, room_hi), max(room_lo, room_hi)
    if 0.1 * near >= max(h, floor):
        e[i] = min(h, 0.1 * near)
        up, down = fn(x + e), fn(x - e)
        return None if up is None or down is None else (up - down) / (2 * e[i])
    e[i] = (1.0 if room_hi >= room_lo else -1.0) * min(floor, 0.4 * far)
    f0 = fn(x) if f0 is None else f0
    f1, f2 = fn(x + e), fn(x + 2 * e)
    if f0 is None or f1 is None or f2 is None:
        return None
    return (-3 * f0 + 4 * f1 - f2) / (2 * e[i])


def _roundoff_floor(f0: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Difference steps ``1e-10 max(1, sum|f0|) max(1, |x|)``, at which roundoff in a sum of the size
    of ``f0`` stays below 1e-5 of a difference."""
    return 1e-10 * max(1.0, float(np.sum(np.abs(f0)))) * np.maximum(1.0, np.abs(x))


def _sided_scores(fn: Callable[[np.ndarray], np.ndarray | None], x: np.ndarray, steps: np.ndarray,
                  room_lo: np.ndarray, room_hi: np.ndarray, floor: np.ndarray | None = None):
    """Like ``_central_scores`` but with the steps of ``_sided_column`` next to a bound.

    ``floor`` defaults to ``_roundoff_floor`` for the values of ``fn`` at ``x``.
    """
    S = None
    f0 = fn(x)
    if f0 is None:
        return None, None
    if floor is None:
        floor = _roundoff_floor(f0, x)
    for i in range(x.size):
        col = _sided_column(fn, x, i, steps[i], room_lo[i], room_hi[i], floor[i], f0)
        if col is None:
            return None, None
        if S is None:
            S = np.zeros((col.size, x.size))
        S[:, i] = col
    return S.sum(axis=0), S


def _curvature_root(H: np.ndarray) -> np.ndarray:
    """``S`` with ``S'S`` the Hessian of -f made positive definite (eigenvalues made positive)."""
    w, V = np.linalg.eigh(-0.5 * (H + H.T))
    w = np.maximum(np.abs(w), 1e-6 * max(float(np.max(np.abs(w))), 1e-300))
    return np.sqrt(w)[:, None] * V.T


def _at_maximum(f: float, g: np.ndarray, H: np.ndarray) -> bool:
    """The test ``_maximise`` applies when no step improves f, with a Newton model: the Hessian
    is negative definite and the gain it predicts is below 1e-7 * max(|f|, 1)."""
    try:
        np.linalg.cholesky(-H)
    except np.linalg.LinAlgError:
        return False
    return 0.5 * float(g @ np.linalg.solve(-H, g)) < 1e-7 * max(abs(f), 1.0)


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
    active_tol: float = 1e-6,
    **maximise_kwargs,
) -> _OptResult:
    """Maximise ``sum(fun_vec(x))`` subject to ``lower < x < upper`` (infinite entries are no bound).

    1. The plain optimiser from ``x0``. An unconstrained maximum strictly inside the bounds is
       the constrained maximum and is returned as is, so bounds that do not bind leave the
       result exactly as without bounds.
    2. Otherwise the adaptive logarithmic barrier of R's ``constrOptim``, which maxLik uses for
       Apollo's inequality constraints. With ``g(x)`` the stacked distances ``x - lower`` and
       ``upper - x`` for the finite bounds and ``A = dg/dx``, each outer iteration maximises

           R(x) = f(x) + mu * sum(g(x_k) * log g(x) - A x)

       from ``x_k``. It stops when an outer iteration improves R by less than
       ``outer_eps * (|R(x_k)| + outer_eps)`` (constrOptim's rule) or when every bound the
       iteration closed in on (gap at least halved) is within ``x_tol * max(1, |bound|)``; the
       second rule keeps f near zero from pushing gaps into the rounding error of x. The barrier
       enters ``_maximise`` as one extra individual with an exact gradient, its Hessian
       ``A' diag(mu g(x_k) / g^2) A`` is added to the BHHH seed, and the inner runs accept a
       stall (see ``_maximise``).
    3. An active-set finish, which removes the barrier's ``mu`` bias and the stiffness it has
       next to a bound. The parameters the barrier left within ``active_tol * max(1, |bound|)``
       of a bound are fixed exactly at it, and the plain optimiser maximises f over the others.
       A free parameter that this takes out of its bounds joins the active set; a fixed one
       whose gradient points into the feasible region (relative size above the convergence
       threshold of ``_maximise``) is released; then the finish is repeated. It ends when the
       maximum over the free parameters is inside the bounds and every fixed parameter's
       gradient points out of them (the Karush-Kuhn-Tucker conditions).

    Code 0 means that the last run of step 3 converged and the KKT conditions hold; a run that
    does not converge ends the search with its code (8, 9 or 10), and an active set that does
    not settle within ``2 n + 2`` rounds gives code 8. A false convergence (8) of a step 3 run
    is accepted only if a Newton model from a finite-difference Hessian predicts no gain, and
    is otherwise retried once with that Hessian as the curvature seed. ``max_iterations`` is
    shared by all runs, and a budget used up gives code 10. Without ``exact_grad_and_scores``
    the derivatives of step 3 are the differences of ``_sided_column``: central, with steps
    shrunk next to a bound and one-sided (second order) at it, so that an optimum closer to a
    bound than the usual step is still seen; the barrier runs use central differences with
    steps kept inside the bounds. Where f is not finite on a bound itself, the parameters
    fixed there sit 1e-12 inside it. The returned ``f`` is the unpenalised objective, and
    ``x`` is kept at least 1e-12 (relative) inside the bounds, so that rescaling it cannot
    round it onto a bound.
    """
    x = np.asarray(x0, dtype=float).copy()
    lower = np.broadcast_to(np.asarray(lower, dtype=float), x.shape)
    upper = np.broadcast_to(np.asarray(upper, dtype=float), x.shape)
    lo, hi = np.isfinite(lower), np.isfinite(upper)
    if not (np.all(x > lower) and np.all(x < upper)):
        raise ValueError("The starting values must lie strictly inside the bounds")
    maximise_kwargs.pop("accept_stall", None)
    max_iterations = int(maximise_kwargs.pop("max_iterations", 500))
    free = _maximise(fun_vec, x, max_iterations=max_iterations, exact_grad_and_scores=exact_grad_and_scores,
                     **maximise_kwargs)
    if not (lo.any() or hi.any()) or (free.code == 0 and np.all(free.x > lower) and np.all(free.x < upper)):
        return free
    used, evals = free.iterations, free.function_evals

    def result(z, code):
        z = _inside(z, lower, upper, 1e-12)
        return _OptResult(z, float(np.sum(fun_vec(z))), code, MESSAGES[code], used, evals)

    if used >= max_iterations:
        # No budget left: the free run's point, brought inside the bounds (or the start).
        xf = _inside(free.x, lower, upper, 1e-12)
        return result(xf if np.all(np.isfinite(fun_vec(xf))) else x, 10)

    n = x.size
    eye = np.eye(n)
    A = np.vstack([eye[lo], -eye[hi]])
    b = np.concatenate([-lower[lo], upper[hi]])
    bound_scale = np.maximum(1.0, np.abs(np.where(lo, lower, np.where(hi, upper, 0.0))))

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
        # Central differences, with each step kept within half the distance to the bounds.
        return _central_scores(fvec, z, np.minimum(1e-5 * np.maximum(1.0, np.abs(z)), 0.5 * gap(z)))

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

    # Step 2: the barrier.
    for _ in range(int(outer_iterations)):
        if used >= max_iterations:
            return result(x, 10)
        r_vec, r_grad = penalised(g(x))
        r_old = float(np.sum(r_vec(x)))
        res = _maximise(r_vec, x, max_iterations=max_iterations - used, exact_grad_and_scores=r_grad,
                        accept_stall=True, **maximise_kwargs)
        used += res.iterations
        evals += res.function_evals
        if res.code in (9, 10):
            return result(res.x, res.code)
        closing = (gap(res.x) < 0.5 * gap(x)) & (lo | hi)
        x = res.x
        if res.code != 0 or abs(res.f - r_old) < outer_eps * (abs(r_old) + outer_eps):
            break
        if closing.any() and np.all(gap(x)[closing] <= x_tol * bound_scale[closing]):
            break

    # Step 3: active-set finish.
    at_lo = lo & (x - lower <= active_tol * bound_scale)
    at_hi = hi & ~at_lo & (upper - x <= active_tol * bound_scale)

    def bound_gradient(z, fixed):
        """Derivatives of f at the fixed parameters, taken from inside the bounds; NaN elsewhere."""
        d = np.full(n, np.nan)
        if exact_grad_and_scores is not None:
            grad, _ = exact_grad_and_scores(z)
            if grad is not None:
                d[fixed] = grad[fixed]
            return d
        f0 = fvec(z)
        if f0 is None:
            return d
        floor = _roundoff_floor(f0, z)
        for i in np.flatnonzero(fixed):
            col = _sided_column(fvec, z, i, 1e-5 * max(1.0, abs(z[i])), z[i] - lower[i], upper[i] - z[i],
                                floor[i], f0)
            if col is not None:
                d[i] = col.sum()
        return d

    def polish(z, idx, budget):
        """Maximise f over the parameters ``idx`` with the others held at ``z``.

        f is taken less its value at the start, spread over the individuals: the tolerances of
        _maximise, relative to max(|f|, 1), then bound the gain still available in absolute
        terms, as a constrained maximum next to a bound needs whatever the level of f. A false
        convergence (8) is checked with a Newton model from a finite-difference Hessian: either
        the point passes the convergence tests of _maximise with that model, or _maximise is run
        again from it with the finite-difference Hessian as its curvature seed. With few
        individuals BHHH is degenerate at a maximum, and far from one it can be so poor that no
        step along its direction improves f.
        """
        v0 = fvec(z)
        shift = float(v0.sum()) / v0.size if v0 is not None else 0.0

        def full(y):
            w = z.copy()
            w[idx] = y
            return w

        def sub_vec(y):
            return np.asarray(fun_vec(full(y)), dtype=float) - shift

        def sub_fvec(y):
            v = sub_vec(y)
            return v if np.all(np.isfinite(v)) else None

        def sub_grad(y):
            if exact_grad_and_scores is not None:
                grad, S = exact_grad_and_scores(full(y))
                return (None, None) if grad is None else (grad[idx], S[:, idx])
            return _sided_scores(sub_fvec, y, 1e-5 * np.maximum(1.0, np.abs(y)), y - lower[idx], upper[idx] - y)

        def grad_and_hessian(y):
            g0, _ = sub_grad(y)
            if g0 is None:
                return None, None
            steps = 1e-4 * np.maximum(1.0, np.abs(y))
            _, H = _sided_scores(lambda w: sub_grad(w)[0], y, steps, y - lower[idx], upper[idx] - y, steps)
            return g0, None if H is None else 0.5 * (H + H.T)

        res = _maximise(sub_vec, z[idx], max_iterations=budget, exact_grad_and_scores=sub_grad, **maximise_kwargs)
        if res.code != 8:
            return res
        g0, H = grad_and_hessian(res.x)
        if H is None:
            return res
        if _at_maximum(res.f, g0, H):
            res.code, res.message = 0, MESSAGES[0]
            return res
        if res.iterations >= budget:
            res.code, res.message = 10, MESSAGES[10]
            return res
        root = _curvature_root(H)

        def seeded(y):
            grad, _ = sub_grad(y)
            return (None, None) if grad is None else (grad, root)

        again = _maximise(sub_vec, res.x, max_iterations=budget - res.iterations, exact_grad_and_scores=seeded,
                          **maximise_kwargs)
        again.iterations += res.iterations
        again.function_evals += res.function_evals
        if again.code == 8:
            g0, H = grad_and_hessian(again.x)
            if H is not None and _at_maximum(again.f, g0, H):
                again.code, again.message = 0, MESSAGES[0]
        return again

    for _ in range(2 * n + 2):
        if used >= max_iterations:
            return result(x, 10)
        fixed = at_lo | at_hi
        z = np.where(at_lo, lower, np.where(at_hi, upper, x))
        if fvec(z) is None:
            # f is not defined on the bound itself (a log, say): the fixed parameters sit 1e-12
            # inside it, where the returned x is clipped to anyway.
            z = np.where(fixed, _inside(z, lower, upper, 1e-12), z)
            if fvec(z) is None:
                return result(x, 8)
        if not fixed.all():
            idx = np.flatnonzero(~fixed)
            res = polish(z, idx, max_iterations - used)
            used += res.iterations
            evals += res.function_evals
            z[idx] = res.x
            if res.code != 0:
                return result(z, res.code)
        out_lo, out_hi = lo & ~fixed & (z <= lower), hi & ~fixed & (z >= upper)
        if out_lo.any() or out_hi.any():
            # Free parameters the maximum takes out of their bounds join the active set.
            at_lo, at_hi = at_lo | out_lo, at_hi | out_hi
            x = np.clip(z, lower, upper)
            continue
        x = z
        grad = bound_gradient(z, fixed)
        if not np.all(np.isfinite(grad[fixed])):
            return result(z, 8)
        # The gradient test of _maximise on the shifted f (|f| <= 1): an absolute threshold.
        size = np.abs(grad) * np.maximum(np.abs(z), 1.0)
        inward = fixed & (size >= 1e-4) & np.where(at_lo, grad > 0, grad < 0)
        if not inward.any():
            return result(z, 0)
        # Release the fixed parameters whose gradient points into the feasible region.
        at_lo, at_hi = at_lo & ~inward, at_hi & ~inward
        x = _inside(z, lower, upper, 1e-8)
    return result(x, 8)
