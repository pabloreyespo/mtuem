"""Numerical derivatives that follow R's numDeriv package (Richardson extrapolation).

Apollo computes the Hessian with ``numDeriv::hessian`` and the per-individual scores with
``numDeriv::jacobian``. Porting the same algorithms keeps standard errors comparable with the
R version of mtuem.
"""

from collections.abc import Callable

import numpy as np

ZERO_TOL = float(np.sqrt(np.finfo(float).eps / 7e-7))


def _richardson(a: np.ndarray, r: int) -> np.ndarray:
    """Richardson-extrapolate along axis 1 (the r successive step halvings)."""
    for m in range(1, r):
        a = (a[:, 1 : r + 1 - m] * 4**m - a[:, 0 : r - m]) / (4**m - 1)
    return a[:, 0]


def hessian(
    func: Callable[[np.ndarray], float],
    x: np.ndarray,
    eps: float = 1e-4,
    d: float = 0.1,
    zero_tol: float = ZERO_TOL,
    r: int = 4,
    v: float = 2.0,
) -> np.ndarray:
    """Hessian of a scalar function, as ``numDeriv::hessian`` (method "Richardson")."""
    x = np.asarray(x, dtype=float)
    n = x.size
    f0 = float(func(x))
    h0 = np.abs(d * x) + eps * (np.abs(x) < zero_tol)
    hdiag = np.zeros(n)
    H = np.zeros((n, n))

    for i in range(n):
        h = h0.copy()
        approx = np.zeros((1, r))
        e = np.zeros(n)
        for k in range(r):
            e[:] = 0.0
            e[i] = h[i]
            f1 = func(x + e)
            f2 = func(x - e)
            approx[0, k] = (f1 - 2 * f0 + f2) / h[i] ** 2
            h = h / v
        hdiag[i] = _richardson(approx, r)[0]
        H[i, i] = hdiag[i]

    for i in range(n):
        for j in range(i):
            h = h0.copy()
            approx = np.zeros((1, r))
            for k in range(r):
                e = np.zeros(n)
                e[i] = h[i]
                e[j] = h[j]
                f1 = func(x + e)
                f2 = func(x - e)
                approx[0, k] = (f1 - 2 * f0 + f2 - hdiag[i] * h[i] ** 2 - hdiag[j] * h[j] ** 2) / (
                    2 * h[i] * h[j]
                )
                h = h / v
            H[i, j] = H[j, i] = _richardson(approx, r)[0]
    return H


def jacobian(
    func: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    eps: float = 1e-4,
    d: float = 1e-4,
    zero_tol: float = ZERO_TOL,
    r: int = 4,
    v: float = 2.0,
) -> np.ndarray:
    """Jacobian of a vector function, as ``numDeriv::jacobian`` (method "Richardson").

    Returns an array of shape (len(func(x)), len(x)).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    f = np.atleast_1d(np.asarray(func(x), dtype=float))
    a = np.zeros((f.size, r, n))
    h = np.abs(d * x) + eps * (np.abs(x) < zero_tol)
    for k in range(r):
        for i in range(n):
            e = np.zeros(n)
            e[i] = h[i]
            a[:, k, i] = (np.atleast_1d(func(x + e)) - np.atleast_1d(func(x - e))) / (2 * h[i])
        h = h / v
    for m in range(1, r):
        a = (a[:, 1 : r + 1 - m, :] * 4**m - a[:, 0 : r - m, :]) / (4**m - 1)
    return a[:, 0, :]


def gradient(func: Callable[[np.ndarray], float], x: np.ndarray, **kwargs) -> np.ndarray:
    """Gradient of a scalar function, as ``numDeriv::grad`` (method "Richardson")."""
    return jacobian(lambda z: np.atleast_1d(func(z)), x, **kwargs)[0]


def central_gradient(func: Callable[[np.ndarray], float], x: np.ndarray, rel_step: float = 1e-6) -> np.ndarray:
    """Cheap two-point central difference gradient, used inside the optimiser."""
    x = np.asarray(x, dtype=float)
    g = np.zeros_like(x)
    for i in range(x.size):
        h = rel_step * max(1.0, abs(x[i]))
        e = np.zeros_like(x)
        e[i] = h
        g[i] = (func(x + e) - func(x - e)) / (2 * h)
    return g
