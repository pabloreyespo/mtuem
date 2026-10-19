"""mtuem: Microeconomic Time-Use-Expenditure Models in Python.

A Python port of the R package mtuem, with a private estimation engine in ``mtuem._engine``.
"""

from . import likelihood, likelihood_jax
from .data import DATASETS, load_dataset, maed
from .io import load
from .lc_profiled import mtuem_get_profiled_sigma, mtuem_lc_profiled_likelihood, mtuem_profile_ic
from .models import MTUEM, LatentClassMTUEM
from .postestimation import cholesky_delta_method, cholesky_to_cov, get_corr, values_of_time
from .results import Result
from .utils import (
    get_additive_quadratic_essentials,
    get_cond_err,
    get_ti_additive_quadratic,
    get_ti_albe,
    get_ti_thph,
    get_tw_additive_quadratic,
    get_tw_albe,
    get_tw_thph,
    get_values_of_time_additive_quadratic,
    get_values_of_time_albe,
    get_values_of_time_thph,
    get_xi_additive_quadratic,
    get_xi_albe,
    get_xi_thph,
)

# Former names of the likelihood functions, kept until the public surface is settled; new code
# uses mtuem.likelihood.theta_phi, .alpha_beta, .additive_quadratic and .residuals.
mtuem_likelihood = likelihood.theta_phi
mtuem_ab_likelihood = likelihood.alpha_beta
additive_quadratic_mtuem_likelihood = likelihood.additive_quadratic
mtuem_residuals = likelihood.residuals
# Former post-estimation names, used by the examples until they move to the declarative API.
mtuem_values_of_time = values_of_time
mtuem_cholesky_deltaMethod = cholesky_delta_method
mtuem_get_corr = get_corr

__all__ = [
    "DATASETS",
    "LatentClassMTUEM",
    "MTUEM",
    "Result",
    "additive_quadratic_mtuem_likelihood",
    "cholesky_delta_method",
    "cholesky_to_cov",
    "get_additive_quadratic_essentials",
    "get_cond_err",
    "get_corr",
    "get_ti_additive_quadratic",
    "get_ti_albe",
    "get_ti_thph",
    "get_tw_additive_quadratic",
    "get_tw_albe",
    "get_tw_thph",
    "get_values_of_time_additive_quadratic",
    "get_values_of_time_albe",
    "get_values_of_time_thph",
    "get_xi_additive_quadratic",
    "get_xi_albe",
    "get_xi_thph",
    "likelihood",
    "likelihood_jax",
    "load",
    "load_dataset",
    "maed",
    "mtuem_ab_likelihood",
    "mtuem_cholesky_deltaMethod",
    "mtuem_get_corr",
    "mtuem_get_profiled_sigma",
    "mtuem_lc_profiled_likelihood",
    "mtuem_likelihood",
    "mtuem_profile_ic",
    "mtuem_residuals",
    "mtuem_values_of_time",
    "values_of_time",
]
