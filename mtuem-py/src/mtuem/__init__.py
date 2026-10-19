"""mtuem: Microeconomic Time-Use-Expenditure Models in Python.

Describe a model (``MTUEM``, ``LatentClassMTUEM``, or ``CustomModel`` for a hand-written
likelihood), take its parameters with ``model.beta()``, and estimate it with ``model.fit``, which
returns a ``Result``. ``mtuem.load`` restores a saved result. The likelihood functions are in
``mtuem.likelihood`` (NumPy) and ``mtuem.likelihood_jax`` (JAX); post-estimation helpers in
``mtuem.postestimation`` and ``mtuem.lc_profiled``, closed-form allocations in ``mtuem.utils``.
"""

from . import likelihood, likelihood_jax
from ._engine import FallbackWarning
from .beta import Beta, Param
from .config import Control, EstimationSettings, MultiStart
from .data import DATASETS, load_dataset, maed
from .io import load
from .models import MTUEM, CustomModel, LatentClassMTUEM
from .results import Result

__all__ = [
    "Beta",
    "Control",
    "CustomModel",
    "DATASETS",
    "EstimationSettings",
    "FallbackWarning",
    "LatentClassMTUEM",
    "MTUEM",
    "MultiStart",
    "Param",
    "Result",
    "likelihood",
    "likelihood_jax",
    "load",
    "load_dataset",
    "maed",
]
