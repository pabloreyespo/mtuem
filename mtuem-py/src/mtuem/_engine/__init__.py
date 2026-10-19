"""Private estimation engine: inputs, optimiser, latent class tools, output and post-estimation.

Not part of the public API. The public likelihood and estimation functions build on it.
"""

from .em import lc_em
from .estimate import Estimate, estimate
from .inputs import Inputs, make_inputs
from .latent import ClassProbs, class_alloc, lc, lc_conditionals
from .output import output_text, save_files
from .postestimation import delta_method, prediction


class FallbackWarning(UserWarning):
    """A model asked for the JAX backend but a likelihood form has no JAX version: NumPy is used."""


__all__ = [
    "ClassProbs",
    "Estimate",
    "FallbackWarning",
    "Inputs",
    "class_alloc",
    "delta_method",
    "estimate",
    "lc",
    "lc_conditionals",
    "lc_em",
    "make_inputs",
    "output_text",
    "prediction",
    "save_files",
]
