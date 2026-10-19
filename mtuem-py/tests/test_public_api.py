"""The public surface of mtuem: exported names, naming rules and the runnable examples."""

import importlib
import inspect
import os
import pkgutil
import subprocess
import sys
from pathlib import Path

import pytest

import mtuem

SRC = Path(mtuem.__file__).parent
EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

PUBLIC = {
    "Control", "EstimationSettings", "MultiStart", "Beta", "Param", "MTUEM", "LatentClassMTUEM",
    "CustomModel", "Result", "load", "likelihood", "likelihood_jax", "FallbackWarning", "DATASETS",
    "load_dataset", "maed",
}

# Names of the former Apollo-style surface that must not come back.
OLD = [
    "mtuem_likelihood", "mtuem_ab_likelihood", "additive_quadratic_mtuem_likelihood", "mtuem_residuals",
    "mtuem_values_of_time", "mtuem_get_corr", "mtuem_cholesky_deltaMethod", "mtuem_lc_profiled_likelihood",
    "mtuem_get_profiled_sigma", "mtuem_profile_ic", "customMultiStart", "apollo",
]


def _modules():
    """mtuem and every submodule, private ones included (their names are checked too)."""
    mods = [mtuem]
    for info in pkgutil.walk_packages(mtuem.__path__, prefix="mtuem."):
        mods.append(importlib.import_module(info.name))
    return mods


def _public_modules():
    return [m for m in _modules() if not any(p.startswith("_") for p in m.__name__.split("."))]


def _names_of(obj, label: str):
    """(label, name) for the attributes of a module or class, and the parameters of functions."""
    out = []
    members = vars(obj).items() if inspect.ismodule(obj) else inspect.getmembers(obj)
    for name, value in members:
        if name.startswith("__") and name.endswith("__"):
            continue
        out.append((f"{label}.{name}", name))
        if inspect.isclass(value) and value.__module__.startswith("mtuem"):
            for attr, member in inspect.getmembers(value):
                if attr.startswith("__") and attr.endswith("__"):
                    continue
                out.append((f"{label}.{name}.{attr}", attr))
                if callable(member):
                    out += _parameters(member, f"{label}.{name}.{attr}")
            out += [(f"{label}.{name}.{f}", f) for f in getattr(value, "model_fields", {})]
        elif inspect.isfunction(value):
            out += _parameters(value, f"{label}.{name}")
    return out


def _parameters(fn, label: str):
    try:
        params = inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return []
    return [(f"{label}({p})", p) for p in params]


def test_all_is_the_public_surface():
    assert set(mtuem.__all__) == PUBLIC
    assert len(mtuem.__all__) == len(set(mtuem.__all__))
    for name in mtuem.__all__:
        assert getattr(mtuem, name) is not None


def test_the_exports_are_the_classes_of_their_modules():
    from mtuem import beta, config, io, models, results
    from mtuem._engine import FallbackWarning

    assert mtuem.Control is config.Control and mtuem.EstimationSettings is config.EstimationSettings
    assert mtuem.MultiStart is config.MultiStart
    assert mtuem.Beta is beta.Beta and mtuem.Param is beta.Param
    assert mtuem.MTUEM is models.MTUEM and mtuem.LatentClassMTUEM is models.LatentClassMTUEM
    assert mtuem.CustomModel is models.CustomModel and mtuem.Result is results.Result
    assert mtuem.load is io.load and mtuem.FallbackWarning is FallbackWarning
    assert issubclass(mtuem.FallbackWarning, UserWarning)


@pytest.mark.parametrize("name", OLD)
def test_the_former_names_are_gone(name):
    assert not hasattr(mtuem, name)


# Helpers that belong to one backend only.
NUMPY_ONLY = {"residuals", "works_in_logs"}
JAX_ONLY = {"require_jax_form"}
LIKELIHOODS = {"theta_phi", "alpha_beta", "additive_quadratic", "class_alloc", "lc", "lc_profiled"}


def _own_functions(module):
    """Public functions defined in ``module`` itself (not imported helpers)."""
    return {n for n, v in vars(module).items()
            if not n.startswith("_") and inspect.isfunction(v) and v.__module__ == module.__name__}


def test_numpy_and_jax_likelihood_modules_expose_the_same_functions():
    numpy_names, jax_names = set(mtuem.likelihood.__all__), set(mtuem.likelihood_jax.__all__)
    assert numpy_names - NUMPY_ONLY == jax_names - JAX_ONLY == LIKELIHOODS
    for module in (mtuem.likelihood, mtuem.likelihood_jax):
        assert all(inspect.isfunction(getattr(module, n)) for n in module.__all__)
        assert _own_functions(module) <= set(module.__all__)


def test_numpy_latent_class_functions_are_the_engine_ones():
    from mtuem import _engine
    from mtuem import lc_profiled

    assert mtuem.likelihood.class_alloc is _engine.class_alloc
    assert mtuem.likelihood.lc is _engine.lc
    for module in (mtuem.likelihood, mtuem.likelihood_jax):
        assert list(inspect.signature(module.lc_profiled).parameters) == ["lc_settings", "functionality", "inputs"]
        assert list(inspect.signature(module.lc).parameters) == ["lc_settings", "inputs", "functionality"]
        assert list(inspect.signature(module.class_alloc).parameters) == ["classAlloc_settings"]
    assert inspect.signature(mtuem.likelihood.lc_profiled) == inspect.signature(lc_profiled.lc_profiled)


def test_load_signature():
    params = list(inspect.signature(mtuem.load).parameters)
    assert params == ["path", "data", "model_class", "model"]
    defaults = {k: p.default for k, p in inspect.signature(mtuem.load).parameters.items()}
    assert defaults["data"] is None and defaults["model_class"] is None and defaults["model"] is None


def test_no_public_function_carries_the_mtuem_prefix():
    hits = [label for m in _public_modules() for label, name in _names_of(m, m.__name__)
            if name.lower().startswith("mtuem_") and not name.startswith("_")]
    assert hits == []


def test_profiled_latent_class_functions_have_plain_names():
    from mtuem import lc_profiled

    for name in ("lc_profiled", "get_profiled_sigma", "profile_ic"):
        assert callable(getattr(lc_profiled, name))


def test_no_public_name_mentions_apollo():
    hits = [label for m in _modules() for label, name in _names_of(m, m.__name__) if "apollo" in name.lower()]
    hits += [m.__name__ for m in _modules() if "apollo" in m.__name__.lower()]
    hits += [str(p.relative_to(SRC)) for p in SRC.rglob("*")
             if "__pycache__" not in p.parts and "apollo" in p.name.lower()]
    assert hits == []


def test_the_name_walk_sees_classes_methods_and_parameters():
    # Guards the test above against walking nothing.
    labels = {label for m in _modules() for label, _ in _names_of(m, m.__name__)}
    assert "mtuem.models.MTUEM.fit" in labels
    assert "mtuem.models.MTUEM.fit(settings)" in labels
    assert "mtuem.models.LatentClassMTUEM.n_classes" in labels
    assert "mtuem.io.load(model_class)" in labels
    assert "mtuem.likelihood.theta_phi(functionality)" in labels


def test_the_engine_names_no_apollo_in_its_messages():
    text = "\n".join(p.read_text(encoding="utf-8") for p in SRC.rglob("*.py"))
    # Apollo may be cited in docstrings and comments, but not as an identifier of the old API.
    for old in ("apollo_inputs", "apollo_probabilities", "apollo_beta", "apollo_fixed", "apolloBeta",
                "apollo_control", "apollo_lcPars", "--extra jax"):
        assert old not in text


EXAMPLE_FILES = sorted(p.name for p in EXAMPLES.glob("maed_*.py"))


def test_every_example_is_listed():
    assert len(EXAMPLE_FILES) == 8


# Each example fits its model twice, once from random starts, in a new interpreter that compiles
# the JAX likelihood again: about 30 s for a single-class example, minutes for a latent class one.
@pytest.mark.slow
@pytest.mark.parametrize("name", EXAMPLE_FILES)
def test_examples_run(name, tmp_path):
    # Run from a temporary folder, so the files the example saves (under output/) land there.
    env = os.environ | {"PYTHONIOENCODING": "utf-8"}
    run = subprocess.run(
        [sys.executable, str(EXAMPLES / name)], cwd=tmp_path, env=env, capture_output=True, text=True,
        encoding="utf-8", timeout=3600,
    )
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    assert "apollo" not in run.stdout.lower()
    assert list((tmp_path / "output").glob("*_model.json"))
