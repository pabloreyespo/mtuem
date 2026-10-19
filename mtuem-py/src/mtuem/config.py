"""Pydantic configuration objects for mtuem estimation."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Control(_Strict):
    model_name: str = "model"
    description: str = ""
    indiv_id: str | None = None
    output_dir: Path = Path(".")


class MultiStart(_Strict):
    n: int = Field(default=100, gt=0)
    seed: int | None = None
    lower: dict[str, float] | None = None
    upper: dict[str, float] | None = None
    first_em: bool = False
    em_iter_max: int = Field(default=5, gt=0)
    non_saddle: bool = True
    normalization: Literal["theta_phi", "alpha_beta"] | None = None


class EstimationSettings(_Strict):
    method: Literal["bfgs", "em"] = "bfgs"
    fallback: Literal["bfgs"] | None = "bfgs"
    hessian: Literal["analytic", "numeric", "none"] | None = None
    max_iterations: int = Field(default=200, gt=0)
    max_function_evals: int = Field(default=100000, gt=0)
    em_max_iterations: int = Field(default=100, gt=0)
    em_tolerance: float = Field(default=1e-5, gt=0)
    saddle_escapes: int = Field(default=3, ge=0)
    backend: Literal["jax", "numpy"] = "jax"
    multistart: MultiStart | None = None
    verbose: bool = True

    @model_validator(mode="after")
    def _check_hessian_backend(self) -> EstimationSettings:
        if self.backend == "numpy" and self.hessian == "analytic":
            raise ValueError("hessian='analytic' requires backend='jax'")
        return self

    def resolved_hessian(self) -> str:
        if self.hessian is not None:
            return self.hessian
        return "analytic" if self.backend == "jax" else "numeric"
