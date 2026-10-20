"""Parameter object: starting values, fixed flags and bounds."""

from __future__ import annotations

import difflib
import math
import re
from collections.abc import Iterable, Mapping

from pydantic import model_validator

from mtuem.config import _Strict


class Param(_Strict):
    value: float
    fixed: bool = False
    lower: float | None = None
    upper: float | None = None

    @model_validator(mode="after")
    def _check_bounds(self) -> Param:
        for field in ("value", "lower", "upper"):
            x = getattr(self, field)
            if x is not None and math.isnan(x):
                raise ValueError(f"{field} must be a number, not NaN")
        if self.lower is not None and self.value <= self.lower:
            raise ValueError(f"value {self.value} must be above lower bound {self.lower}")
        if self.upper is not None and self.value >= self.upper:
            raise ValueError(f"value {self.value} must be below upper bound {self.upper}")
        return self


class Beta(_Strict):
    params: dict[str, Param]
    locked: dict[str, str] = {}

    @classmethod
    def from_dict(cls, values: Mapping[str, float], fixed: str | Iterable[str] = ()) -> Beta:
        beta = cls(params={k: Param(value=float(v)) for k, v in values.items()})
        beta.fix([fixed] if isinstance(fixed, str) else list(fixed))
        return beta

    def _resolve(self, name: str) -> list[str]:
        if name in self.params:
            return [name]
        pattern = re.compile(rf"^{re.escape(name)}_\d+$")
        hits = [k for k in self.params if pattern.match(k)]
        if not hits:
            close = difflib.get_close_matches(name, list(self.params), n=5, cutoff=0.5)
            raise KeyError(f"unknown parameter {name!r}; close names: {close}")
        return hits

    def _resolve_all(self, names: str | list[str]) -> list[str]:
        if isinstance(names, str):
            names = [names]
        out: list[str] = []
        for n in names:
            for k in self._resolve(n):
                if k not in out:
                    out.append(k)
        return out

    def _check_unlocked(self, key: str) -> None:
        if key in self.locked:
            raise ValueError(f"parameter {key!r} is locked: {self.locked[key]}")

    def fix(self, names: str | list[str], value: float | None = None) -> Beta:
        keys = self._resolve_all(names)
        if value is not None:
            for k in keys:
                if k in self.locked and value != self.params[k].value:
                    self._check_unlocked(k)
        new = {k: self.params[k].model_copy(update={"fixed": True}) for k in keys}
        if value is not None:
            new = {k: Param(**{**p.model_dump(), "value": float(value)}) for k, p in new.items()}
        self.params.update(new)
        return self

    def free(self, names: str | list[str]) -> Beta:
        keys = self._resolve_all(names)
        for k in keys:
            self._check_unlocked(k)
        for k in keys:
            self.params[k].fixed = False
        return self

    def set(self, **values: float) -> Beta:
        keys = {}
        for n, v in values.items():
            for k in self._resolve(n):
                self._check_unlocked(k)
                keys[k] = float(v)
        new = {k: Param(**{**self.params[k].model_dump(), "value": v}) for k, v in keys.items()}
        self.params.update(new)
        return self

    def lower(self, names: str | list[str], bound: float) -> Beta:
        keys = self._resolve_all(names)
        new = {k: Param(**{**self.params[k].model_dump(), "lower": float(bound)}) for k in keys}
        self.params.update(new)
        return self

    def upper(self, names: str | list[str], bound: float) -> Beta:
        keys = self._resolve_all(names)
        new = {k: Param(**{**self.params[k].model_dump(), "upper": float(bound)}) for k in keys}
        self.params.update(new)
        return self

    @property
    def names(self) -> list[str]:
        return list(self.params)

    @property
    def free_names(self) -> list[str]:
        return [k for k, p in self.params.items() if not p.fixed]

    @property
    def fixed_names(self) -> list[str]:
        return [k for k, p in self.params.items() if p.fixed]

    def values(self) -> dict[str, float]:
        return {k: p.value for k, p in self.params.items()}

    def bounds(self) -> dict[str, tuple[float | None, float | None]]:
        return {k: (p.lower, p.upper) for k, p in self.params.items()}

    def __getitem__(self, name: str) -> Param:
        return self.params[name]

    def __contains__(self, name: object) -> bool:
        return name in self.params
