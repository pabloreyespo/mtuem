import pytest

from mtuem.beta import Beta, Param


def make() -> Beta:
    return Beta.from_dict({"Theta_1": 0.5, "Theta_2": 0.6, "chk21_1": 1.0, "chk2_1": 2.0, "gamma": 0.3})


def test_fix_base_name_covers_all_classes():
    b = make().fix("Theta")
    assert b.fixed_names == ["Theta_1", "Theta_2"]
    assert "Theta_1" not in b.free_names


def test_fix_exact_name():
    b = make().fix("gamma")
    assert b.fixed_names == ["gamma"]


def test_fix_with_value():
    b = make().fix("gamma", 0.9)
    assert b["gamma"].value == 0.9 and b["gamma"].fixed


def test_base_name_does_not_prefix_match():
    b = Beta.from_dict({"chk21_1": 1.0})
    with pytest.raises(KeyError):
        b.fix("chk2")
    assert b.fixed_names == []


def test_unknown_name_suggests():
    b = Beta.from_dict({"Phi_1": 1.0, "Theta_1": 0.5})
    with pytest.raises(KeyError, match="Phi_1"):
        b.fix("Phi_9")


def test_lower_bound_replaces_constraint():
    b = make().lower("Theta", 0.0).lower("Theta", 0.1)
    assert b.bounds()["Theta_1"] == (0.1, None)
    assert b.bounds()["gamma"] == (None, None)
    b.upper("Theta", 2.0)
    assert b.bounds()["Theta_2"] == (0.1, 2.0)


def test_start_outside_bounds_raises():
    b = make()
    with pytest.raises(ValueError):
        b.lower("Theta", 0.55)
    with pytest.raises(ValueError):
        b.upper("Theta", 0.5)
    assert b.bounds()["Theta_1"] == (None, None)
    b.lower("gamma", 0.0)
    with pytest.raises(ValueError):
        b.set(gamma=-1.0)


def test_set_and_values():
    b = make().set(gamma=0.4)
    assert b.values()["gamma"] == 0.4
    assert "gamma" in b and "zzz" not in b


def test_from_dict_marks_fixed():
    b = Beta.from_dict({"a_1": 1.0, "b_1": 2.0}, fixed=["b"])
    assert b.fixed_names == ["b_1"] and b.free_names == ["a_1"]


def test_json_round_trip():
    b = make().fix("Theta").lower("gamma", 0.0)
    b.locked["gamma"] = "normalisation"
    b2 = Beta.model_validate_json(b.model_dump_json())
    assert b2 == b
    assert b2.locked == {"gamma": "normalisation"}


def test_locked_parameter_cannot_be_freed():
    b = make()
    b.fix("gamma")
    b.locked["gamma"] = "identification"
    with pytest.raises(ValueError, match="identification"):
        b.free("gamma")
    with pytest.raises(ValueError, match="identification"):
        b.set(gamma=1.0)
    b.fix("gamma", 0.3)
    with pytest.raises(ValueError, match="identification"):
        b.fix("gamma", 0.8)
    assert b.fixed_names == ["gamma"]


def test_from_dict_fixed_accepts_a_single_name():
    b = Beta.from_dict({"Theta_1": 0.5, "Theta_2": 0.6, "gamma": 0.3}, fixed="Theta")
    assert b.fixed_names == ["Theta_1", "Theta_2"]


@pytest.mark.parametrize("field", ["value", "lower", "upper"])
def test_param_rejects_nan(field):
    kwargs = {"value": 0.0} | {field: float("nan")}
    with pytest.raises(ValueError, match="NaN"):
        Param(**kwargs)


def test_beta_rejects_nan_start_and_bounds():
    with pytest.raises(ValueError, match="NaN"):
        Beta.from_dict({"a": float("nan")})
    b = Beta.from_dict({"a": 0.0})
    with pytest.raises(ValueError, match="NaN"):
        b.lower("a", float("nan"))
    with pytest.raises(ValueError, match="NaN"):
        b.set(a=float("nan"))
    assert b.bounds()["a"] == (None, None) and b.values()["a"] == 0.0


def test_param_allows_infinite_bounds():
    p = Param(value=0.0, lower=float("-inf"), upper=float("inf"))
    assert p.lower == float("-inf") and p.upper == float("inf")
