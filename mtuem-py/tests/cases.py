"""The eight README models, built exactly as tests/reference/make_reference.R builds them in R."""

from pathlib import Path

import mtuem
from mtuem import _engine as engine
from mtuem.config import Control

CASES = ["eq1", "eq3_profiled", "eq3_sig", "eq3_full", "eq3_cholesky", "lc_cholesky", "lc_profiled", "lc_shared"]
CHK = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]


def database():
    db = mtuem.maed()
    return db[db["EcI"] > 0].reset_index(drop=True)


def start_values(case: str) -> tuple[dict, list[str]]:
    if case == "eq1":
        return {"Theta": 1, "Phi": 1, "thw": 0}, ["Theta"]
    base3 = {"Theta": 1, "Phi": 1, "thw": 0, "theta_1": 0.5, "phi_1": 0.5}
    sig = {"sig_1": 100, "sig_2": 100, "sig_3": 100}
    if case == "eq3_profiled":
        return base3, ["Theta"]
    if case == "eq3_sig":
        return base3 | sig, ["Theta"]
    if case == "eq3_full":
        return base3 | sig | {"rho_12": 0, "rho_13": 0, "rho_23": 0}, ["Theta"]
    if case == "eq3_cholesky":
        return base3 | dict(zip(CHK, [100, 0, 100, 0, 0, 100])), ["Theta"]

    base = {"Theta": 1, "Phi": 1, "thw": 0, "th1": 0.3, "ph1": 0.3}
    memb = {"asc": 0, "x_female": 0, "x_older45": 0}
    chk = dict(zip(CHK, [10, 0, 20, 0, 0, 20]))

    def sfx(d, s):
        return {f"{k}_{s}": v for k, v in d.items()}

    if case == "lc_cholesky":
        beta = sfx(base, 1) | sfx(chk, 1) | sfx(memb, 1) | sfx(base, 2) | sfx(chk, 2) | sfx(memb, 2)
    elif case == "lc_profiled":
        beta = sfx(base, 1) | sfx(memb, 1) | sfx(base, 2) | sfx(memb, 2)
    else:
        beta = sfx(base, 1) | sfx(memb, 1) | sfx(base, 2) | sfx(memb, 2) | dict(zip(CHK, [6, 0, 7, 0, 0, 30]))
    return beta, ["Theta_1", "Theta_2", "asc_1", "x_female_1", "x_older45_1"]


def build(case: str, autodiff: bool = False):
    """Return (beta, fixed, probabilities, inputs).

    With ``autodiff=True`` the same model is written with ``mtuem.likelihood_jax`` (log-likelihoods,
    exact derivatives).
    """
    from mtuem import likelihood_jax as lj

    if autodiff:
        loglik, alloc, mix = lj.theta_phi, lj.class_alloc, lj.lc
    else:
        loglik, alloc, mix = mtuem.likelihood.theta_phi, engine.class_alloc, engine.lc
    beta, fixed = start_values(case)
    lc = case.startswith("lc")
    control = Control(model_name=case, description=case, indiv_id="PeID", output_dir=Path("output"))

    if not lc:
        def probabilities(params, inputs, functionality="estimate"):
            b = inputs.attach(params)
            s = {
                "work_times": ["Tw"],
                "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
                "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
            }
            if case != "eq1":
                s |= {"free_times": ["Tf1"], "free_goods": ["Ef1"], "goods_cost": {"Ef1": 1},
                      "times_elasticities": {"theta_1": b.theta_1}, "goods_elasticities": {"phi_1": b.phi_1}}
            if case in ("eq3_sig", "eq3_full"):
                s["sig"] = {"sig_1": b.sig_1, "sig_2": b.sig_2, "sig_3": b.sig_3}
            if case == "eq3_full":
                s["rho"] = {"rho_12": b.rho_12, "rho_13": b.rho_13, "rho_23": b.rho_23}
            if case == "eq3_cholesky":
                s["cholesky"] = {k: getattr(b, k) for k in CHK}
            P = {"model": loglik(s, functionality, inputs)}
            return P

        inputs = engine.make_inputs(database(), beta, fixed, control, autodiff=autodiff)
        return beta, fixed, probabilities, inputs

    def lc_pars(b, inputs):
        lcpars = {p: [getattr(b, f"{p}_1"), getattr(b, f"{p}_2")] for p in ("Theta", "Phi", "thw", "th1", "ph1")}
        if case == "lc_cholesky":
            lcpars |= {p: [getattr(b, f"{p}_1"), getattr(b, f"{p}_2")] for p in CHK}
        V = {
            "class_1": b.asc_1 + b.female * b.x_female_1 + b.older45 * b.x_older45_1,
            "class_2": b.asc_2 + b.female * b.x_female_2 + b.older45 * b.x_older45_2,
        }
        lcpars["pi_values"] = alloc({"classes": {"class_1": 1, "class_2": 2}, "utilities": V})
        return lcpars

    def probabilities(params, inputs, functionality="estimate"):
        b = inputs.attach(params)
        P, class_settings = {}, []
        for s in range(len(b.pi_values)):
            st = {
                "work_times": ["Tw"], "free_times": ["Tf1"], "free_goods": ["Ef1"], "goods_cost": {"Ef1": 1},
                "work_elasticities": {"Theta": b.Theta[s], "Phi": b.Phi[s], "thw": b.thw[s]},
                "times_elasticities": {"th1": b.th1[s]}, "goods_elasticities": {"ph1": b.ph1[s]},
                "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168, "componentName": f"class_{s + 1}",
            }
            if case == "lc_cholesky":
                st["cholesky"] = [getattr(b, k)[s] for k in CHK]
            if case == "lc_shared":
                st["cholesky"] = [getattr(b, k) for k in CHK]
            if case == "lc_profiled":
                class_settings.append(st)
            else:
                P[f"class_{s + 1}"] = loglik(st, functionality, inputs)
        if case == "lc_profiled":
            profiled = lj.lc_profiled if autodiff else mtuem.mtuem_lc_profiled_likelihood
            P["model"] = profiled({"class_settings": class_settings, "classProb": b.pi_values}, functionality, inputs)
        else:
            P["model"] = mix({"inClassProb": P, "classProb": b.pi_values}, inputs, functionality)
        return P

    inputs = engine.make_inputs(database(), beta, fixed, control, lc_pars, autodiff=autodiff)
    return beta, fixed, probabilities, inputs
