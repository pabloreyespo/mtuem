"""Example 3c: latent class model, 2 classes x 3 equations, one Cholesky factor shared by both.

The shared parameters have no class suffix and are not listed in apollo_lcPars.
Python port of examples/maed_lc_2class_shared.r.
"""

import numpy as np
import pandas as pd

import mtuem
from mtuem import apollo

est_set = {
    "writeIter": False,
    "silent": False,
    "maxIterations": 500,
    "scaleHessian": False,
    "scaleAfterConvergence": False,
    "estimationRoutine": "bgw",
    "hessianRoutine": "numDeriv",
}

apollo_control = {
    "modelName": "maed-lc-2class-shared-cov",
    "modelDescr": "Latent class model with 2 classes sharing a single error covariance",
    "indivID": "PeID",
    "outputDirectory": "output",
    "noValidation": True,
    "noDiagnostics": True,
}

database = mtuem.maed()
database = database[database["EcI"] > 0]

cholesky_names = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]
eq_names = ["Tw", "Tf1", "Ef1"]

### Parameters for 2 classes, with a SINGLE shared covariance structure
apollo_beta = {}
for s in (1, 2):
    apollo_beta |= {
        f"Theta_{s}": 1, f"Phi_{s}": 1, f"thw_{s}": 0, f"th1_{s}": 0.3, f"ph1_{s}": 0.3,
        f"asc_{s}": 0, f"x_female_{s}": 0, f"x_older45_{s}": 0,
    }
apollo_beta |= {"chk11": 6, "chk21": 0, "chk22": 7, "chk31": 0, "chk32": 0, "chk33": 30}

apollo_fixed = ["Theta_1", "Theta_2", "asc_1", "x_female_1", "x_older45_1"]


def apollo_lcPars(b, apollo_inputs):
    lcpars = {p: [getattr(b, f"{p}_1"), getattr(b, f"{p}_2")] for p in ["Theta", "Phi", "thw", "th1", "ph1"]}

    V = {
        "class_1": b.asc_1 + b.female * b.x_female_1 + b.older45 * b.x_older45_1,
        "class_2": b.asc_2 + b.female * b.x_female_2 + b.older45 * b.x_older45_2,
    }
    lcpars["pi_values"] = apollo.class_alloc({"classes": {"class_1": 1, "class_2": 2}, "utilities": V})
    return lcpars


apollo_inputs = apollo.validate_inputs(database, apollo_beta, apollo_fixed, apollo_control, apollo_lcPars)


def apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"):
    b = apollo_inputs.attach(apollo_beta)

    P = {}
    chk = {k: getattr(b, k) for k in cholesky_names}
    for s in range(len(b.pi_values)):
        mtuem_settings = {
            "work_times": ["Tw"],
            "free_times": ["Tf1"],
            "free_goods": ["Ef1"],
            "goods_cost": {"Ef1": 1},
            "work_elasticities": {"Theta": b.Theta[s], "Phi": b.Phi[s], "thw": b.thw[s]},
            "times_elasticities": {"th1": b.th1[s]},
            "goods_elasticities": {"ph1": b.ph1[s]},
            "cholesky": chk,
            "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
            "componentName": f"class_{s + 1}",
        }
        P[f"class_{s + 1}"] = mtuem.mtuem_likelihood(mtuem_settings, functionality, apollo_inputs)

    lc_settings = {"inClassProb": P, "classProb": b.pi_values}
    P["model"] = apollo.lc(lc_settings, apollo_inputs, functionality)
    return apollo.prepare_prob(P, apollo_inputs, functionality)


if __name__ == "__main__":
    ### Single estimation from the starting values above
    model = apollo.estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, est_set)
    apollo.model_output(model)

    ### Starting values module: best of several random starts
    multi_starts = mtuem.customMultiStart(
        apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
        customMultistart_settings={
            "apolloBetaMax": {k: v + 1 for k, v in apollo_beta.items()},
            "apolloBetaMin": {k: v - 1 for k, v in apollo_beta.items()},
            "nCandidates": 5,
        },
        estimate_settings=est_set,
        first_em=False,
        non_saddle=True,
        normalization={
            "normalization": "theta_phi",
            "work_elasticities": ["Theta", "Phi", "thw"],
            "times_elasticities": ["th1"],
            "goods_elasticities": ["ph1"],
            "cholesky": cholesky_names,
            "covs": ["asc", "x_female", "x_older45"],
        },
        nClass=2,
        seed=42,
    )

    if multi_starts["best_model"] is None:
        raise RuntimeError(
            "No candidate converged to a non-saddle optimum. Increase nCandidates, widen "
            "apolloBetaMin/apolloBetaMax, or set non_saddle=False to inspect the saddle solutions."
        )
    model = multi_starts["best_model"]
    apollo.model_output(model)
    apollo.save_output(model)

    # Decoding the shared Cholesky factor
    Sigma = pd.DataFrame(mtuem.cholesky_to_cov(model.estimate, cholesky_names), index=eq_names, columns=eq_names)
    sd = np.sqrt(np.diag(Sigma))
    print(Sigma)
    print(pd.Series(sd, index=eq_names))
    print(Sigma / np.outer(sd, sd))
    print(mtuem.mtuem_cholesky_deltaMethod(model, cholesky_names, eq_names))

    ### Values of time by class, averaged with the posterior class membership weights
    vot = mtuem.mtuem_values_of_time(model, apollo_inputs, apollo_probabilities, nClass=2, Ec="EcI")
