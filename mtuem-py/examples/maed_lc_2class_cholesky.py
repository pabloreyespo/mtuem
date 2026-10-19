"""Example 3a: latent class model, 2 classes x 3 equations, per-class Cholesky covariance.

Estimated with an EM warm start (apollo.lc_em), then with customMultiStart(first_em=True).
Python port of examples/maed_lc_2class_cholesky.r.
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
    "hessianRoutine": "numDeriv",
}

apollo_control = {
    "modelName": "maed-lc-2class-3eq",
    "modelDescr": "Latent class model with 2 classes and 3 equations each",
    "indivID": "PeID",
    "outputDirectory": "output",
    "noValidation": True,
    "noDiagnostics": True,
}

database = mtuem.maed()
database = database[database["EcI"] > 0]

CHK = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]

apollo_beta = {}
for s in (1, 2):
    apollo_beta |= {
        # Class s structural
        f"Theta_{s}": 1, f"Phi_{s}": 1, f"thw_{s}": 0, f"th1_{s}": 0.3, f"ph1_{s}": 0.3,
        # Class s error structure
        f"chk11_{s}": 10,
        f"chk21_{s}": 0, f"chk22_{s}": 20,
        f"chk31_{s}": 0, f"chk32_{s}": 0, f"chk33_{s}": 20,
        # Class s membership
        f"asc_{s}": 0, f"x_female_{s}": 0, f"x_older45_{s}": 0,
    }

apollo_fixed = ["Theta_1", "Theta_2", "asc_1", "x_female_1", "x_older45_1"]


### Class membership and class-specific parameters
def apollo_lcPars(b, apollo_inputs):
    lcpars = {p: [getattr(b, f"{p}_1"), getattr(b, f"{p}_2")] for p in ["Theta", "Phi", "thw", "th1", "ph1", *CHK]}

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
    for s in range(len(b.pi_values)):
        mtuem_settings = {
            "work_times": ["Tw"],
            "free_times": ["Tf1"],
            "free_goods": ["Ef1"],
            "work_elasticities": {"Theta": b.Theta[s], "Phi": b.Phi[s], "thw": b.thw[s]},
            "times_elasticities": {"th1": b.th1[s]},
            "goods_elasticities": {"ph1": b.ph1[s]},
            "cholesky": {k: getattr(b, k)[s] for k in CHK},
            "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
            "componentName": f"class_{s + 1}",
        }
        P[f"class_{s + 1}"] = mtuem.mtuem_likelihood(mtuem_settings, functionality, apollo_inputs)

    lc_settings = {"inClassProb": P, "classProb": b.pi_values}
    P["model"] = apollo.lc(lc_settings, apollo_inputs, functionality)
    return apollo.prepare_prob(P, apollo_inputs, functionality)


if __name__ == "__main__":
    eq_names = ["Tw", "Tf1", "Ef1"]

    ### Single estimation with an EM warm start
    model = apollo.lc_em(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                         lcEM_settings={"EMmaxIterations": 10}, estimate_settings=est_set)
    apollo.model_output(model)

    ### Starting values module: best of several random starts, each warm started with lc_em
    multi_starts = mtuem.customMultiStart(
        apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
        customMultistart_settings={
            "apolloBetaMax": {k: v + 1 for k, v in apollo_beta.items()},
            "apolloBetaMin": {k: v - 1 for k, v in apollo_beta.items()},
            "nCandidates": 5,
        },
        estimate_settings=est_set,
        first_em=True,
        em_iter_max=10,
        non_saddle=True,
        normalization={
            "normalization": "theta_phi",
            "work_elasticities": ["Theta", "Phi", "thw"],
            "times_elasticities": ["th1"],
            "goods_elasticities": ["ph1"],
            "cholesky": CHK,
            "covs": ["asc", "x_female", "x_older45"],
        },
        nClass=2,
        seed=42,
    )

    ### best_model is None when no candidate reaches a non-saddle optimum
    if multi_starts["best_model"] is None:
        raise RuntimeError(
            "No candidate converged to a non-saddle optimum. Increase nCandidates, widen "
            "apolloBetaMin/apolloBetaMax, or set non_saddle=False to inspect the saddle solutions."
        )
    model = multi_starts["best_model"]
    apollo.model_output(model)
    apollo.save_output(model)

    # ---------------------------------------------------------------------------------------
    # Decoding the Cholesky factor, class by class
    # ---------------------------------------------------------------------------------------
    for s in (1, 2):
        names = [f"{k}_{s}" for k in CHK]
        Sigma = pd.DataFrame(mtuem.cholesky_to_cov(model.estimate, names), index=eq_names, columns=eq_names)
        sd = np.sqrt(np.diag(Sigma))
        print(f"\n--- Class {s} ---")
        print(Sigma)
        print(pd.Series(sd, index=eq_names))
        print(Sigma / np.outer(sd, sd))
        ### Same quantities with standard errors
        print(mtuem.mtuem_cholesky_deltaMethod(model, names, eq_names))

    ### Values of time by class, averaged with the posterior class membership weights
    vot = mtuem.mtuem_values_of_time(model, apollo_inputs, apollo_probabilities, nClass=2, Ec="EcI")
