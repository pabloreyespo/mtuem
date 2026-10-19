"""Example 2d: three equations (Tw, Tf1, Ef1), Cholesky parameterisation (recommended).

Python port of examples/maed_3eq_cholesky.r.
"""

import numpy as np
import pandas as pd

import mtuem
from mtuem import apollo

apollo_control = {
    "modelName": "maed-3eq-cholesky",
    "modelDescr": "Three equations with Cholesky covariance estimation",
    "indivID": "PeID",
    "outputDirectory": "output",
}

est_set = {
    "writeIter": False,
    "maxIterations": 500,
    "estimationRoutine": "bgw",
    "hessianRoutine": "numDeriv",
}

database = mtuem.maed()
database = database[database["EcI"] > 0]

### Cholesky parameters are given row-wise over the lower triangle
apollo_beta = {
    "Theta": 1, "Phi": 1, "thw": 0,
    "theta_1": 0.5, "phi_1": 0.5,
    "chk11": 100,
    "chk21": 0, "chk22": 100,
    "chk31": 0, "chk32": 0, "chk33": 100,
}
apollo_fixed = ["Theta"]

cholesky_names = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]
eq_names = ["Tw", "Tf1", "Ef1"]

apollo_inputs = apollo.validate_inputs(database, apollo_beta, apollo_fixed, apollo_control)


def apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"):
    b = apollo_inputs.attach(apollo_beta)

    P = {}

    mtuem_settings = {
        "work_times": ["Tw"],
        "free_times": ["Tf1"],
        "free_goods": ["Ef1"],
        "goods_cost": {"Ef1": 1},
        "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
        "times_elasticities": {"theta_1": b.theta_1},
        "goods_elasticities": {"phi_1": b.phi_1},
        "cholesky": {k: getattr(b, k) for k in cholesky_names},
        "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
    }

    P["model"] = mtuem.mtuem_likelihood(mtuem_settings, functionality, apollo_inputs)
    return apollo.prepare_prob(P, apollo_inputs, functionality)


if __name__ == "__main__":
    ### Single estimation from the starting values above
    model = apollo.estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, est_set)
    apollo.model_output(model)

    ### Starting values module: same model, best of several random starts
    multi_start = mtuem.customMultiStart(
        apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
        customMultistart_settings={
            "apolloBetaMax": {k: v + 1 for k, v in apollo_beta.items()},
            "apolloBetaMin": {k: v - 1 for k, v in apollo_beta.items()},
            "nCandidates": 5,
        },
        estimate_settings=est_set,
        non_saddle=True,
        normalization={
            "normalization": "theta_phi",
            "work_elasticities": ["Theta", "Phi", "thw"],
            "times_elasticities": ["theta_1"],
            "goods_elasticities": ["phi_1"],
            "cholesky": cholesky_names,
        },
        seed=42,
    )
    model = multi_start["best_model"]
    apollo.model_output(model)
    apollo.save_output(model)

    # -------------------------------------------------------------------------------------
    # Decoding the Cholesky factor
    # -------------------------------------------------------------------------------------
    # chk11...chk33 are the non-zero entries of the lower triangular L, read row by row.
    # Sigma = L L', the standard deviations are sqrt(diag(Sigma)) and the correlations
    # follow from dividing by them.
    Sigma = pd.DataFrame(mtuem.cholesky_to_cov(model.estimate, cholesky_names), index=eq_names, columns=eq_names)
    sd = np.sqrt(np.diag(Sigma))
    print(Sigma)
    print(pd.Series(sd, index=eq_names))
    print(Sigma / np.outer(sd, sd))

    ### Same quantities with standard errors, via the delta method
    print(mtuem.mtuem_cholesky_deltaMethod(model, cholesky_names, eq_names))

    ### Values of time:
    ###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
    ###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi
    vot = mtuem.mtuem_values_of_time(model, apollo_inputs, Ec="EcI")

    ### Predicted allocations
    pred = apollo.prediction(model, apollo_probabilities, apollo_inputs)
    print(pred.head())
