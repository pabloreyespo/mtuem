"""Example 2c: three equations (Tw, Tf1, Ef1), full estimation of sigma and rho.

Python port of examples/maed_3eq_full.r.
"""

import mtuem
from mtuem import apollo

apollo_control = {
    "modelName": "maed-3eq-full",
    "modelDescr": "Three equations with full covariance estimation",
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

### Structural + sigma + rho parameters
### rho ordering for 3 equations: rho_12, rho_13, rho_23 (horizontal upper triangle)
apollo_beta = {
    "Theta": 1, "Phi": 1, "thw": 0,
    "theta_1": 0.5, "phi_1": 0.5,
    "sig_1": 100, "sig_2": 100, "sig_3": 100,
    "rho_12": 0, "rho_13": 0, "rho_23": 0,
}
apollo_fixed = ["Theta"]

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
        "sig": {"sig_1": b.sig_1, "sig_2": b.sig_2, "sig_3": b.sig_3},
        "rho": {"rho_12": b.rho_12, "rho_13": b.rho_13, "rho_23": b.rho_23},
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
            "sig": ["sig_1", "sig_2", "sig_3"],
            "rho": ["rho_12", "rho_13", "rho_23"],
        },
        seed=42,
    )
    model = multi_start["best_model"]
    apollo.model_output(model)
    apollo.save_output(model)

    ### Error covariance implied by the residuals
    corr_result = mtuem.mtuem_get_corr(model, apollo_probabilities, apollo_inputs)
    print(corr_result["covar"])
    print(corr_result["corr"])

    ### Values of time:
    ###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
    ###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi
    vot = mtuem.mtuem_values_of_time(model, apollo_inputs, Ec="EcI")

    ### Predicted allocations
    pred = apollo.prediction(model, apollo_probabilities, apollo_inputs)
    print(pred.head())
