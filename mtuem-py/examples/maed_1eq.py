"""Example 1: single equation (work only). Theta is fixed to 1 for identification.

Python port of examples/maed_1eq.r.
"""

import mtuem
from mtuem import apollo

apollo_control = {
    "modelName": "maed-1eq",
    "modelDescr": "Single equation MTUEM (work only)",
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

apollo_beta = {"Theta": 1, "Phi": 1, "thw": 0}
apollo_fixed = ["Theta"]

apollo_inputs = apollo.validate_inputs(database, apollo_beta, apollo_fixed, apollo_control)


def apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"):
    b = apollo_inputs.attach(apollo_beta)

    P = {}

    work_elasticities = {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw}

    mtuem_settings = {
        "work_times": ["Tw"],
        "work_elasticities": work_elasticities,
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
        },
        seed=42,
    )
    model = multi_start["best_model"]
    apollo.model_output(model)
    apollo.save_output(model)

    ### Values of time:
    ###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
    ###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi
    vot = mtuem.mtuem_values_of_time(model, apollo_inputs, Ec="EcI")

    ### Predicted allocations
    pred = apollo.prediction(model, apollo_probabilities, apollo_inputs)
    print(pred.head())
