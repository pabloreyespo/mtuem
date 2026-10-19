"""Example 3b: latent class model, 2 classes x 3 equations, covariance profiled out.

mtuem_lc_profiled_likelihood concentrates the class-specific covariance out of the likelihood
with an inner fixed point, so apollo_beta holds only structural and membership parameters. It
already returns the class mixture, so apollo.lc is bypassed and EM cannot be used.
Python port of examples/maed_lc_2class_profiled.r.
"""

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
    "modelName": "maed-lc-2class-profiled-cov",
    "modelDescr": "Latent class model with 2 classes and a profiled-out error covariance",
    "indivID": "PeID",
    "outputDirectory": "output",
    "noValidation": True,
    "noDiagnostics": True,
}

database = mtuem.maed()
database = database[database["EcI"] > 0]

apollo_beta = {}
for s in (1, 2):
    apollo_beta |= {
        # Class s structural
        f"Theta_{s}": 1, f"Phi_{s}": 1, f"thw_{s}": 0, f"th1_{s}": 0.3, f"ph1_{s}": 0.3,
        # Class s membership
        f"asc_{s}": 0, f"x_female_{s}": 0, f"x_older45_{s}": 0,
    }

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
    class_settings = []
    for s in range(len(b.pi_values)):
        class_settings.append({
            "work_times": ["Tw"],
            "free_times": ["Tf1"],
            "free_goods": ["Ef1"],
            "goods_cost": {"Ef1": 1},
            "work_elasticities": {"Theta": b.Theta[s], "Phi": b.Phi[s], "thw": b.thw[s]},
            "times_elasticities": {"th1": b.th1[s]},
            "goods_elasticities": {"ph1": b.ph1[s]},
            "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
            "componentName": f"class_{s + 1}",
        })

    lc_settings = {"class_settings": class_settings, "classProb": b.pi_values}
    P["model"] = mtuem.mtuem_lc_profiled_likelihood(lc_settings, functionality, apollo_inputs)
    return apollo.prepare_prob(P, apollo_inputs, functionality)


if __name__ == "__main__":
    ### Single estimation from the starting values above
    model = apollo.estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, est_set)
    apollo.model_output(model)

    ### Starting values module (first_em=False: apollo.lc is bypassed, so there is no EM structure)
    multi_starts = mtuem.customMultiStart(
        apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
        customMultistart_settings={
            "apolloBetaMax": {k: v + 1 for k, v in apollo_beta.items()},
            "apolloBetaMin": {k: v - 1 for k, v in apollo_beta.items()},
            "nCandidates": 5,
        },
        estimate_settings=est_set,
        first_em=False,
        verbose=False,
        non_saddle=True,
        normalization={
            "normalization": "theta_phi",
            "work_elasticities": ["Theta", "Phi", "thw"],
            "times_elasticities": ["th1"],
            "goods_elasticities": ["ph1"],
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

    ### Re-evaluate at the converged estimates so the profiled covariance matches them
    apollo_probabilities(model.estimate.to_dict(), apollo_inputs, "estimate")

    profiled = mtuem.mtuem_get_profiled_sigma()
    for s in (0, 1):
        print(f"\n--- Profiled-out error covariance (Class {s + 1}) ---")
        print("sig:\n", profiled["sig"][s])
        print("rho:\n", profiled["rho"][s])
    print("\nClass shares:\n", profiled["shares"])

    ### AIC/BIC must count the profiled-out covariance parameters: nClass * nEq * (nEq + 1) / 2
    print(mtuem.mtuem_profile_ic(model, nClass=2, nEq=3))

    ### Values of time by class, averaged with the posterior class membership weights
    vot = mtuem.mtuem_values_of_time(model, apollo_inputs, apollo_probabilities, nClass=2, Ec="EcI")
