"""Example 3b: latent class model, 2 classes x 3 equations, covariance profiled out.

With covariance="profiled" the class-specific covariances are concentrated out of the
likelihood by an inner fixed point, so the parameters are only the structural and membership
ones. The likelihood mixes the classes itself, so EM cannot be used.
Python version of examples/maed_lc_2class_profiled.r of the R package.
"""

import mtuem

database = mtuem.maed()
database = database[database["EcI"] > 0]

model = mtuem.LatentClassMTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], goods_cost={"Ef1": 1},
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="profiled",
    n_classes=2,
    allocation=["female", "older45"],
)

control = mtuem.Control(
    model_name="maed-lc-2class-profiled-cov",
    description="Latent class model with 2 classes and a profiled-out error covariance",
    indiv_id="PeID",
    output_dir="output",
)

beta = model.beta().set(Phi=1, thw=0, th1=0.3, ph1=0.3, asc=0, x_female=0, x_older45=0)

settings = mtuem.EstimationSettings(max_iterations=500)


if __name__ == "__main__":
    ### Single estimation from the starting values above
    result = model.fit(database, beta, control, settings)
    result.summary()

    ### Best of several random starts (no EM warm start: the profiled likelihood has no class
    ### components)
    multistart = mtuem.MultiStart(n=5, seed=42, normalization="theta_phi")
    result = model.fit(
        database, beta, control, settings.model_copy(update={"multistart": multistart, "verbose": False})
    )
    result.summary()
    result.save()

    ### The profiled-out covariances at the estimates
    profiled = result.profiled_covariance()
    for s in (0, 1):
        print(f"\n--- Profiled-out error covariance (Class {s + 1}) ---")
        print("sig:\n", profiled["sig"][s])
        print("rho:\n", profiled["rho"][s])
    print("\nClass shares:\n", profiled["shares"])

    ### AIC/BIC counting the profiled-out covariance parameters: n_classes * M * (M + 1) / 2
    print({k: profiled[k] for k in ("k", "AIC", "BIC")})

    ### Values of time by class, averaged with the posterior class membership weights
    vot = result.values_of_time(silent=False)
