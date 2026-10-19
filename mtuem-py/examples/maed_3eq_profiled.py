"""Example 2a: three equations (Tw, Tf1, Ef1), profiled covariance matrix.

Python version of examples/maed_3eq_profiled.r of the R package.
"""

import mtuem

database = mtuem.maed()
database = database[database["EcI"] > 0]

### The covariance is profiled out of the likelihood: only structural parameters
model = mtuem.MTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], goods_cost={"Ef1": 1},
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="profiled",
)

control = mtuem.Control(
    model_name="maed-3eq-profiled",
    description="Three equations with profiled covariance",
    indiv_id="PeID",
    output_dir="output",
)

### Theta (fixed at 1), Phi, thw, th1 (Tf1) and ph1 (Ef1)
beta = model.beta().set(Phi=1, thw=0, th1=0.5, ph1=0.5)

settings = mtuem.EstimationSettings(max_iterations=500)


if __name__ == "__main__":
    ### Single estimation from the starting values above
    result = model.fit(database, beta, control, settings)
    result.summary()

    ### Same model, best of several random starts (start values plus or minus 1)
    multistart = mtuem.MultiStart(n=5, seed=42, normalization="theta_phi")
    result = model.fit(database, beta, control, settings.model_copy(update={"multistart": multistart}))
    result.summary()
    result.save()

    ### Error covariance implied by the residuals
    errors = result.residual_covariance()
    print(errors["covar"])
    print(errors["corr"])

    ### Values of time:
    ###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
    ###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi
    vot = result.values_of_time(silent=False)

    ### Predicted allocations
    print(result.predict().head())
