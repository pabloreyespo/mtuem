"""Example 2b: three equations (Tw, Tf1, Ef1), explicit standard deviations, no correlations.

Python version of examples/maed_3eq_sig.r of the R package.
"""

import mtuem

database = mtuem.maed()
database = database[database["EcI"] > 0]

### One standard deviation per equation: sig1 (Tw), sig2 (Tf1), sig3 (Ef1)
model = mtuem.MTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], goods_cost={"Ef1": 1},
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="diagonal",
)

control = mtuem.Control(
    model_name="maed-3eq-sig",
    description="Three equations with explicit sigma",
    indiv_id="PeID",
    output_dir="output",
)

beta = model.beta().set(Phi=1, thw=0, th1=0.5, ph1=0.5, sig1=100, sig2=100, sig3=100)

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
