"""Example 2c: three equations (Tw, Tf1, Ef1), full estimation of sigma and rho.

Python version of examples/maed_3eq_full.r of the R package.
"""

import mtuem

database = mtuem.maed()
database = database[database["EcI"] > 0]

### Standard deviations sig1..sig3 and correlations rho12, rho13, rho23 (upper triangle, by row)
model = mtuem.MTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], goods_cost={"Ef1": 1},
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="full",
)

control = mtuem.Control(
    model_name="maed-3eq-full",
    description="Three equations with full covariance estimation",
    indiv_id="PeID",
    output_dir="output",
)

beta = model.beta().set(
    Phi=1, thw=0, th1=0.5, ph1=0.5,
    sig1=100, sig2=100, sig3=100,
    rho12=0, rho13=0, rho23=0,
)

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
