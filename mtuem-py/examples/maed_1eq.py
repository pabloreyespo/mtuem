"""Example 1: single equation (work only). Theta is fixed to 1 for identification.

Python version of examples/maed_1eq.r of the R package.
"""

import mtuem

database = mtuem.maed()
database = database[database["EcI"] > 0]

### The model: one work equation, Theta-Phi utility
model = mtuem.MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", tau=168)

control = mtuem.Control(
    model_name="maed-1eq",
    description="Single equation MTUEM (work only)",
    indiv_id="PeID",
    output_dir="output",
)

### Parameters Theta (fixed at 1 by model.beta()), Phi and thw
beta = model.beta().set(Phi=1, thw=0)

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

    ### Values of time:
    ###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
    ###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi
    vot = result.values_of_time(silent=False)

    ### Predicted allocations
    print(result.predict().head())
