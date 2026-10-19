"""Example 2d: three equations (Tw, Tf1, Ef1), Cholesky parameterisation (recommended).

Python version of examples/maed_3eq_cholesky.r of the R package.
"""

import numpy as np
import pandas as pd

import mtuem
from mtuem.postestimation import cholesky_to_cov

database = mtuem.maed()
database = database[database["EcI"] > 0]

### The covariance is L L', with L lower triangular: chk11, chk21, chk22, chk31, chk32, chk33
model = mtuem.MTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], goods_cost={"Ef1": 1},
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="cholesky",
)

control = mtuem.Control(
    model_name="maed-3eq-cholesky",
    description="Three equations with Cholesky covariance estimation",
    indiv_id="PeID",
    output_dir="output",
)

beta = model.beta().set(
    Phi=1, thw=0, th1=0.5, ph1=0.5,
    chk11=100,
    chk21=0, chk22=100,
    chk31=0, chk32=0, chk33=100,
)

settings = mtuem.EstimationSettings(max_iterations=500)

cholesky_names = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]
eq_names = ["Tw", "Tf1", "Ef1"]


if __name__ == "__main__":
    ### Single estimation from the starting values above
    result = model.fit(database, beta, control, settings)
    result.summary()

    ### Same model, best of several random starts (start values plus or minus 1)
    multistart = mtuem.MultiStart(n=5, seed=42, normalization="theta_phi")
    result = model.fit(database, beta, control, settings.model_copy(update={"multistart": multistart}))
    result.summary()
    result.save()

    # -------------------------------------------------------------------------------------
    # Decoding the Cholesky factor
    # -------------------------------------------------------------------------------------
    # chk11...chk33 are the non-zero entries of the lower triangular L, read row by row.
    # Sigma = L L', the standard deviations are sqrt(diag(Sigma)) and the correlations
    # follow from dividing by them.
    Sigma = pd.DataFrame(cholesky_to_cov(result.estimates, cholesky_names), index=eq_names, columns=eq_names)
    sd = np.sqrt(np.diag(Sigma))
    print(Sigma)
    print(pd.Series(sd, index=eq_names))
    print(Sigma / np.outer(sd, sd))

    ### Same quantities with standard errors, via the delta method
    print(result.correlations())

    ### Values of time:
    ###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
    ###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi
    vot = result.values_of_time(silent=False)

    ### Predicted allocations
    print(result.predict().head())
