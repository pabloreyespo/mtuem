"""Example 3c: latent class model, 2 classes x 3 equations, one Cholesky factor shared by both.

Parameters listed in shared= have no class suffix and are common to every class.
Python version of examples/maed_lc_2class_shared.r of the R package.
"""

import numpy as np
import pandas as pd

import mtuem
from mtuem.postestimation import cholesky_to_cov

database = mtuem.maed()
database = database[database["EcI"] > 0]

cholesky_names = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]
eq_names = ["Tw", "Tf1", "Ef1"]

### Two classes with their own Theta-Phi parameters and a SINGLE shared covariance structure
model = mtuem.LatentClassMTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], goods_cost={"Ef1": 1},
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="cholesky",
    n_classes=2,
    allocation=["female", "older45"],
    shared=cholesky_names,
)

control = mtuem.Control(
    model_name="maed-lc-2class-shared-cov",
    description="Latent class model with 2 classes sharing a single error covariance",
    indiv_id="PeID",
    output_dir="output",
)

beta = model.beta().set(
    Phi=1, thw=0, th1=0.3, ph1=0.3, asc=0, x_female=0, x_older45=0,
    chk11=6, chk21=0, chk22=7, chk31=0, chk32=0, chk33=30,
)

settings = mtuem.EstimationSettings(max_iterations=500)


if __name__ == "__main__":
    ### Single estimation from the starting values above
    result = model.fit(database, beta, control, settings)
    result.summary()

    ### Best of several random starts
    multistart = mtuem.MultiStart(n=5, seed=42, normalization="theta_phi")
    result = model.fit(database, beta, control, settings.model_copy(update={"multistart": multistart}))
    result.summary()
    result.save()

    # Decoding the shared Cholesky factor
    Sigma = pd.DataFrame(cholesky_to_cov(result.estimates, cholesky_names), index=eq_names, columns=eq_names)
    sd = np.sqrt(np.diag(Sigma))
    print(Sigma)
    print(pd.Series(sd, index=eq_names))
    print(Sigma / np.outer(sd, sd))
    print(result.correlations())

    ### Values of time by class, averaged with the posterior class membership weights
    vot = result.values_of_time(silent=False)
