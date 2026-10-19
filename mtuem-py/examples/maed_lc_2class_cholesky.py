"""Example 3a: latent class model, 2 classes x 3 equations, per-class Cholesky covariance.

Estimated with an EM warm start, then from several random starts, each warm started with a
few EM iterations. Python version of examples/maed_lc_2class_cholesky.r of the R package.
"""

import numpy as np
import pandas as pd

import mtuem
from mtuem.postestimation import cholesky_to_cov

database = mtuem.maed()
database = database[database["EcI"] > 0]

### Two classes, each with its own Theta-Phi parameters and Cholesky factor (suffix _1, _2).
### Class membership is a logit on female and older45, with class 1 as reference.
model = mtuem.LatentClassMTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"],
    Tc="Tc", Ec="EcI", w="w", tau=168,
    covariance="cholesky",
    n_classes=2,
    allocation=["female", "older45"],
)

control = mtuem.Control(
    model_name="maed-lc-2class-3eq",
    description="Latent class model with 2 classes and 3 equations each",
    indiv_id="PeID",
    output_dir="output",
)

### model.beta() fixes Theta_1, Theta_2 and the class 1 allocation (asc_1, x_female_1,
### x_older45_1). A base name such as chk22 sets the entry of every class.
beta = model.beta().set(
    Phi=1, thw=0, th1=0.3, ph1=0.3,
    chk11=10,
    chk21=0, chk22=20,
    chk31=0, chk32=0, chk33=20,
    asc=0, x_female=0, x_older45=0,
)

CHK = ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]
eq_names = ["Tw", "Tf1", "Ef1"]


if __name__ == "__main__":
    ### Single estimation: 10 EM iterations, then BFGS
    em = mtuem.EstimationSettings(method="em", em_max_iterations=10, max_iterations=500)
    result = model.fit(database, beta, control, em)
    result.summary()

    ### Best of several random starts, each warm started with 10 EM iterations. fit raises a
    ### RuntimeError when no candidate reaches a non-saddle optimum: then raise n, or set
    ### non_saddle=False to inspect the saddle solutions (in the error's .attempts).
    multistart = mtuem.MultiStart(n=5, seed=42, first_em=True, em_iter_max=10, normalization="theta_phi")
    settings = mtuem.EstimationSettings(max_iterations=500, multistart=multistart)
    result = model.fit(database, beta, control, settings)
    result.summary()
    result.save()

    # ---------------------------------------------------------------------------------------
    # Decoding the Cholesky factor, class by class
    # ---------------------------------------------------------------------------------------
    for s in (1, 2):
        names = [f"{k}_{s}" for k in CHK]
        Sigma = pd.DataFrame(cholesky_to_cov(result.estimates, names), index=eq_names, columns=eq_names)
        sd = np.sqrt(np.diag(Sigma))
        print(f"\n--- Class {s} ---")
        print(Sigma)
        print(pd.Series(sd, index=eq_names))
        print(Sigma / np.outer(sd, sd))

    ### Same quantities with standard errors, one block per class
    print(result.correlations())

    ### Class shares and posterior class membership
    print(result.class_shares())
    print(result.posterior().head())

    ### Values of time by class, averaged with the posterior class membership weights
    vot = result.values_of_time(silent=False)
