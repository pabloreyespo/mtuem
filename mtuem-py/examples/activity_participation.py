"""Activity participation: a latent class MTUEM whose classes differ in the activities they do.

Python version of the model of the article on latent classes of activity participation
(R/models/hobbies_cholesky.R and R/models/hobbies_rest.R of the article's repository). Each class
takes part in some free time activities and not in others: in the 2015 model, class 1 has no
hobbies (free time Tf2) and class 2 does everything. Individuals who spend hours on an activity
are unlikely to belong to a class that does not do it, and a penalty in the class likelihood
says so.

The file shows how to customize a latent class model in mtuem: subclass ``LatentClassMTUEM``,
add the fields you need and override ``class_loglik``. Everything else (EM, BFGS, standard
errors, posteriors, values of time, saving) keeps working.

Run it from mtuem-py with ``uv run python examples/activity_participation.py``: it fits the
2015 model at penalization 5 on the data of the article's fixture, as the R script does (EM
from the start values, then BFGS), and prints the estimation report. EM takes about 17
iterations and most of the run time; BFGS then needs a few iterations.
Tests import it as a module (``tests/test_activity_participation.py`` adds examples/ to
``sys.path``).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
from pydantic import model_validator

import mtuem
from mtuem import Beta, Control, LatentClassMTUEM

FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "reference" / "article"


class ActivityParticipationMTUEM(LatentClassMTUEM):
    """A latent class MTUEM with the participation penalty of the article.

    ``participation`` (a field of ``LatentClassMTUEM``) has one row per class over
    ``free_times``: ``True`` when the class takes part in that free time. A class that does not
    take part in a free time has no equation for it, and the model fixes its elasticity and its
    covariance terms at 0 (``model.beta()`` locks them).

    Dropping the equation alone is not enough. Without it, the class likelihood covers fewer
    equations, so it would happily take individuals who spend many hours on the dropped activity.
    The article multiplies the class likelihood by ``exp(-penalization * hours)``, where
    ``hours`` is the time the individual spends on the activities the class does not do. A
    person with 0 hours on them is not penalized; each hour costs a factor ``exp(-penalization)``.

    ``likelihood_floor`` is the optional floor of the R script: when positive, the class
    likelihood is ``max(likelihood, likelihood_floor)``, so a heavy penalty cannot make it
    underflow to 0 for every class at once.
    """

    penalization: float
    likelihood_floor: float = 0.0

    @model_validator(mode="after")
    def _check_penalty(self) -> ActivityParticipationMTUEM:
        if self.penalization < 0:
            raise ValueError(f"penalization must not be negative, not {self.penalization}")
        if self.likelihood_floor < 0:
            raise ValueError(f"likelihood_floor must not be negative, not {self.likelihood_floor}")
        return self

    def class_loglik(self, s, params, data):
        """Per-row log-likelihood of class ``s``: the standard one, penalized and floored.

        ``super().class_loglik`` is the MTUEM log-likelihood over the equations of the free times
        class ``s`` takes part in (``params`` holds the class parameters under their base names:
        ``Phi``, ``th1``, ``chk21``). mtuem works in logs, so the R product
        ``likelihood * exp(-penalization * hours)`` becomes a subtraction:
        ``log(likelihood) - penalization * hours``. Working in logs also means a heavily
        penalized row gives a large negative number instead of a density that underflows to 0.

        The code runs on NumPy or, compiled, on JAX: ``data["Tf2"]`` is a column and
        ``data[["Tf2", "Tf5"]]`` a matrix of the backend's arrays, and ``data.xp`` is ``numpy``
        or ``jax.numpy``. The Python ``if`` statements below depend on the model, not on
        parameter values, so they are fine under ``jax.jit``.
        """
        ll = super().class_loglik(s, params, data)

        # The free times this class does not take part in, from the participation row.
        absent = self.absent_times(s)
        if absent:
            # Hours on those activities, row by row, then the penalty in logs.
            hours = data[absent].sum(axis=1)
            ll = ll - self.penalization * hours

        # R's pmax(P_s, floor) on densities is max(log P_s, log floor) in logs. It applies to
        # every class, with or without absent activities, as in the R script.
        if self.likelihood_floor > 0:
            ll = data.xp.maximum(ll, math.log(self.likelihood_floor))
        return ll


def article_model(
    work_times: str | list[str],
    free_times: list[str],
    participation: list[list[bool]],
    allocation: list[str],
    penalization: float,
    likelihood_floor: float = 0.0,
) -> ActivityParticipationMTUEM:
    """The model of the article on the ENUT data (columns ``Tc``, ``Ec``, ``w``, a 168 hour week).

    One class per row of ``participation``; the error covariance is a Cholesky factor per class
    (``chk{i}{j}_s``, equations numbered 1 = work, then ``free_times``); class membership is a
    logit on the ``allocation`` covariates with class 1 as the reference. Work time is optimal
    (``optimal_tw=True``), as in the R script.
    """
    return ActivityParticipationMTUEM(
        work_times=[work_times] if isinstance(work_times, str) else list(work_times),
        free_times=list(free_times),
        Tc="Tc", Ec="Ec", w="w", tau=168,
        covariance="cholesky",
        n_classes=len(participation),
        participation=[list(row) for row in participation],
        allocation=list(allocation),
        penalization=penalization,
        likelihood_floor=likelihood_floor,
    )


def article_beta(model: ActivityParticipationMTUEM) -> Beta:
    """The start values and fixed parameters of the article's R scripts.

    ``model.beta()`` already has the structure: ``Theta_s`` fixed at 1, the Cholesky diagonal at
    10 and the rest at 0, the parameters of absent activities fixed at 0 and locked, and the
    allocation of class 1 (the reference) fixed at 0. Only the elasticities change: a class with
    hobbies (the second free time, Tf2) starts from the hobbies class of the two class model, a
    class without them from its no hobbies class. This is the rule of hobbies_rest.R, and it
    gives the start values of hobbies_cholesky.R for two classes.
    """
    if not 2 <= len(model.free_times) <= 4:
        raise ValueError("the article's start values cover two to four free times (Tf1, Tf2, Tf3 and Tf5)")
    beta = model.beta()
    for s in range(1, model.n_classes + 1):
        hobbies = model.participation is None or model.participation[s - 1][1]
        start = {
            "Phi": 0.9 if hobbies else 0.8,
            "thw": -1.0 if hobbies else -0.5,
            "th1": 0.3 if hobbies else 0.2,
            "th2": 0.1,
            "th3": 0.1 if hobbies else 0.5,
            "th4": 0.2,
        }
        # Locked parameters (absent activities) stay at 0; th4 exists only with four free times.
        names = {f"{n}_{s}": v for n, v in start.items()}
        beta.set(**{n: v for n, v in names.items() if n in beta and n not in beta.locked})
    return beta


def from_fixture(path: str | Path) -> tuple[ActivityParticipationMTUEM, Beta, pd.DataFrame, Control]:
    """The model, start values, data and control of a fixture of the article.

    ``path`` is a json written by tests/reference/make_article_reference.R; the data is the
    parquet file next to it. The start values and fixed parameters are those of the R fit
    (``beta_start`` and ``fixed``); the locks come from the model.
    """
    path = Path(path)
    ref = json.loads(path.read_text(encoding="utf-8"))
    expected = {"Tc": "Tc", "Ec": "Ec", "w": "w", "tau": 168, "optimal_tw": True}
    differ = [k for k, v in expected.items() if ref[k] != v]
    if differ:
        raise ValueError(f"the fixture differs from the article model in: {', '.join(differ)}")
    free_times = ref["free_times"]
    model = article_model(
        work_times=ref["work_times"],
        free_times=free_times,
        participation=[[mix[t] for t in free_times] for mix in ref["act_mix"]],
        allocation=ref["allocation"],
        penalization=ref["penalization"],
        likelihood_floor=ref["class_likelihood_floor"],
    )

    beta = model.beta()
    if set(beta.names) != set(ref["beta_start"]):
        raise ValueError("the fixture's start values do not name the parameters of the model")
    for name, reason in beta.locked.items():
        if name not in ref["fixed"] or ref["beta_start"][name] != 0:
            raise ValueError(f"the fixture must fix {name} at 0: {reason}")
    beta.set(**{n: v for n, v in ref["beta_start"].items() if n not in beta.locked})
    beta.free([n for n in beta.fixed_names if n not in ref["fixed"]])
    beta.fix(ref["fixed"])

    data = pd.read_parquet(path.with_suffix(".parquet"))
    control = Control(model_name=ref["model_name"], indiv_id=ref["indiv_id"], output_dir=Path("output"))
    return model, beta, data, control


if __name__ == "__main__":
    ### The 2015 model at penalization 5: two classes over Tf1, Tf2 (hobbies) and Tf3; class 1
    ### does not take part in Tf2. The fixture holds the data and the start values of the R run.
    model, beta, data, control = from_fixture(FIXTURES / "2015_p5.json")
    print(f"Classes: {model.n_classes}; class 1 does not take part in {', '.join(model.absent_times(1))}")
    print(f"Locked at 0: {', '.join(beta.locked)}")

    ### As in R: EM from the start values (at most 100 iterations, tolerance 1e-5), then BFGS
    ### from the EM estimates (at most 500 iterations).
    settings = mtuem.EstimationSettings(
        method="em", em_max_iterations=100, em_tolerance=1e-5, max_iterations=500,
    )
    result = model.fit(data, beta, control, settings)
    result.summary()
    print(result.class_shares())
