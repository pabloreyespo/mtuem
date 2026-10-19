"""The printed report and the saved result files."""

import getpass
import io
import math
import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd

from .estimate import Estimate

try:
    __version__ = version("mtuem")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "dev"


def _r_num(x, digits: int = 7) -> str:
    """Format a number the way R's ``cat`` prints it (up to 7 significant digits)."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, (float, np.floating)) and math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    s = f"{x:.{digits}g}"
    if "e" in s:
        mant, exp = s.split("e")
        s = f"{mant}e{exp[0]}{exp[1:].zfill(2)}"
    return s


def _format_column(values: np.ndarray, n_signif: int = 4, width: int = 11) -> list[str]:
    """Port of the column formatting in ``apollo_print``."""
    out, scie, ndec = [], [], []
    for v in values:
        if v is None or not np.isfinite(v):
            m = "NA".rjust(width) if (v is None or np.isnan(v)) else ("Inf" if v > 0 else "-Inf").rjust(width)
        else:
            m = f"%#{width}.{n_signif}g" % v
        is_sci = "e" in m
        nd = len(m) - m.index(".") - 1 if "." in m else 0
        if not is_sci and nd > 6 and np.isfinite(v):
            m = f"%{width}.{n_signif}e" % v
            is_sci = True
        out.append(m)
        scie.append(is_sci)
        ndec.append(nd)
    non_sci = [d for d, s in zip(ndec, scie) if not s]
    obj = min(max(non_sci), 6) if non_sci else 6
    for i, v in enumerate(values):
        if not scie[i] and ndec[i] < obj and v is not None and np.isfinite(v):
            out[i] = f"%{width}.{obj}f" % v
    w = max(len(m) for m in out)
    return [m.rjust(w) for m in out]


def print_matrix(M: pd.DataFrame, n_signif: int = 4, width: int = 11) -> str:
    """Render a numeric table the way R prints a data.frame of formatted numbers."""
    cols = [_format_column(M[c].to_numpy(dtype=float), n_signif, width) for c in M.columns]
    rn = [str(r) for r in M.index]
    rw = max(len(r) for r in rn) if rn else 0
    widths = [max(len(str(c)), len(col[0]) if col else 0) for c, col in zip(M.columns, cols)]
    lines = [" " * rw + "".join(" " + str(c).rjust(w) for c, w in zip(M.columns, widths))]
    for i, r in enumerate(rn):
        lines.append(r.ljust(rw) + "".join(" " + col[i].rjust(w) for col, w in zip(cols, widths)))
    return "\n".join(lines)


def _hms(t: float) -> str:
    h = int(t // 3600)
    m = int((t - h * 3600) // 60)
    s = round(t - h * 3600 - m * 60, 2)
    return f"{h:02d}:{m:02d}:{_r_num(s)}"


def _estimates_table(est: Estimate) -> pd.DataFrame:
    with np.errstate(invalid="ignore", divide="ignore"):
        return pd.DataFrame(
            {
                "Estimate": est.estimate,
                "s.e.": est.se,
                "t.rat.(0)": est.estimate / est.se,
                "Rob.s.e.": est.robse,
                "Rob.t.rat.(0)": est.estimate / est.robse,
            }
        )


def _signif(x: float, digits: int) -> float:
    if x == 0 or not np.isfinite(x):
        return x
    return round(x, digits - int(math.floor(math.log10(abs(x)))) - 1)


def output_text(est: Estimate, print_covar: bool = False, print_corr: bool = False) -> str:
    """The text report of an estimation."""
    buf = io.StringIO()
    w = buf.write
    ctrl = est.control
    w(f"mtuem {__version__}: Microeconomic Time-Use-Expenditure Models\n")
    w(f"Model run by {getpass.getuser()} on Python {platform.python_version()} for {platform.system()}.\n\n")
    w(f"Model name                                  : {ctrl.model_name}\n")
    w(f"Model description                           : {ctrl.description}\n")
    w(f"Model run at                                : {est.startTime}\n")
    w(f"Estimation method                           : {est.routine}\n")
    w(f"Estimation diagnosis                        : {est.message}\n")
    eig = est.hessianEigenValue
    if eig is None:
        w("Optimisation diagnosis                      : Unknown (no Hessian requested by user)\n")
    else:
        if np.iscomplexobj(eig):
            txt = ["Non-symmetrical hessian", "Complex eigenvalues", "NA"]
        else:
            if np.all(eig < 0):
                txt = ["Maximum found", "Negative definite"]
            elif np.any(eig == 0):
                txt = ["Inconclusive test", "Some eigenvalues are zero"]
            elif np.all(eig > 0):
                txt = ["Minimum found", "Positive definite"]
            else:
                txt = ["Saddle point found", "Some eigenvalues are positive and others negative"]
            txt.append(_r_num(round(float(np.max(eig)), 6), 15))
        if not np.iscomplexobj(eig) and np.all(eig < 0):
            txt.append(_r_num(_signif(float(np.max(eig) / np.min(eig)), 6), 15))
        else:
            txt.append("not calculated (Hessian is not negative definite)")
        w(f"Optimisation diagnosis                      : {txt[0]}\n")
        w(f"     hessian properties                     : {txt[1]}\n")
        w(f"     maximum eigenvalue                     : {txt[2]}\n")
        w(f"     reciprocal of condition number         : {txt[3]}\n")
    w(f"Number of individuals                       : {est.nIndivs}\n")
    w(f"Number of rows in database                  : {est.nObs}\n")
    w(f"Number of modelled outcomes                 : {est.nObs}\n\n")
    w("Model without mixing\n\n")
    multi = len(est.LLout) > 1
    w(f"LL(start)                                   : {_r_num(round(est.LLStart, 2))}\n")
    if multi:
        w("LL (whole model) at equal shares, LL(0)     : NA\n")
        w("LL (whole model) at observed shares, LL(C)  : NA\n")
        w(f"LL(final, whole model)                      : {_r_num(round(est.maximum, 2))}\n")
    else:
        w("LL at equal shares, LL(0)                   : NA\n")
        w("LL at observed shares, LL(C)                : NA\n")
        w(f"LL(final)                                   : {_r_num(round(est.maximum, 2))}\n")
    for label in ("Rho-squared vs equal shares", "Adj.Rho-squared vs equal shares",
                  "Rho-squared vs observed shares", "Adj.Rho-squared vs observed shares"):
        w(f"{label.ljust(45)}:  Not applicable \n")
    w(f"AIC                                         :  {_r_num(round(est.AIC, 2))} \n")
    w(f"BIC                                         :  {_r_num(round(est.BIC, 2))} \n")
    if multi:
        w("\n")
        for nam, val in list(est.LLout.items())[1:]:
            sp1 = " " * (33 - len(nam) - 6) if 6 + len(nam) < 33 else ""
            sp2 = " " * (33 - len(nam) - 10) if 9 + len(nam) < 33 else ""
            w(f"LL(0,{nam}){sp1}: Not applicable\n")
            w(f"LL(final,{nam}){sp2}: {_r_num(round(val, 2))}\n")
    w("\n")
    w(f"Estimated parameters                        : {est.nFreeParams}\n")
    w(f"Time taken (hh:mm:ss)                       :  {_hms(est.timeTaken)} \n")
    w(f"     pre-estimation                         :  {_hms(est.timePre)} \n")
    w(f"     estimation                             :  {_hms(est.timeEst)} \n")
    w(f"     post-estimation                        :  {_hms(est.timePost)} \n")
    tmp = "" if est.successfulEstimation else f"({est.message})"
    w(f"Iterations                                  :  {est.nIter} {tmp} \n\n")
    w("Unconstrained optimisation.\n\n")
    w("Estimates:\n")
    w(print_matrix(_estimates_table(est)) + "\n\n")
    if est.classShares is not None:
        w("\nSummary of class allocation for model component :\n")
        shares = pd.DataFrame({"Mean prob.": est.classShares.round(4)})
        w(shares.to_string() + "\n\n")
    if print_covar:
        w("\nClassical covariance matrix:\n" + print_matrix(est.varcov) + "\n")
        w("\nRobust covariance matrix:\n" + print_matrix(est.robvarcov) + "\n")
    if print_corr:
        w("\nClassical correlation matrix:\n" + print_matrix(est.corrmat) + "\n")
        w("\nRobust correlation matrix:\n" + print_matrix(est.robcorrmat) + "\n")
    return buf.getvalue()


def save_files(est: Estimate, folder: Path, name: str) -> dict[str, Path]:
    """Write ``<name>_output.txt``, ``<name>_estimates.csv`` and the covariance and correlation csv files."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    files = {}
    files["output"] = folder / f"{name}_output.txt"
    files["output"].write_text(output_text(est, True, True), encoding="utf-8")
    table = _estimates_table(est)
    table.columns = ["Estimate", "Std.err.", "t-ratio(0)", "Rob.std.err.", "Rob.t-ratio(0)"]
    files["estimates"] = folder / f"{name}_estimates.csv"
    table.to_csv(files["estimates"], na_rep="NA")
    for key, mat in (("covar", est.varcov), ("robcovar", est.robvarcov), ("corr", est.corrmat), ("robcorr", est.robcorrmat)):
        files[key] = folder / f"{name}_{key}.csv"
        mat.to_csv(files[key], na_rep="NA")
    return files
