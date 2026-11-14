"""Wall-clock cost of the Python methods on the SUpDEq subset (one laptop CPU).

Each method starts from the measured context HRTFs and ends with the HRTFs at all targets; loading the
data and the layout mapping are shared by all methods and not timed. SUpDEq times come from run_supdeq.m.
"""
import json
import time

import numpy as np
import pandas as pd

from align import advance, toa
from mca import BUDGETS, LAYOUTS, SUBJECTS
from reconstruct import TOA_ORDER, prepare, setup
from sh import n_fit, sh_fit, sh_matrix
from sofa import load_subject
from sphere import angles, voronoi_weights

REPEATS = 5


def nn(d, c, t):
    return d["H"][c[np.argmin(angles(d["xyz"][t], d["xyz"][c]), axis=1)]]


def sh(d, c, t, N, lam=0.0, aligned=False):
    w = voronoi_weights(d["xyz"][c])
    H_c = d["H"][c]
    if aligned:
        tau = toa(d["ir"][c], d["fs"])  # onsets of the measured context HRIRs only
        tau_t = sh_matrix(TOA_ORDER, d["az"][t], d["el"][t]) @ sh_fit(
            sh_matrix(TOA_ORDER, d["az"][c], d["el"][c]), tau, w)
        H_c = advance(H_c, d["f"], tau)
    estimate = np.tensordot(sh_matrix(N, d["az"][t], d["el"][t]),
                            sh_fit(sh_matrix(N, d["az"][c], d["el"][c]), H_c, w, lam), axes=1)
    return advance(estimate, d["f"], -tau_t) if aligned else estimate


def seconds(run):
    times = []
    for _ in range(REPEATS):
        start = time.perf_counter()
        run()
        times.append(time.perf_counter() - start)
    return np.median(times)


frozen = json.loads(open("results/frozen_config.json").read())["aligned_ridge_sh"]
rows = []
for database, subject in SUBJECTS:
    d = prepare(load_subject(database, subject))
    for M in BUDGETS:
        for layout in LAYOUTS:
            s = setup(d, layout, M)
            c, t = s["context"], s["targets"]
            N, lam = frozen[layout][str(M)]["N"], frozen[layout][str(M)]["lambda"]
            runs = {"NN": lambda: nn(d, c, t),
                    "SH": lambda: sh(d, c, t, n_fit(M)),
                    "SH-aligned": lambda: sh(d, c, t, n_fit(M), aligned=True),
                    "SH-aligned-ridge": lambda: sh(d, c, t, N, lam, aligned=True)}
            for method, run in runs.items():
                rows.append(dict(database=database, subject=str(subject), layout=layout, M=M,
                                 targets=len(t), method=method, seconds=seconds(run)))
    print(database, subject, "done")
table = pd.DataFrame(rows)
table.to_csv("results/timing.csv", index=False)
print(table.groupby(["method", "M"]).seconds.median().mul(1e3).round(2).unstack("M").to_string())
