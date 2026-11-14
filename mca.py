"""SUpDEq / MCA reference on a subset.

python mca.py export   -> results/mca/in/*.mat (context HRTFs and target grid for MATLAB)
matlab -batch run_supdeq
python mca.py score    -> results/mca_targets.csv.gz, scored exactly like the Python methods
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat, savemat

from reconstruct import cue_errors, ear_lsd, lsd_bands, prepare, setup
from sh import n_fit
from sofa import HUTUBS_TEST, SADIE_SUBJECTS, load_subject

SUBJECTS = [("hutubs", s) for s in HUTUBS_TEST[:5]] + [("sadie", s) for s in SADIE_SUBJECTS[::4]]
LAYOUTS = ["uniform", "gap"]
BUDGETS = [25, 50, 100]
FOLDER = Path("results/mca")
DEMO_SUBJECTS = [("hutubs", 7), ("sadie", "H3")]
DEMO_BUDGETS = [25, 50]


def conditions():
    for database, subject in SUBJECTS:
        d = prepare(load_subject(database, subject))
        for M in BUDGETS:
            for layout in LAYOUTS:
                yield database, subject, layout, M, d, setup(d, layout, M)


def save_input(name, d, context, points, M):
    # SUpDEq grids are [azimuth, colatitude] in degrees
    savemat(FOLDER / "in" / f"{name}.mat", dict(
        HRTF_L=d["H"][context, 0], HRTF_R=d["H"][context, 1], f=d["f"], Nmax=n_fit(M),
        samplingGrid=np.column_stack([d["az"][context], 90 - d["el"][context]]),
        ipGrid=np.column_stack([d["az"][points], 90 - d["el"][points]])))


def export():
    (FOLDER / "in").mkdir(parents=True, exist_ok=True)
    for database, subject, layout, M, d, s in conditions():
        save_input(f"{database}_{subject}_{layout}_{M}", d, s["context"], s["targets"], M)
    # listening demo and Figure of the horizontal plane: reconstruct every measured direction
    for database, subject in DEMO_SUBJECTS:
        d = prepare(load_subject(database, subject))
        for M in DEMO_BUDGETS:
            save_input(f"demo_{database}_{subject}_uniform_{M}", d, setup(d, "uniform", M)["context"],
                       np.arange(len(d["az"])), M)


def score():
    frames = []
    for database, subject, layout, M, d, s in conditions():
        out = loadmat(FOLDER / "out" / f"{database}_{subject}_{layout}_{M}.mat")
        for method, key in [("SUpDEq-OBTA", "off"), ("SUpDEq-MCA", "on")]:
            estimate = np.stack([out[f"HRTF_L_{key}"], out[f"HRTF_R_{key}"]], axis=1)
            itd_error, ild_error = cue_errors(d, s, estimate)
            frames.append(pd.DataFrame(dict(
                database=database, subject=str(subject), layout=layout, M=M, method=method,
                target=s["targets"], az=d["az"][s["targets"]], el=d["el"][s["targets"]],
                d_min=s["d_min"], itd_error_us=itd_error, ild_error_db=ild_error,
                **lsd_bands(d, s, estimate), **ear_lsd(d, s, estimate))))
    pd.concat(frames).to_csv("results/mca_targets.csv.gz", index=False, float_format="%.6g")


if __name__ == "__main__":
    {"export": export, "score": score}[sys.argv[1]]()
