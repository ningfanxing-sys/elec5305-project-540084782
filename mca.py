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

from reconstruct import cue_errors, lsd_bands, prepare, setup
from sh import n_fit
from sofa import HUTUBS_TEST, SADIE_SUBJECTS, load_subject

SUBJECTS = [("hutubs", s) for s in HUTUBS_TEST[:5]] + [("sadie", s) for s in SADIE_SUBJECTS[::4]]
LAYOUTS = ["uniform", "gap"]
BUDGETS = [25, 50, 100]
FOLDER = Path("results/mca")


def conditions():
    for database, subject in SUBJECTS:
        d = prepare(load_subject(database, subject))
        for M in BUDGETS:
            for layout in LAYOUTS:
                yield database, subject, layout, M, d, setup(d, layout, M)


def export():
    (FOLDER / "in").mkdir(parents=True, exist_ok=True)
    for database, subject, layout, M, d, s in conditions():
        c, t = s["context"], s["targets"]
        # SUpDEq grids are [azimuth, colatitude] in degrees
        savemat(FOLDER / "in" / f"{database}_{subject}_{layout}_{M}.mat", dict(
            HRTF_L=d["H"][c, 0], HRTF_R=d["H"][c, 1], f=d["f"], Nmax=n_fit(M),
            samplingGrid=np.column_stack([d["az"][c], 90 - d["el"][c]]),
            ipGrid=np.column_stack([d["az"][t], 90 - d["el"][t]])))


def score():
    frames = []
    for database, subject, layout, M, d, s in conditions():
        out = loadmat(FOLDER / "out" / f"{database}_{subject}_{layout}_{M}.mat")
        for method, key in [("SUpDEq-OBTA", "off"), ("SUpDEq-MCA", "on")]:
            estimate = np.stack([out[f"HRTF_L_{key}"], out[f"HRTF_R_{key}"]], axis=1)
            itd_error, ild_error = cue_errors(d, s, estimate)
            frames.append(pd.DataFrame(dict(
                database=database, subject=str(subject), layout=layout, M=M, method=method,
                target=s["targets"], d_min=s["d_min"], itd_error_us=itd_error, ild_error_db=ild_error,
                **lsd_bands(d, s, estimate))))
    pd.concat(frames).to_csv("results/mca_targets.csv.gz", index=False, float_format="%.6g")


if __name__ == "__main__":
    {"export": export, "score": score}[sys.argv[1]]()
