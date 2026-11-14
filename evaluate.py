"""Apply the frozen configuration, unchanged, to the HUTUBS held-out half and to SADIE II."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from reconstruct import BUDGETS, LAYOUTS, cue_errors, ear_lsd, lsd_bands, prepare, setup, sh_estimate
from sh import n_fit, n_max
from sofa import HUTUBS_TEST, SADIE_SUBJECTS, load_subject

RESULTS = Path("results")


def evaluate_subject(database, subject, config):
    d = prepare(load_subject(database, subject))
    summary, targets = [], []
    for M in BUDGETS:
        for layout in LAYOUTS:
            s = setup(d, layout, M)
            frozen = config[layout][str(M)]
            runs = {("NN", 0, 0.0): (d["H"][s["nearest"]], np.nan, np.nan)}
            for N in range(1, n_max(M) + 1):
                runs[("SH", N, 0.0)] = sh_estimate(d, s, N)
                runs[("SH-aligned", N, 0.0)] = sh_estimate(d, s, N, aligned=True)
            runs[("SH-aligned-ridge", frozen["N"], frozen["lambda"])] = sh_estimate(
                d, s, frozen["N"], frozen["lambda"], aligned=True)
            # follow-up: alignment with the smoothed onsets of the context directions, same orders and lambda
            runs[("SH-aligned-smooth", n_fit(M), 0.0)] = sh_estimate(d, s, n_fit(M), aligned=True, smooth=True)
            runs[("SH-aligned-ridge-smooth", frozen["N"], frozen["lambda"])] = sh_estimate(
                d, s, frozen["N"], frozen["lambda"], aligned=True, smooth=True)

            for (method, N, lam), (estimate, kappa, coef_norm) in runs.items():
                errors = lsd_bands(d, s, estimate)
                energy_db = 10 * np.log10(np.sum(np.abs(estimate) ** 2) / np.sum(np.abs(d["H"][s["targets"]]) ** 2))
                summary.append(dict(database=database, subject=subject, layout=layout, M=M, method=method,
                                    N=N, lam=lam, kappa=kappa, coef_norm=coef_norm, energy_db=energy_db,
                                    finite=bool(np.isfinite(estimate).all()),
                                    **{f"median_{band}": np.median(e) for band, e in errors.items()}))
                if method in ("SH", "SH-aligned") and N != n_fit(M):
                    continue  # per-target rows only for the orders compared across methods
                itd_error, ild_error = cue_errors(d, s, estimate)
                targets.append(pd.DataFrame(dict(
                    database=database, subject=subject, layout=layout, M=M, method=method,
                    target=s["targets"], az=d["az"][s["targets"]], el=d["el"][s["targets"]],
                    d_min=s["d_min"], itd_error_us=itd_error, ild_error_db=ild_error, **errors,
                    **ear_lsd(d, s, estimate))))
    return summary, targets


def main():
    config = json.loads((RESULTS / "frozen_config.json").read_text())["aligned_ridge_sh"]
    summary, targets = [], []
    for database, subjects in [("hutubs", HUTUBS_TEST), ("sadie", SADIE_SUBJECTS)]:
        for subject in subjects:
            rows, frames = evaluate_subject(database, subject, config)
            summary += rows
            targets += frames
            print(database, subject, "done")
    pd.DataFrame(summary).to_csv(RESULTS / "eval_summary.csv", index=False)
    pd.concat(targets).astype({"subject": str}).to_csv(RESULTS / "eval_targets.csv.gz", index=False,
                                                       float_format="%.6g")


if __name__ == "__main__":
    main()
