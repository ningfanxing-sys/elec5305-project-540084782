"""Choose the SH order and ridge parameter of time-aligned SH on the HUTUBS selection half, then freeze."""
import datetime
import json
from pathlib import Path

import numpy as np
import pandas as pd

from reconstruct import BUDGETS, LAYOUTS, TOA_ORDER, lsd_bands, prepare, setup, sh_estimate
from sh import n_max
from sofa import HUTUBS_SELECT, HUTUBS_TEST, SADIE_SUBJECTS, load_subject
from sphere import EL_MIN

LAMBDAS = [0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
RESULTS = Path("results")

TRANSFER_CRITERIA = [
    "kappa: for ordinary and time-aligned SH (lambda = 0), fit median LSD (200 Hz-16 kHz, mean over subjects) "
    "against log10(kappa) over all (layout, M, N) with kappa > 10. The relationship transfers if the SADIE II "
    "slope lies inside the 95% bootstrap interval (subjects resampled, 2000 times) of the HUTUBS held-out slope.",
    "d_min: for NN and the frozen time-aligned ridge SH, fit target LSD against d_min for targets with "
    "d_min <= 35 degrees, separately per subject. The relationship transfers if the mean SADIE II slope lies "
    "inside the 95% bootstrap interval of the mean HUTUBS held-out slope.",
    "ranking: for each of the 12 (layout, M) conditions, compare the frozen time-aligned ridge SH with NN by "
    "the median over subjects of the per-subject median LSD. The ranking transfers if it agrees between HUTUBS "
    "held-out and SADIE II in at least 10 of 12 conditions.",
]


def main():
    rows = []
    for subject in HUTUBS_SELECT:
        d = prepare(load_subject("hutubs", subject))
        for M in BUDGETS:
            for layout in LAYOUTS:
                s = setup(d, layout, M)
                for N in range(1, n_max(M) + 1):
                    for lam in LAMBDAS:
                        estimate, kappa, _ = sh_estimate(d, s, N, lam, aligned=True)
                        rows.append(dict(subject=subject, layout=layout, M=M, N=N, lam=lam, kappa=kappa,
                                         median_lsd=np.median(lsd_bands(d, s, estimate)["200-16k"])))
        print("selection subject", subject, "done")
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "selection.csv", index=False)

    mean = df.groupby(["layout", "M", "N", "lam"]).median_lsd.mean().reset_index()
    best = mean.loc[mean.groupby(["layout", "M"]).median_lsd.idxmin()]
    config = {
        "frozen_on": datetime.date.today().isoformat(),
        "selection_subjects": [f"pp{s}" for s in HUTUBS_SELECT],
        "held_out_subjects": [f"pp{s}" for s in HUTUBS_TEST],
        "external_subjects": SADIE_SUBJECTS,
        "toa": "4th-order 1 kHz high-pass, first sample above 20% of the peak, 8x upsampling",
        "toa_order": TOA_ORDER,
        "directions_used": f"elevation >= {EL_MIN} degrees in both databases",
        "lambda_grid": LAMBDAS,
        "aligned_ridge_sh": {layout: {str(r.M): {"N": int(r.N), "lambda": float(r.lam)}
                                      for r in best[best.layout == layout].itertuples()} for layout in LAYOUTS},
        "transfer_criteria": TRANSFER_CRITERIA,
    }
    (RESULTS / "frozen_config.json").write_text(json.dumps(config, indent=2) + "\n")
    print(best.to_string(index=False))


if __name__ == "__main__":
    main()
