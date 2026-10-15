"""Post hoc: how much per-target error do budget, local distance and conditioning explain?

Linear models of target LSD (200 Hz-16 kHz) are fitted on HUTUBS held-out and applied unchanged to SADIE II.
"""
import numpy as np
import pandas as pd

from sh import n_fit

MODELS = {"M only": [], "M + log10 kappa": ["logk"], "M + d_min": ["d", "d2"],
          "M + d_min + log10 kappa": ["d", "d2", "logk"]}

targets = pd.read_csv("results/eval_targets.csv.gz", dtype={"subject": str})
summary = pd.read_csv("results/eval_summary.csv", dtype={"subject": str})

# condition number of the fit that produced each target row
used = summary[((summary.method == "SH") & (summary.N == summary.M.map(n_fit))) |
               (summary.method == "SH-aligned-ridge")]
kappa = used.set_index(["database", "subject", "layout", "M", "method"]).kappa
df = targets.join(kappa, on=["database", "subject", "layout", "M", "method"])
df = df.assign(logk=np.log10(df.kappa), d=df.d_min, d2=df.d_min ** 2, y=df["200-16k"])


def design(d, columns):
    budget = [(d.M == M).astype(float) for M in (25, 50, 100)]
    return np.column_stack(budget + [d[c] for c in columns])


def r2(y, prediction):
    return 1 - np.sum((y - prediction) ** 2) / np.sum((y - y.mean()) ** 2)


rows = []
for method in ["NN", "SH", "SH-aligned-ridge"]:
    hutubs = df[(df.method == method) & (df.database == "hutubs")]
    sadie = df[(df.method == method) & (df.database == "sadie")]
    for name, columns in MODELS.items():
        if method == "NN" and "logk" in columns:
            continue
        beta = np.linalg.lstsq(design(hutubs, columns), hutubs.y, rcond=None)[0]
        own = np.linalg.lstsq(design(sadie, columns), sadie.y, rcond=None)[0]
        rows.append(dict(method=method, model=name,
                         r2_hutubs=r2(hutubs.y, design(hutubs, columns) @ beta),
                         r2_sadie_transferred=r2(sadie.y, design(sadie, columns) @ beta),
                         r2_sadie_refit=r2(sadie.y, design(sadie, columns) @ own)))
table = pd.DataFrame(rows)
table.to_csv("results/prediction.csv", index=False)
print(table.round(2).to_string(index=False))
