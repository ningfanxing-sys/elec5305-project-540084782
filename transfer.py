"""Check the transfer criteria written into results/frozen_config.json before SADIE II was evaluated."""
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path("results")
BOOTSTRAP = 2000
rng = np.random.default_rng(0)

summary = pd.read_csv(RESULTS / "eval_summary.csv")
targets = pd.read_csv(RESULTS / "eval_targets.csv.gz", dtype={"subject": str})


def kappa_slope(df):
    g = df.groupby(["layout", "M", "N"])[["kappa", "median_200-16k"]].mean()
    g = g[g.kappa > 10]
    return np.polyfit(np.log10(g.kappa), g["median_200-16k"], 1)[0]


def bootstrap(df, stat):
    groups = [g for _, g in df.groupby("subject")]
    values = [stat(pd.concat([groups[i] for i in rng.integers(len(groups), size=len(groups))]))
              for _ in range(BOOTSTRAP)]
    return np.percentile(values, [2.5, 97.5])


rows = []
for method in ["SH", "SH-aligned"]:
    df = summary[summary.method == method]
    hutubs, sadie = df[df.database == "hutubs"], df[df.database == "sadie"]
    low, high = bootstrap(hutubs, kappa_slope)
    rows.append(dict(criterion="kappa slope (dB per decade)", method=method, hutubs=kappa_slope(hutubs),
                     ci_low=low, ci_high=high, sadie=kappa_slope(sadie)))

near = targets[targets.d_min <= 35]
for method in ["NN", "SH-aligned-ridge"]:
    slopes = {}
    for database in ["hutubs", "sadie"]:
        df = near[(near.method == method) & (near.database == database)]
        slopes[database] = df.groupby("subject")[["d_min", "200-16k"]].apply(
            lambda g: 10 * np.polyfit(g.d_min, g["200-16k"], 1)[0])
    h = slopes["hutubs"].values
    low, high = np.percentile([rng.choice(h, len(h)).mean() for _ in range(BOOTSTRAP)], [2.5, 97.5])
    rows.append(dict(criterion="d_min slope (dB per 10 degrees)", method=method, hutubs=h.mean(),
                     ci_low=low, ci_high=high, sadie=slopes["sadie"].mean()))

per_subject = summary[summary.method.isin(["NN", "SH-aligned-ridge"])]
medians = per_subject.groupby(["database", "layout", "M", "method"])["median_200-16k"].median().unstack("method")
ridge_wins = (medians["SH-aligned-ridge"] < medians["NN"]).unstack("database")
agree = int((ridge_wins.hutubs == ridge_wins.sadie).sum())
rows.append(dict(criterion="ranking vs NN (conditions agreeing of 12)", method="SH-aligned-ridge",
                 hutubs=int(ridge_wins.hutubs.sum()), sadie=int(ridge_wins.sadie.sum()), agree=agree))

table = pd.DataFrame(rows)
table["transfers"] = [(r.ci_low <= r.sadie <= r.ci_high) if "slope" in r.criterion else r.agree >= 10
                      for r in table.itertuples()]
table.to_csv(RESULTS / "transfer.csv", index=False)
print(table.round(2).to_string(index=False))
print("\nConditions where time-aligned ridge SH beats NN (median over subjects):")
print(ridge_wins.to_string())
