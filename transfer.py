"""Compare the error relationships of HUTUBS held-out and SADIE II: slopes with subject-bootstrap intervals and the NN/ridge ranking."""
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path("results")
BOOTSTRAP = 2000
rng = np.random.default_rng(0)
rng_sadie = np.random.default_rng(1)  # separate stream, so the SADIE II intervals leave the HUTUBS ones unchanged

summary = pd.read_csv(RESULTS / "eval_summary.csv")
targets = pd.read_csv(RESULTS / "eval_targets.csv.gz", dtype={"subject": str})


def kappa_slope(df):
    g = df.groupby(["layout", "M", "N"])[["kappa", "median_200-16k"]].mean()
    g = g[g.kappa > 10]
    return np.polyfit(np.log10(g.kappa), g["median_200-16k"], 1)[0]


def bootstrap(df, stat, generator=rng):
    groups = [g for _, g in df.groupby("subject")]
    values = [stat(pd.concat([groups[i] for i in generator.integers(len(groups), size=len(groups))]))
              for _ in range(BOOTSTRAP)]
    return np.percentile(values, [2.5, 97.5])


rows = []
for method in ["SH", "SH-aligned"]:
    df = summary[summary.method == method]
    hutubs, sadie = df[df.database == "hutubs"], df[df.database == "sadie"]
    low, high = bootstrap(hutubs, kappa_slope)
    sadie_low, sadie_high = bootstrap(sadie, kappa_slope, rng_sadie)
    rows.append(dict(relationship="kappa slope (dB per decade)", method=method, hutubs=kappa_slope(hutubs),
                     ci_low=low, ci_high=high, sadie=kappa_slope(sadie), sadie_ci_low=sadie_low, sadie_ci_high=sadie_high))

near = targets[targets.d_min <= 35]
for method in ["NN", "SH-aligned-ridge"]:
    slopes = {}
    for database in ["hutubs", "sadie"]:
        df = near[(near.method == method) & (near.database == database)]
        slopes[database] = df.groupby("subject")[["d_min", "200-16k"]].apply(
            lambda g: 10 * np.polyfit(g.d_min, g["200-16k"], 1)[0])
    h, a = slopes["hutubs"].values, slopes["sadie"].values
    low, high = np.percentile([rng.choice(h, len(h)).mean() for _ in range(BOOTSTRAP)], [2.5, 97.5])
    sadie_low, sadie_high = np.percentile([rng_sadie.choice(a, len(a)).mean() for _ in range(BOOTSTRAP)], [2.5, 97.5])
    rows.append(dict(relationship="d_min slope (dB per 10 degrees)", method=method, hutubs=h.mean(),
                     ci_low=low, ci_high=high, sadie=a.mean(), sadie_ci_low=sadie_low, sadie_ci_high=sadie_high))

per_subject = summary[summary.method.isin(["NN", "SH-aligned-ridge"])]
medians = per_subject.groupby(["database", "layout", "M", "method"])["median_200-16k"].median().unstack("method")
ridge_wins = (medians["SH-aligned-ridge"] < medians["NN"]).unstack("database")
agree = int((ridge_wins.hutubs == ridge_wins.sadie).sum())
rows.append(dict(relationship="ranking vs NN (conditions agreeing of 12)", method="SH-aligned-ridge",
                 hutubs=int(ridge_wins.hutubs.sum()), sadie=int(ridge_wins.sadie.sum()), agree=agree))

table = pd.DataFrame(rows)
slope = table.relationship.str.contains("slope")
table.loc[slope, "change_pct"] = 100 * (table.sadie[slope] / table.hutubs[slope] - 1)
table.to_csv(RESULTS / "transfer.csv", index=False)
print(table.round(2).to_string(index=False))
print("\nConditions where time-aligned ridge SH beats NN (median over subjects):")
print(ridge_wins.to_string())
