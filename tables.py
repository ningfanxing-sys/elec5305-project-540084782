"""Print the numbers behind the report's tables (all values are medians unless stated)."""
import numpy as np
import pandas as pd

pd.set_option("display.width", 160)
targets = pd.read_csv("results/eval_targets.csv.gz", dtype={"subject": str})
summary = pd.read_csv("results/eval_summary.csv", dtype={"subject": str})
mca = pd.read_csv("results/mca_targets.csv.gz", dtype={"subject": str})
timing = pd.read_csv("results/timing.csv", dtype={"subject": str})
supdeq = pd.read_csv("results/supdeq_times.csv")
COLUMNS = {"200-16k": "LSD", "8-16k": "LSD 8-16k", "itd_error_us": "ITD us", "ild_error_db": "ILD dB"}

print("Overall comparison, all layouts and budgets (follow-up rows: smoothed onsets)")
order = ["NN", "SH", "SH-aligned", "SH-aligned-ridge", "SH-aligned-smooth", "SH-aligned-ridge-smooth"]
table = targets.groupby(["database", "method"])[list(COLUMNS)].median().rename(columns=COLUMNS)
print(table.reindex(order, level="method").round(2).to_string(), "\n")

print("Frozen ridge against lambda = 0 at the same order (mean of per-subject medians, dB)")
ridge = summary[summary.method == "SH-aligned-ridge"][["database", "subject", "layout", "M", "N", "median_200-16k"]]
plain = summary[(summary.method == "SH-aligned")][["database", "subject", "layout", "M", "N", "median_200-16k"]]
pair = ridge.merge(plain, on=["database", "subject", "layout", "M", "N"], suffixes=("_ridge", "_plain"))
print((pair["median_200-16k_plain"] - pair["median_200-16k_ridge"]).groupby(pair.database).mean().round(2).to_string(), "\n")

print("Matched distance d_min (all layouts and budgets)")
binned = targets.assign(bin=pd.cut(targets.d_min, [5, 10, 15, 20, 30]))
print(binned.groupby(["bin", "method", "database"], observed=True)["200-16k"].median().unstack(["method", "database"])
      [["NN", "SH", "SH-aligned", "SH-aligned-ridge"]].round(1).to_string(), "\n")

print("SUpDEq subset, identical targets")
subset = targets.merge(mca[["database", "subject", "layout", "M"]].drop_duplicates())
both = pd.concat([subset[subset.method.isin(["NN", "SH-aligned", "SH-aligned-ridge"])], mca])
lateral = np.rad2deg(np.arcsin(np.cos(np.deg2rad(both.el)) * np.sin(np.deg2rad(both.az))))
ears = pd.concat([both.assign(rel=lateral, ear=both.lsd_left), both.assign(rel=-lateral, ear=both.lsd_right)])
mca_table = both.groupby(["database", "method"])[list(COLUMNS)].median().rename(columns=COLUMNS)
mca_table["contra LSD"] = ears[ears.rel < -45].groupby(["database", "method"]).ear.median()
mca_table["ipsi LSD"] = ears[ears.rel > 45].groupby(["database", "method"]).ear.median()
far = both[(both.d_min >= 25) & (both.d_min < 35)]
mca_table["LSD d25-35"] = far.groupby(["database", "method"])["200-16k"].median()
order = ["NN", "SH-aligned", "SUpDEq-OBTA", "SUpDEq-MCA", "SH-aligned-ridge"]
print(mca_table.reindex(order, level="method").round(2).to_string(), "\n")

print("Run time per reconstruction (median, ms; Python on the SUpDEq subset; MATLAB SUpDEq in s)")
print(timing.groupby(["database", "method"]).seconds.median().mul(1e3).unstack("database").round(1).to_string())
supdeq = supdeq[~supdeq.file.str.startswith("demo")]
supdeq["database"] = supdeq.file.str.split("_").str[0]
print(supdeq.groupby("database")[["obta_s", "mca_s"]].median().round(2).to_string())
print("targets per reconstruction:", timing.groupby("database").targets.median().to_dict())
