import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figstyle import save
from metrics import BANDS

BUDGETS = [25, 50, 100]
LAYOUTS = ["uniform", "random1", "random2", "gap"]
MARKERS = dict(zip(LAYOUTS, ["o", "s", "^", "D"]))
METHODS = {"NN": "#0072B2", "SH": "#E69F00", "SH-aligned": "#009E73", "SH-aligned-ridge": "#CC79A7",
           "SUpDEq-OBTA": "#56B4E9", "SUpDEq-MCA": "#D55E00"}
MAIN = ["NN", "SH", "SH-aligned", "SH-aligned-ridge"]
M_COLORS = dict(zip(BUDGETS, ["#0072B2", "#E69F00", "#009E73"]))
DATABASES = {"hutubs": "HUTUBS held-out", "sadie": "SADIE II"}

selection = pd.read_csv("results/selection.csv")
summary = pd.read_csv("results/eval_summary.csv").rename(columns={"median_200-16k": "lsd"})
targets = pd.read_csv("results/eval_targets.csv.gz", dtype={"subject": str})
mca = pd.read_csv("results/mca_targets.csv.gz", dtype={"subject": str})


def dmin_curve(ax, df, method):
    groups = df.groupby(pd.cut(df.d_min, np.arange(0, 70, 5)), observed=True)["200-16k"]
    medians = groups.median()[groups.size() >= 50]
    ax.plot([b.mid for b in medians.index], medians.values, marker="o", markersize=3,
            color=METHODS[method], label=method)


def m_and_layout_legend(ax, **kwargs):
    handles = [plt.Line2D([], [], color=c, marker="o", linestyle="none", label=f"M = {M}") for M, c in M_COLORS.items()]
    handles += [plt.Line2D([], [], color="gray", marker=m, linestyle="none", label=l) for l, m in MARKERS.items()]
    ax.legend(handles=handles, ncol=2, **kwargs)


def grouped_bars(ax, values, labels):
    """values: DataFrame with methods as rows and groups as columns."""
    x = np.arange(len(labels))
    width = 0.8 / len(values)
    for i, (method, row) in enumerate(values.iterrows()):
        ax.bar(x + (i - (len(values) - 1) / 2) * width, row.values, width, color=METHODS[method], label=method)
    ax.set_xticks(x, labels)
    ax.grid(axis="x", visible=False)


# 1. ridge on the selection half: preferred lambda and its benefit against kappa
mean = selection.groupby(["layout", "M", "N", "lam"])[["kappa", "median_lsd"]].mean().reset_index()
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.6))
for (layout, M, N), g in mean.groupby(["layout", "M", "N"]):
    best = g.loc[g.median_lsd.idxmin()]
    gain = g[g.lam == 0].median_lsd.iloc[0] - best.median_lsd
    axes[0].scatter(np.log10(best.kappa), max(best.lam, 1e-6), color=M_COLORS[M], marker=MARKERS[layout], s=12)
    axes[1].scatter(np.log10(best.kappa), gain, color=M_COLORS[M], marker=MARKERS[layout], s=12)
axes[0].set(yscale="log", xlabel="log$_{10}$ κ", ylabel="Preferred λ (λ = 0 at 10$^{-6}$)",
            title="(a) Preferred ridge parameter")
axes[1].set(xlabel="log$_{10}$ κ", ylabel="LSD reduction from ridge (dB)", title="(b) Benefit of ridge")
m_and_layout_legend(axes[1], loc="upper left")
save(fig, "lambda_vs_kappa")

# 2. error against condition number, HUTUBS held-out against SADIE II
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), sharey=True)
for ax, letter, method in zip(axes, "ab", ["SH", "SH-aligned"]):
    for database, fill, style in [("hutubs", True, "-"), ("sadie", False, "--")]:
        g = summary[(summary.method == method) & (summary.database == database)]
        g = g.groupby(["layout", "M", "N"])[["kappa", "lsd"]].mean().reset_index()
        for r in g.itertuples():
            ax.scatter(np.log10(r.kappa), r.lsd, marker=MARKERS[r.layout], s=11, linewidths=0.7,
                       facecolor=M_COLORS[r.M] if fill else "none", edgecolor=M_COLORS[r.M])
        high = g[g.kappa > 10]
        slope, intercept = np.polyfit(np.log10(high.kappa), high.lsd, 1)
        x = np.linspace(1, np.log10(high.kappa.max()), 10)
        ax.plot(x, slope * x + intercept, style, color="black", linewidth=1,
                label=f"{DATABASES[database]}: {slope:.1f} dB/decade")
    ax.set(title=f"({letter}) {method}, λ = 0", xlabel="log$_{10}$ κ")
    ax.legend(loc="upper left")
axes[0].set_ylabel("Median LSD, 200 Hz–16 kHz (dB)")
handles = [plt.Line2D([], [], color=c, marker="o", linestyle="none", label=f"M = {M}") for M, c in M_COLORS.items()]
handles += [plt.Line2D([], [], color="gray", marker=m, linestyle="none", label=l) for l, m in MARKERS.items()]
handles += [plt.Line2D([], [], color="gray", marker="o", linestyle="none", label="HUTUBS held-out (filled)"),
            plt.Line2D([], [], color="gray", marker="o", markerfacecolor="none", linestyle="none", label="SADIE II (open)")]
fig.legend(handles=handles, loc="outside lower center", ncol=5)
save(fig, "error_vs_kappa")

# 3. error against distance to the nearest context direction
fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.4), sharey=True, sharex="col")
for i, (row, database) in enumerate(zip(axes, DATABASES)):
    for j, (ax, M) in enumerate(zip(row, BUDGETS)):
        for method in MAIN:
            dmin_curve(ax, targets[(targets.database == database) & (targets.M == M) & (targets.method == method)], method)
        ax.set(title=f"({'abcdef'[3 * i + j]}) {DATABASES[database]}, M = {M}",
               xlabel="$d_{min}$ (degrees)")
        ax.label_outer()
    row[0].set_ylabel("Median LSD (dB)")
axes[0, 0].legend()
save(fig, "error_vs_dmin")

# 4. covering radius against median and 95th-percentile error (per-subject values averaged over subjects)
per_condition = targets.groupby(["database", "subject", "layout", "M", "method"]).agg(
    r_max=("d_min", "max"), median=("200-16k", "median"), p95=("200-16k", lambda x: x.quantile(0.95)))
per_condition = per_condition.groupby(["database", "layout", "M", "method"]).mean().reset_index()
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.6))
for ax, letter, column, label in zip(axes, "ab", ["median", "p95"], ["Median", "95th-percentile"]):
    for r in per_condition[per_condition.method.isin(["NN", "SH-aligned-ridge"])].itertuples():
        ax.scatter(r.r_max, getattr(r, column), marker=MARKERS[r.layout], s=12, linewidths=0.7,
                   facecolor=METHODS[r.method] if r.database == "hutubs" else "none", edgecolor=METHODS[r.method])
    ax.set(xlabel="Covering radius r$_{max}$ (degrees)", ylabel=f"{label} LSD (dB)", title=f"({letter}) {label} error")
handles = [plt.Line2D([], [], color=METHODS[m], marker="o", linestyle="none", label=m) for m in ["NN", "SH-aligned-ridge"]]
handles += [plt.Line2D([], [], color="gray", marker=m, linestyle="none", label=l) for l, m in MARKERS.items()]
axes[1].legend(handles=handles, loc="upper left", ncol=2)
save(fig, "error_vs_rmax")

# 5. spectral error per band and binaural cue errors, all layouts and budgets
fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.2), gridspec_kw={"width_ratios": [3, 1, 1]})
for i, (row, database) in enumerate(zip(axes, DATABASES)):
    df = targets[targets.database == database]
    grouped_bars(row[0], df.groupby("method")[list(BANDS)].median().loc[MAIN], [b.replace("-", "–") for b in BANDS])
    grouped_bars(row[1], df.groupby("method")[["itd_error_us"]].median().loc[MAIN], ["ITD"])
    grouped_bars(row[2], df.groupby("method")[["ild_error_db"]].median().loc[MAIN], ["ILD"])
    row[0].set(ylabel="Median LSD (dB)", xlabel="Band (Hz)", title=f"({'ad'[i]}) {DATABASES[database]}: LSD by band")
    row[1].set(ylabel="Median |ITD error| (µs)", title=f"({'be'[i]}) ITD")
    row[2].set(ylabel="Median |ILD error| (dB)", title=f"({'cf'[i]}) ILD")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center", ncol=4)
save(fig, "band_and_cue_errors")

# 6. directional error maps, gap layout, M = 50
fig, axes = plt.subplots(2, 3, figsize=(6.5, 3.6), sharex=True, sharey=True)
panels = [("d_min", "NN", "$d_{min}$ (degrees)", "viridis", 60),
          ("200-16k", "NN", "NN LSD (dB)", "magma", 20),
          ("200-16k", "SH-aligned-ridge", "ridge SH LSD (dB)", "magma", 20)]
for i, (row, (database, subject)) in enumerate(zip(axes, [("hutubs", "7"), ("sadie", "H3")])):
    df = targets[(targets.database == database) & (targets.subject == subject) &
                 (targets.layout == "gap") & (targets.M == 50)]
    name = f"pp{subject}" if database == "hutubs" else subject
    for j, (ax, (column, method, title, cmap, vmax)) in enumerate(zip(row, panels)):
        d = df[df.method == method]
        sc = ax.scatter((d.az + 180) % 360 - 180, d.el, c=d[column], cmap=cmap, vmin=0, vmax=vmax, s=2)
        ax.set(title=f"({'abcdef'[3 * i + j]}) {name}: {title}", xlim=(-185, 185), ylim=(-70, 95))
        ax.grid(False)
        fig.colorbar(sc, ax=ax)
for ax in axes[1]:
    ax.set_xlabel("Azimuth (degrees)")
for ax in axes[:, 0]:
    ax.set_ylabel("Elevation (degrees)")
save(fig, "error_maps")

# 7. SUpDEq / MCA reference on the subset, scored on the same targets
subset = targets.merge(mca[["database", "subject", "layout", "M"]].drop_duplicates())
both = pd.concat([subset[subset.method.isin(["NN", "SH-aligned", "SH-aligned-ridge"])], mca])
order = ["NN", "SH-aligned", "SUpDEq-OBTA", "SUpDEq-MCA", "SH-aligned-ridge"]
fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.6))
for i, (row, database) in enumerate(zip(axes, DATABASES)):
    df = both[both.database == database]
    grouped_bars(row[0], df.groupby("method")[list(BANDS)].median().loc[order], [b.replace("-", "–") for b in BANDS])
    row[0].set(ylabel="Median LSD (dB)", xlabel="Band (Hz)", title=f"({'ac'[i]}) {DATABASES[database]}: LSD by band")
    for method in order:
        dmin_curve(row[1], df[df.method == method], method)
    row[1].set(xlabel="$d_{min}$ (degrees)", ylabel="Median LSD (dB)",
               title=f"({'bd'[i]}) {DATABASES[database]}: LSD against $d_{{min}}$")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center", ncol=5)
save(fig, "mca_reference")
