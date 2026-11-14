"""Result figures of the report (vector PDF in figures/pdf, PNG for the web page)."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from align import toa
from figstyle import save
from metrics import BANDS
from sofa import load_sofa, load_subject

BUDGETS = [25, 50, 100]
LAYOUTS = ["uniform", "random1", "random2", "gap"]
MARKERS = dict(zip(LAYOUTS, ["o", "s", "^", "D"]))
METHODS = {"NN": ("#0072B2", "o"), "SH": ("#E69F00", "s"), "SH-aligned": ("#009E73", "^"),
           "SH-aligned-ridge": ("#CC79A7", "D"), "SUpDEq-OBTA": ("#56B4E9", "v"),
           "SUpDEq-MCA": ("#D55E00", "P"), "SH-aligned-smooth": ("#882255", "X")}
MAIN = ["NN", "SH", "SH-aligned", "SH-aligned-ridge"]
M_COLORS = dict(zip(BUDGETS, ["#0072B2", "#E69F00", "#009E73"]))
DATABASES = {"hutubs": "HUTUBS held-out", "sadie": "SADIE II"}
BAND_LABELS = ["0.2–16", "0.2–1", "1–4", "4–8", "8–16"]

selection = pd.read_csv("results/selection.csv")
summary = pd.read_csv("results/eval_summary.csv").rename(columns={"median_200-16k": "lsd"})
targets = pd.read_csv("results/eval_targets.csv.gz", dtype={"subject": str})
mca = pd.read_csv("results/mca_targets.csv.gz", dtype={"subject": str})
roughness = pd.read_csv("results/toa_roughness.csv")
layouts = json.loads(open("results/layouts.json").read())


def color(method):
    return METHODS[method][0]


def binned_median(ax, x, y, bins, method, linestyle="-", filled=True, min_count=50, label=None):
    groups = pd.Series(y.values).groupby(pd.cut(x.values, bins), observed=True)
    medians = groups.median()[groups.size() >= min_count]
    c, m = METHODS[method]
    ax.plot([b.mid for b in medians.index], medians.values, linestyle=linestyle, marker=m, markersize=3.5,
            color=c, markerfacecolor=c if filled else "white", label=label)


def grouped_bars(ax, values, labels):
    """values: DataFrame with methods as rows and groups as columns."""
    x = np.arange(len(labels))
    width = 0.8 / len(values)
    for i, (method, row) in enumerate(values.iterrows()):
        ax.bar(x + (i - (len(values) - 1) / 2) * width, row.values, width, color=color(method), label=method)
    ax.set_xticks(x, labels)
    ax.grid(axis="x", visible=False)


def vector_colorbar(fig, mappable, ax, label):
    bar = fig.colorbar(mappable, ax=ax, label=label)
    bar.solids.set_rasterized(False)


def condition_legend(fig, extra=()):
    handles = [plt.Line2D([], [], color=c, marker="o", linestyle="none", label=f"M = {M}") for M, c in M_COLORS.items()]
    handles += [plt.Line2D([], [], color="gray", marker=m, linestyle="none", label=l) for l, m in MARKERS.items()]
    fig.legend(handles=handles + list(extra), loc="outside lower center", ncol=4)


# 1. directional maps for the gap layout (M = 50): error follows the distance to the nearest context direction
fig, axes = plt.subplots(2, 3, figsize=(6.5, 3.05), sharex=True, sharey=True)
panels = [("d_min", "NN", "$d_{min}$ (°)", "viridis", 60),
          ("200-16k", "NN", "NN LSD (dB)", "magma", 20),
          ("200-16k", "SH-aligned-ridge", "SH-aligned-ridge LSD (dB)", "magma", 20)]
for i, (row, (database, subject, grid, name)) in enumerate(zip(axes, [("hutubs", "7", "hutubs", "HUTUBS pp7"),
                                                                    ("sadie", "H3", "sadie_h3", "SADIE II H3")])):
    df = targets[(targets.database == database) & (targets.subject == subject) &
                 (targets.layout == "gap") & (targets.M == 50)]
    data = load_subject(database, int(subject) if database == "hutubs" else subject)
    context = layouts[grid]["gap"]["50"]
    for j, (ax, (column, method, title, cmap, vmax)) in enumerate(zip(row, panels)):
        d = df[df.method == method]
        sc = ax.scatter((d.az + 180) % 360 - 180, d.el, c=d[column], cmap=cmap, vmin=0, vmax=vmax,
                        s=9 if database == "hutubs" else 2.5, linewidths=0)
        ax.scatter((data["az"][context] + 180) % 360 - 180, data["el"][context], marker="x", s=10,
                   color="white", linewidths=1.6)
        ax.scatter((data["az"][context] + 180) % 360 - 180, data["el"][context], marker="x", s=8,
                   color="black", linewidths=0.8)
        ax.set(title=f"({'abcdef'[3 * i + j]}) {title}", xlim=(-185, 185), ylim=(-70, 95),
               xticks=[-180, -90, 0, 90, 180])
        ax.grid(False)
        vector_colorbar(fig, sc, ax, "")
    row[0].set_ylabel(f"{name}\nElevation (°)")
for ax in axes[1]:
    ax.set_xlabel("Azimuth (°, left = +90)")
save(fig, "error_maps")

# 2. error against distance to the nearest context direction: one panel per method, both databases
fig, axes = plt.subplots(1, 4, figsize=(6.5, 2.3), sharey=True)
bins = np.arange(0, 65, 5)
for ax, letter, method in zip(axes, "abcd", MAIN):
    for database, style, filled in [("hutubs", "-", True), ("sadie", "--", False)]:
        df = targets[(targets.database == database) & (targets.method == method)]
        binned_median(ax, df.d_min, df["200-16k"], bins, method, style, filled, label=DATABASES[database])
    ax.set(title=f"({letter}) {method}", xlabel="$d_{min}$ (°)", xticks=[0, 20, 40, 60], ylim=(0, 30))
axes[0].set_ylabel("Median LSD (dB)")
axes[0].legend(loc="upper left")
save(fig, "error_vs_dmin")

# 3. conditioning: error against kappa, energy amplification and the benefit of ridge regularisation
fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.05))
fits = []
for ax, letter, method in zip(axes[0], "ab", ["SH", "SH-aligned"]):
    for database, fill, style in [("hutubs", True, "-"), ("sadie", False, "--")]:
        g = summary[(summary.method == method) & (summary.database == database)]
        g = g.groupby(["layout", "M", "N"])[["kappa", "lsd"]].mean().reset_index()
        for r in g.itertuples():
            ax.scatter(np.log10(r.kappa), r.lsd, marker=MARKERS[r.layout], s=12, linewidths=0.7,
                       facecolor=M_COLORS[r.M] if fill else "none", edgecolor=M_COLORS[r.M])
        high = g[g.kappa > 10]
        slope, intercept = np.polyfit(np.log10(high.kappa), high.lsd, 1)
        x = np.linspace(1, np.log10(high.kappa.max()), 10)
        ax.plot(x, slope * x + intercept, style, color="black", linewidth=1,
                label=f"{DATABASES[database]}: {slope:.1f} dB/decade")
    ax.set(title=f"({letter}) {method}, λ = 0", xlabel="log$_{10}$ κ", ylabel="Median LSD (dB)")
    ax.legend(loc="upper left")
ax = axes[1, 0]
for method, edge, marker in [("SH", "black", "v"), ("SH-aligned-ridge", "#CC79A7", "P")]:
    for database, fill in [("hutubs", True), ("sadie", False)]:
        g = summary[(summary.method == method) & (summary.database == database)]
        g = g.groupby(["layout", "M", "N"])[["kappa", "energy_db"]].median().reset_index()
        ax.scatter(np.log10(g.kappa), g.energy_db, s=12, marker=marker, linewidths=0.7,
                   facecolor=edge if fill else "none", edgecolor=edge, label=f"{method}, {DATABASES[database]}")
ax.axhline(0, color="0.4", linewidth=0.8)
ax.set(title="(c) Energy of the reconstruction", xlabel="log$_{10}$ κ",
       ylabel="Energy re measured (dB)", ylim=(-15, 135))
ax.legend(loc="upper left", fontsize=6.8, handletextpad=0.3, borderaxespad=0.2)
ax = axes[1, 1]
mean = selection.groupby(["layout", "M", "N", "lam"])[["kappa", "median_lsd"]].mean().reset_index()
for (layout, M, N), g in mean.groupby(["layout", "M", "N"]):
    plain = g[g.lam == 0].iloc[0]
    ax.scatter(np.log10(plain.kappa), plain.median_lsd - g.median_lsd.min(), color=M_COLORS[M],
               marker=MARKERS[layout], s=12)
ax.set(title="(d) Benefit of ridge (HUTUBS selection half)", xlabel="log$_{10}$ κ at λ = 0",
       ylabel="LSD reduction at best λ (dB)")
condition_legend(fig)
save(fig, "conditioning")

# 4. covering radius against median and 95th-percentile error (per-subject values averaged over subjects)
per_condition = targets.groupby(["database", "subject", "layout", "M", "method"]).agg(
    r_max=("d_min", "max"), median=("200-16k", "median"), p95=("200-16k", lambda x: x.quantile(0.95)))
per_condition = per_condition.groupby(["database", "layout", "M", "method"]).mean().reset_index()
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.15))
for ax, letter, column, label in zip(axes, "ab", ["median", "p95"], ["Median", "95th-percentile"]):
    for r in per_condition[per_condition.method.isin(["NN", "SH-aligned-ridge"])].itertuples():
        ax.scatter(r.r_max, getattr(r, column), marker=MARKERS[r.layout], s=13, linewidths=0.7,
                   facecolor=color(r.method) if r.database == "hutubs" else "none", edgecolor=color(r.method))
    ax.set(xlabel="Covering radius $r_{max}$ (°)", ylabel=f"{label} LSD (dB)", title=f"({letter}) {label} error")
handles = [plt.Line2D([], [], color=color(m), marker="o", linestyle="none", label=m) for m in ["NN", "SH-aligned-ridge"]]
handles += [plt.Line2D([], [], color="gray", marker=m, linestyle="none", label=l) for l, m in MARKERS.items()]
fig.legend(handles=handles, loc="outside right center", ncol=1)
save(fig, "error_vs_rmax")

# 5. time alignment: onset detection, onset roughness and what alignment buys in each database
fig, axes = plt.subplots(1, 4, figsize=(6.5, 2.55), gridspec_kw={"width_ratios": [1.25, 0.8, 1, 1]})
ax = axes[0]
data = load_sofa("data/raw/sadie/H3_44K_16bit_256tap_FIR_SOFA.sofa")
index = int(np.argmin(np.abs(data["az"] - 90) + np.abs(data["el"])))
h, fs = data["ir"][index], data["fs"]
t = np.arange(135) / fs * 1e3
onset_hp = toa(h, fs) * 1e3
for ear, (label, c) in enumerate([("near ear", "#0072B2"), ("far ear", "#D55E00")]):
    x = h[ear] / np.abs(h[ear]).max()
    ax.plot(t, x[:135], color=c, linewidth=0.9, label=label)
    raw = np.argmax(np.abs(x) > 0.1) / fs * 1e3
    ax.plot(raw, 0.1 * np.sign(x[np.argmax(np.abs(x) > 0.1)]), "x", color=c, markersize=5, mew=1.3)
    ax.axvline(onset_hp[ear], color=c, linestyle=":", linewidth=1.1)
ax.axhline(0.1, color="0.5", linestyle="--", linewidth=0.7)
ax.set(title="(a) SADIE II H3, source left", xlabel="Time (ms)", ylabel="Amplitude (peak = 1)")
ax.legend(loc="upper left", fontsize=7, handlelength=1.2)
ax = axes[1]
for k, (database, c) in enumerate([("hutubs", "#0072B2"), ("sadie", "#D55E00")]):
    r = roughness[roughness.database == database]
    jitter = np.random.default_rng(k).uniform(-0.15, 0.15, len(r))
    ax.scatter(k + jitter, r.median_us, s=7, color=c)
ax.set(xticks=[0, 1], xticklabels=["HUTUBS", "SADIE II"], xlim=(-0.5, 1.5), ylim=(0, 35),
       ylabel="Median |TOA − SH fit| (µs)", title="(b) Onset roughness")
ax.grid(axis="x", visible=False)
order = ["NN", "SH", "SH-aligned", "SH-aligned-smooth"]
for ax, letter, column, label, title in [(axes[2], "c", "8-16k", "Median LSD, 8–16 kHz (dB)", "Spectrum above 8 kHz"),
                                         (axes[3], "d", "itd_error_us", "Median |ITD error| (µs)", "ITD")]:
    values = targets[targets.method.isin(order)].groupby(["method", "database"])[column].median().unstack()
    grouped_bars(ax, values.loc[order], ["HUTUBS", "SADIE II"])
    ax.set(ylabel=label, title=f"({letter}) {title}")
fig.legend(*axes[2].get_legend_handles_labels(), loc="outside lower center", ncol=4)
save(fig, "alignment")

# 6. magnitude correction on the SUpDEq subset: frequency bands, coverage and ear side, on identical targets
subset = targets.merge(mca[["database", "subject", "layout", "M"]].drop_duplicates())
both = pd.concat([subset[subset.method.isin(["NN", "SH-aligned", "SH-aligned-ridge"])], mca])
order = ["NN", "SH-aligned", "SUpDEq-OBTA", "SUpDEq-MCA", "SH-aligned-ridge"]
lateral = np.rad2deg(np.arcsin(np.cos(np.deg2rad(both.el)) * np.sin(np.deg2rad(both.az))))  # + = left
ears = pd.concat([both.assign(rel=lateral, ear_lsd=both.lsd_left), both.assign(rel=-lateral, ear_lsd=both.lsd_right)])
fig, axes = plt.subplots(2, 3, figsize=(6.5, 3.95), gridspec_kw={"width_ratios": [1.25, 1, 1]})
for i, (row, database) in enumerate(zip(axes, DATABASES)):
    df = both[both.database == database]
    grouped_bars(row[0], df.groupby("method")[list(BANDS)].median().loc[order], BAND_LABELS)
    row[0].set(ylabel="Median LSD (dB)", xlabel="Band (kHz)", title=f"({'ad'[i]}) {DATABASES[database]}: bands")
    for method in order:
        d = df[df.method == method]
        binned_median(row[1], d.d_min, d["200-16k"], np.arange(0, 50, 5), method, min_count=30)
        e = ears[(ears.database == database) & (ears.method == method)]
        binned_median(row[2], e.rel, e.ear_lsd, np.arange(-90, 91, 15), method, min_count=30)
    row[1].set(xlabel="$d_{min}$ (°)", ylabel="Median LSD (dB)", title=f"({'be'[i]}) Coverage")
    row[2].set(xlabel="Lateral angle re ear (°)\ncontralateral ← → ipsilateral", ylabel="Median LSD of that ear (dB)",
               title=f"({'cf'[i]}) Ear side", xticks=[-90, -45, 0, 45, 90])
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside lower center", ncol=5)
save(fig, "mca_reference")

# Appendix: preferred ridge parameter against kappa (selection half)
fig, ax = plt.subplots(figsize=(2.3, 2.5))
for (layout, M, N), g in mean.groupby(["layout", "M", "N"]):
    best = g.loc[g.median_lsd.idxmin()]
    ax.scatter(np.log10(g[g.lam == 0].kappa.iloc[0]), max(best.lam, 1e-6), color=M_COLORS[M],
               marker=MARKERS[layout], s=12)
ax.set(yscale="log", xlabel="log$_{10}$ κ at λ = 0", ylabel="Preferred λ (0 plotted at 10$^{-6}$)",
       title="(c) Preferred ridge parameter")
save(fig, "lambda_vs_kappa")
