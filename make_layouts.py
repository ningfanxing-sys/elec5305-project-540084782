import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from figstyle import save
from sh import n_fit, n_max, sh_matrix, singular_values
from sofa import load_subject
from sphere import angles, map_to_grid, nominal_layouts, to_xyz, voronoi_weights

M_VALUES = (25, 50, 100)
LAYOUTS = ("uniform", "random1", "random2", "gap")
GRIDS = {"hutubs": ("hutubs", 2), "sadie_h3": ("sadie", "H3"), "sadie_h10": ("sadie", "H10")}
RESULTS = Path("results")
FIGURES = Path("figures")


def nearest_spacing(xyz):
    distance = angles(xyz, xyz)
    np.fill_diagonal(distance, np.inf)
    return distance.min(axis=1)


def make_figures(data, layouts):
    az, el = data["az"], data["el"]
    wrapped = (az + 180) % 360 - 180
    xyz = to_xyz(az, el)
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 3.8), sharex=True, sharey=True)
    for letter, ax, layout in zip("abcd", axes.flat, LAYOUTS):
        index = layouts[layout]["50"]
        ax.scatter(wrapped, el, s=3, color="0.8", alpha=0.5)
        ax.scatter(wrapped[index], el[index], s=9, color="tab:blue")
        ax.set(title=f"({letter}) {layout}", xlabel="Azimuth (degrees; left = +90)",
               ylabel="Elevation (degrees)", xlim=(-180, 180), ylim=(-90, 90))
        ax.label_outer()
    save(fig, "layouts_M50")

    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.2), sharey=True)
    for letter, ax, M in zip("abc", axes, M_VALUES):
        spacing = [nearest_spacing(xyz[layouts[name][str(M)]])
                   for name in ("uniform", "random1")]
        bins = np.linspace(0, max(values.max() for values in spacing) + 1, 14)
        for name, values in zip(("uniform", "random1"), spacing):
            ax.hist(values, bins=bins, histtype="step", linewidth=1.8, label=name)
        ax.set(title=f"({letter}) M = {M}", xlabel="NN spacing (degrees)")
        if letter == "a":
            ax.set_ylabel("Context directions")
            ax.legend()
    save(fig, "nn_spacing")

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.6))
    for ax, settings in zip(axes, ([(M, n_max(M)) for M in M_VALUES],
                                   [(50, N) for N in range(2, 7)])):
        for M, N in settings:
            index = layouts["uniform"][str(M)]
            w = voronoi_weights(xyz[index])
            sv = singular_values(sh_matrix(N, az[index], el[index]), w)
            ax.semilogy(np.arange(1, len(sv) + 1), sv, label=f"M={M}, N={N}")
        ax.set(xlabel="Singular value index", ylabel="Singular value", ylim=(1e-3, 3))
        ax.legend()
    axes[0].set_title("(a) Uniform layouts at their highest order")
    axes[1].set_title("(b) Uniform layout, M = 50, orders 2 to 6")
    save(fig, "singular_values")


def main():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    nominal = {M: nominal_layouts(M) for M in M_VALUES}
    layouts, rows = {}, []
    hutubs = None
    for grid, (database, subject) in GRIDS.items():
        data = load_subject(database, subject)
        az, el = data["az"], data["el"]
        xyz = to_xyz(az, el)
        layouts[grid] = {name: {} for name in LAYOUTS}
        for name in LAYOUTS:
            for M in M_VALUES:
                index, error = map_to_grid(to_xyz(*nominal[M][name]), xyz)
                layouts[grid][name][str(M)] = index.tolist()
                target = np.ones(len(xyz), dtype=bool)
                target[index] = False
                covering = angles(xyz[target], xyz[index]).min(axis=1).max()
                spacing = nearest_spacing(xyz[index])
                Y, w = sh_matrix(n_max(M), az[index], el[index]), voronoi_weights(xyz[index])
                kappa = {}
                for N in (n_fit(M), n_max(M)):
                    sv = singular_values(Y[:, :(N + 1) ** 2], w)
                    kappa[N] = sv[0] / sv[-1]
                rows.append(dict(grid=grid, layout=name, M=M,
                                 mapping_mean=error.mean(), mapping_median=np.median(error),
                                 mapping_max=error.max(), covering_radius=covering,
                                 nn_median=np.median(spacing), nn_max=spacing.max(),
                                 N_fit=n_fit(M), kappa_fit=kappa[n_fit(M)],
                                 N_max=n_max(M), kappa_max=kappa[n_max(M)]))
        if grid == "hutubs":
            hutubs = data
    (RESULTS / "layouts.json").write_text(json.dumps(layouts, indent=2) + "\n")
    with (RESULTS / "geometry.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print("grid       layout     M  map_mean map_med map_max covering nn_med nn_max  kappa(N_fit) kappa(N_max)")
    for r in rows:
        print(f"{r['grid']:<10} {r['layout']:<8} {r['M']:3} "
              f"{r['mapping_mean']:8.2f} {r['mapping_median']:7.2f} {r['mapping_max']:7.2f} "
              f"{r['covering_radius']:8.2f} {r['nn_median']:6.2f} {r['nn_max']:6.2f} "
              f"{r['kappa_fit']:9.1f} (N={r['N_fit']}) {r['kappa_max']:9.1e} (N={r['N_max']})")
    print("All mapping, covering and spacing values are in degrees.")
    make_figures(hutubs, layouts["hutubs"])


if __name__ == "__main__":
    main()
