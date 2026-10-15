import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from align import toa
from figstyle import save
from metrics import ild_db, itd_us
from sh import sh_fit, sh_matrix
from sofa import HUTUBS_TEST, SADIE_SUBJECTS, load_sofa, load_subject
from sphere import EL_MIN, angles, to_xyz, voronoi_weights

FILES = {
    "HUTUBS pp2": "data/raw/hutubs/pp2_HRIRs_measured.sofa",
    "SADIE H3": "data/raw/sadie/H3_44K_16bit_256tap_FIR_SOFA.sofa",
    "SADIE H10": "data/raw/sadie/H10_44K_16bit_256tap_FIR_SOFA.sofa",
}
RESULTS = Path("results")
FIGURES = Path("figures")
LANDMARKS = {"Front": (0, 0), "Left": (90, 0), "Back": (-180, 0),
             "Right": (-90, 0), "Top": (0, 90)}

RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.4), sharey=True)
ir_fig, ir_axes = plt.subplots(1, 3, figsize=(6.5, 2.1), sharey=True)
itd_fig, itd_axes = plt.subplots(1, 3, figsize=(6.5, 2.1), sharey=True)
rows = []
print(f"{'Dataset':<12} {'Source':<6} {'10% onset':>10} {'HP onset':>9} {'xcorr':>7} {'ILD dB':>7} {'Sign':>5}")
for letter, ax, ir_ax, itd_ax, (name, path) in zip("abc", axes, ir_axes, itd_axes, FILES.items()):
    data = load_sofa(path)
    az, el, ir, fs = data["az"], data["el"], data["ir"], data["fs"]
    xyz = to_xyz(az, el)
    ax.scatter((az + 180) % 360 - 180, el, s=1.5, alpha=0.6, color="tab:blue")
    for label, point in LANDMARKS.items():
        offset = (0, 6) if label == "Top" else (3, 3)
        ax.annotate(label, point, xytext=offset, textcoords="offset points",
                    ha="center" if label == "Top" else "left", fontsize=6,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
    ax.set(title=f"({letter}) {name}: {len(az)} directions", xlabel="Azimuth (degrees)",
           xlim=(-188, 188), ylim=(-98, 108),
           xticks=[-180, -90, 0, 90, 180], yticks=[-90, -45, 0, 45, 90])
    ax.grid(alpha=0.2)

    # ITD from the alignment onset rule over the whole sphere, against lateral angle (+ = left)
    lateral = np.rad2deg(np.arcsin(xyz[:, 1]))
    tau = toa(ir, fs)
    itd_ax.scatter(lateral, (tau[:, 1] - tau[:, 0]) * 1e6, s=1, alpha=0.4)
    itd_ax.set(title=f"({'def'['abc'.index(letter)]}) {name}", xlabel="Lateral angle (degrees)",
               xticks=[-90, -45, 0, 45, 90])
    itd_ax.grid(alpha=0.2)

    for side, source_az, sign in [("left", 90, 1), ("right", 270, -1)]:
        distances = angles(to_xyz([source_az], [0]), xyz)[0]
        index = int(np.argmin(distances))
        h = ir[index:index + 1]
        onset = np.argmax(np.abs(h[0]) > 0.1 * np.abs(h[0]).max(axis=1)[:, None], axis=1)
        itd_onset = float((onset[1] - onset[0]) / fs * 1e6)  # rule from the proposal
        itd_hp = float((tau[index, 1] - tau[index, 0]) * 1e6)  # rule used for alignment
        itd = float(itd_us(h, fs)[0])  # cross-correlation, used for scoring
        ild = float(ild_db(h)[0])
        ok = bool(sign * itd_hp > 0 and sign * itd > 0 and sign * ild > 0)
        if side == "left":
            t = np.arange(160) / fs * 1e3
            for ear, label in enumerate(["left ear", "right ear"]):
                ir_ax.plot(t, h[0, ear, :160] / np.abs(h[0, ear]).max(), label=label)
            ir_ax.axhline(0.1, color="0.5", linestyle="--", linewidth=0.8, label="10% onset threshold")
            ir_ax.set(title=f"({letter}) {name}", xlabel="Time (ms)")
            if letter == "a":
                ir_ax.legend()
        rows.append({"dataset": name, "source": side, "az": float(az[index]), "el": float(el[index]),
                     "itd_onset10_us": itd_onset, "itd_hp_onset_us": itd_hp, "itd_xcorr_us": itd,
                     "ild_db": ild, "sign_ok": ok})
        print(f"{name:<12} {side:<6} {itd_onset:10.0f} {itd_hp:9.0f} {itd:7.0f} {ild:7.2f} "
              f"{'PASS' if ok else 'FAIL':>5}")

itd_axes[0].set_ylabel("ITD from onset rule (µs)")
ir_axes[0].set_ylabel("Amplitude (peak = 1)")
axes[0].set_ylabel("Elevation (degrees)")
save(fig, "measured_directions")
save(ir_fig, "ir_left_source")
save(itd_fig, "itd_vs_lateral")
with (RESULTS / "coordinate_check.csv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

# SADIE II: where are the implausible onsets? Decides EL_MIN (impulse responses only, no reconstruction errors).
ring_rows = []
for subject in SADIE_SUBJECTS:
    data = load_sofa(f"data/raw/sadie/{subject}_44K_16bit_256tap_FIR_SOFA.sofa")
    sin_lateral = to_xyz(data["az"], data["el"])[:, 1]
    tau = toa(data["ir"], data["fs"])
    itd = (tau[:, 1] - tau[:, 0]) * 1e6
    A = np.column_stack([sin_lateral, np.ones_like(sin_lateral)])
    outlier = np.abs(itd - A @ np.linalg.lstsq(A, itd, rcond=None)[0]) > 400
    low = data["el"] < EL_MIN
    ring_rows.append(dict(subject=subject, directions=len(itd), outliers=int(outlier.sum()),
                          outliers_below_el_min=int((outlier & low).sum()), directions_below_el_min=int(low.sum())))
ring_table = pd.DataFrame(ring_rows)
ring_table.to_csv(RESULTS / "sadie_onset_outliers.csv", index=False)
print("\nSADIE II directions whose onset ITD is more than 400 us from a sin(lateral) fit")
print(ring_table.to_string(index=False))

# How smooth are the onset estimates across directions? Residual of an order-3 SH fit to each ear's TOA
# over all directions used (elevation >= EL_MIN), per subject.
rough_rows = []
for database, subjects in [("hutubs", HUTUBS_TEST), ("sadie", SADIE_SUBJECTS)]:
    for subject in subjects:
        data = load_subject(database, subject)
        tau = toa(data["ir"], data["fs"]) * 1e6
        Y = sh_matrix(3, data["az"], data["el"])
        residual = np.abs(tau - Y @ sh_fit(Y, tau, voronoi_weights(to_xyz(data["az"], data["el"]))))
        rough_rows.append(dict(database=database, subject=subject, median_us=np.median(residual),
                               p95_us=np.percentile(residual, 95)))
rough = pd.DataFrame(rough_rows)
rough.to_csv(RESULTS / "toa_roughness.csv", index=False)
print("\nTOA residual from an order-3 SH fit (mean over subjects, us)")
print(rough.groupby("database")[["median_us", "p95_us"]].mean().round(1).to_string())

# Onset rule: compare candidate rules with the cross-correlation ITD on the horizontal plane
from scipy.signal import butter, resample_poly, sosfilt  # noqa: E402

def onset_itd(ir, fs, highpass, threshold, up=8):
    x = sosfilt(butter(4, 1000, "highpass", fs=fs, output="sos"), ir, axis=-1) if highpass else ir
    x = np.abs(resample_poly(x, up, 1, axis=-1))
    t = np.argmax(x > threshold * x.max(axis=-1, keepdims=True), axis=-1) / (fs * up)
    return (t[:, 1] - t[:, 0]) * 1e6

rule_rows = []
for name, path in FILES.items():
    data = load_sofa(path)
    plane = np.abs(data["el"]) < 0.1
    ir, fs = data["ir"][plane], data["fs"]
    reference = itd_us(ir, fs)
    for highpass, threshold in [(False, 0.1), (True, 0.1), (True, 0.2)]:
        itd = onset_itd(ir, fs, highpass, threshold)
        lateral = np.abs(reference) > 200
        rule_rows.append(dict(dataset=name, rule=f"{'1 kHz high-pass, ' if highpass else ''}{threshold:.0%} threshold",
                              sign_errors=int(np.sum(np.sign(itd[lateral]) != np.sign(reference[lateral]))),
                              median_abs_diff_us=np.median(np.abs(itd - reference))))
rules = pd.DataFrame(rule_rows)
rules.to_csv(RESULTS / "onset_rule_check.csv", index=False)
print("\nOnset rules against the cross-correlation ITD, horizontal plane")
print(rules.round(1).to_string(index=False))
