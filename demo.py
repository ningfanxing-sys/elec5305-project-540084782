"""Listening demo and horizontal-plane figure for one held-out HUTUBS subject and one SADIE II subject.

Every measured direction is reconstructed from the uniform layout. Pink-noise bursts then step around the
listener on the horizontal plane in 10-degree steps, rendered with the measured HRIRs and with each
reconstruction. Run `python mca.py export` and `matlab -batch run_supdeq` first (for MCA).
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, wavfile

from figstyle import save
from mca import DEMO_BUDGETS, DEMO_SUBJECTS, FOLDER
from reconstruct import prepare, setup, sh_estimate
from sh import n_fit
from sofa import load_subject

AUDIO = Path("demo/audio")
STEP, BURST, GAP, RAMP = 10, 0.15, 0.05, 0.01  # degrees, seconds
FROZEN = json.loads(Path("results/frozen_config.json").read_text())["aligned_ridge_sh"]["uniform"]
NAMES = {"hutubs": "HUTUBS pp{}", "sadie": "SADIE II {}"}


def reconstructions(database, subject, M):
    d = prepare(load_subject(database, subject))
    s = setup(d, "uniform", M, targets=np.arange(len(d["az"])))
    mca = loadmat(FOLDER / "out" / f"demo_{database}_{subject}_uniform_{M}.mat")
    return d, {"Measured": d["H"],
               "NN": d["H"][s["nearest"]],
               "SH": sh_estimate(d, s, n_fit(M))[0],
               "SH-aligned-ridge": sh_estimate(d, s, FROZEN[str(M)]["N"], FROZEN[str(M)]["lambda"], aligned=True)[0],
               "SUpDEq-MCA": np.stack([mca["HRTF_L_on"], mca["HRTF_R_on"]], axis=1)}


def horizontal_plane(d):
    plane = np.flatnonzero(np.abs(d["el"]) < 0.5)
    return plane[np.argsort(d["az"][plane])]


def bursts(fs, count):
    rng = np.random.default_rng(0)
    n, ramp = int(BURST * fs), int(RAMP * fs)
    window = np.ones(n)
    window[:ramp] = window[-ramp:][::-1] = 0.5 - 0.5 * np.cos(np.pi * np.arange(ramp) / ramp)
    out = []
    for _ in range(count):
        spectrum = np.fft.rfft(rng.standard_normal(n))
        x = np.fft.irfft(spectrum / np.sqrt(np.maximum(np.arange(len(spectrum)), 1)), n)  # pink
        out.append(window * x / np.abs(x).max())
    return out


def render(d, H, plane):
    fs = d["fs"]
    steps = np.arange(0, 360, STEP)  # counter-clockwise from the front: 90 degrees is the left ear
    hop = int((BURST + GAP) * fs)
    y = np.zeros((hop * len(steps) + 512, 2))
    for k, (target, burst) in enumerate(zip(steps, bursts(fs, len(steps)))):
        i = plane[np.argmin(np.abs((d["az"][plane] - target + 180) % 360 - 180))]
        hrir = np.fft.irfft(H[i], n=d["ir"].shape[-1], axis=-1)
        for ear in range(2):
            out = np.convolve(burst, hrir[ear])
            y[k * hop:k * hop + len(out), ear] += out
    return y


AUDIO.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(2, 5, figsize=(6.5, 3.0), sharex=True, sharey=True)
for i, (row, (database, subject)) in enumerate(zip(axes, DEMO_SUBJECTS)):
    name = NAMES[database].format(subject)
    for M in DEMO_BUDGETS:
        d, estimates = reconstructions(database, subject, M)
        plane = horizontal_plane(d)
        signals = {method: render(d, H, plane) for method, H in estimates.items()}
        gain = 0.9 / np.abs(signals["Measured"]).max()  # one gain per subject keeps level differences
        for method, y in signals.items():
            tag = f"{database}_{subject}_M{M}_{method}".replace("hutubs_", "hutubs_pp")
            wavfile.write(AUDIO / f"{tag}.wav", int(d["fs"]), np.int16(np.clip(gain * y, -1, 1) * 32767))
        if M != 50:
            continue
        f = d["f"]
        band = (f >= 200) & (f <= 16000)
        levels = {method: 20 * np.log10(np.abs(H[plane][:, 0][:, band]) + 1e-12) for method, H in estimates.items()}
        top = np.ceil(np.percentile(levels["Measured"], 99) / 5) * 5  # colour range set by the measured HRTFs
        for j, (ax, (method, level)) in enumerate(zip(row, levels.items())):
            mesh = ax.pcolormesh(d["az"][plane], f[band] / 1e3, level.T, shading="nearest", cmap="viridis",
                                 vmin=top - 45, vmax=top, linewidth=0, antialiased=False)
            ax.set(yscale="log", xticks=[0, 90, 180, 270], yticks=[0.5, 1, 2, 4, 8, 16],
                   yticklabels=["0.5", "1", "2", "4", "8", "16"])
            ax.set_title(f"({'abcdefghij'[5 * i + j]}) {method.replace('SUpDEq-', '')}", fontsize=7.5)
            ax.grid(False)
        row[0].set_ylabel(f"{name}\nFrequency (kHz)")
        colorbar = fig.colorbar(mesh, ax=row, shrink=0.9, pad=0.02, label="Left ear (dB)")
        colorbar.solids.set_rasterized(False)
for ax in axes[1]:
    ax.set_xlabel("Azimuth (°)")
save(fig, "hrtf_horizontal_plane")
print("wrote", len(list(AUDIO.glob("*.wav"))), "WAV files to", AUDIO)
