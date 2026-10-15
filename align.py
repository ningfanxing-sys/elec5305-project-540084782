import numpy as np
from scipy.signal import butter, resample_poly, sosfilt


def toa(ir, fs, up=8):
    # Onset of the direct sound in seconds. The 1 kHz high-pass removes the slow low-frequency
    # rise that SADIE II HRIRs show before the direct sound; then first sample above 20% of peak.
    x = sosfilt(butter(4, 1000, "highpass", fs=fs, output="sos"), ir, axis=-1)
    x = np.abs(resample_poly(x, up, 1, axis=-1))
    return np.argmax(x > 0.2 * x.max(axis=-1, keepdims=True), axis=-1) / (fs * up)


def advance(H, f, tau):
    # remove a delay tau (seconds) from H; a negative tau adds the delay back
    return H * np.exp(2j * np.pi * f * tau[..., None])
