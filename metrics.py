import numpy as np

BANDS = {"200-16k": (200, 16000), "200-1k": (200, 1000),
         "1-4k": (1000, 4000), "4-8k": (4000, 8000), "8-16k": (8000, 16000)}


def to_hrtf(ir, fs):
    return np.fft.rfft(ir, axis=-1), np.fft.rfftfreq(ir.shape[-1], 1 / fs)


def lsd(H_ref, H_est, f, band):
    bins = (f >= band[0]) & (f <= band[1])
    ratio = (np.abs(H_ref[..., bins]) + 1e-12) / (np.abs(H_est[..., bins]) + 1e-12)
    return np.sqrt(np.mean((20 * np.log10(ratio)) ** 2, axis=-1)).mean(axis=1)


def itd_us(ir, fs, up=8):
    # lag of the interaural cross-correlation (right relative to left), 8x interpolated.
    # Kept separate from the onset rule used for alignment, so alignment is not scored by its own estimator.
    n = ir.shape[-1]
    X = np.fft.rfft(ir, n=2 * n, axis=-1)
    cross = X[:, 1] * np.conj(X[:, 0])
    cross[:, -1] /= 2  # the old Nyquist bin becomes an interior bin after zero-padding
    c = np.fft.irfft(cross, n=2 * n * up, axis=-1)
    lag = np.argmax(c, axis=-1)
    lag = np.where(lag > n * up, lag - 2 * n * up, lag)
    return lag / (fs * up) * 1e6


def ild_db(ir):
    energy = np.sum(ir ** 2, axis=-1)
    return 10 * np.log10(energy[:, 0] / energy[:, 1])
