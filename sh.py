import numpy as np

try:
    from scipy.special import sph_harm_y as _sph_harm
    _new_signature = True
except ImportError:
    from scipy.special import sph_harm as _sph_harm
    _new_signature = False


def n_max(M):
    # highest order with (N+1)^2 <= M
    return min(10, int(np.sqrt(M)) - 1)


def n_fit(M):
    # fixed rule chosen before looking at errors: about two context directions per coefficient
    return int(np.sqrt(M / 2)) - 1


def sh_matrix(N, az, el):
    az = np.deg2rad(az)
    polar = np.deg2rad(90 - np.asarray(el))
    Y = np.empty((len(az), (N + 1) ** 2))
    for n in range(N + 1):
        for m in range(-n, n + 1):
            value = (_sph_harm(n, abs(m), polar, az) if _new_signature
                     else _sph_harm(abs(m), n, az, polar))
            real = value.imag if m < 0 else value.real
            Y[:, n * n + n + m] = real * (np.sqrt(2) * (-1) ** m if m else 1)
    return Y


def sh_fit(Y, h, w, lam=0.0):
    root_w = np.sqrt(w)[:, None]
    A = root_w * Y
    b = root_w * h.reshape(len(h), -1)
    if lam:
        A = np.vstack((A, np.sqrt(lam) * np.eye(Y.shape[1])))
        b = np.vstack((b, np.zeros((Y.shape[1], b.shape[1]), dtype=h.dtype)))
    coef = np.linalg.lstsq(A, b, rcond=None)[0]
    return coef.reshape((Y.shape[1],) + h.shape[1:])


def singular_values(Y, w):
    return np.linalg.svd(np.sqrt(w)[:, None] * Y, compute_uv=False)
