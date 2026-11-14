import numpy as np

from align import advance, toa
from metrics import BANDS, ild_db, itd_us, lsd, lsd_ears, to_hrtf
from sh import n_max, sh_fit, sh_matrix, singular_values
from sphere import angles, map_to_grid, nominal_layouts, to_xyz, voronoi_weights

BUDGETS = [25, 50, 100]
LAYOUTS = ["uniform", "random1", "random2", "gap"]
TOA_ORDER = 3  # SH order used to interpolate the times of arrival


def prepare(data):
    H, f = to_hrtf(data["ir"], data["fs"])
    return dict(data, H=H, f=f, tau=toa(data["ir"], data["fs"]), xyz=to_xyz(data["az"], data["el"]))


def setup(d, layout, M, targets=None):
    """Context and target directions plus everything the SH fits share.

    By default every measured direction outside the layout is a target."""
    context = map_to_grid(to_xyz(*nominal_layouts(M)[layout]), d["xyz"])[0]
    if targets is None:
        targets = np.setdiff1d(np.arange(len(d["xyz"])), context)
    distances = angles(d["xyz"][targets], d["xyz"][context])
    Y = sh_matrix(n_max(M), d["az"], d["el"])
    w = voronoi_weights(d["xyz"][context])
    K = (TOA_ORDER + 1) ** 2
    tau_coef = sh_fit(Y[context, :K], d["tau"][context], w)
    return dict(context=context, targets=targets, d_min=distances.min(axis=1),
                nearest=context[np.argmin(distances, axis=1)], Y=Y, w=w,
                tau_targets=Y[targets, :K] @ tau_coef, tau_context_fit=Y[context, :K] @ tau_coef)


def sh_estimate(d, s, N, lam=0.0, aligned=False, smooth=False):
    """SH reconstruction at the targets; returns the estimate, kappa and coefficient norm.

    smooth=True removes the order-3 SH fit of the context onsets instead of the onsets themselves
    (follow-up experiment on onset roughness)."""
    K = (N + 1) ** 2
    Y_c, Y_t = s["Y"][s["context"], :K], s["Y"][s["targets"], :K]
    H_c = d["H"][s["context"]]
    if aligned:
        H_c = advance(H_c, d["f"], s["tau_context_fit"] if smooth else d["tau"][s["context"]])
    a = sh_fit(Y_c, H_c, s["w"], lam)
    estimate = np.tensordot(Y_t, a, axes=1)
    if aligned:
        estimate = advance(estimate, d["f"], -s["tau_targets"])
    sv = singular_values(Y_c, s["w"])
    return estimate, sv[0] / sv[-1], np.linalg.norm(a)


def lsd_bands(d, s, estimate):
    return {band: lsd(d["H"][s["targets"]], estimate, d["f"], limits) for band, limits in BANDS.items()}


def ear_lsd(d, s, estimate):
    # broadband LSD of the left and right ear, for the ipsilateral/contralateral analysis
    ears = lsd_ears(d["H"][s["targets"]], estimate, d["f"], BANDS["200-16k"])
    return dict(lsd_left=ears[:, 0], lsd_right=ears[:, 1])


def cue_errors(d, s, estimate):
    ir_true = d["ir"][s["targets"]]
    ir_est = np.fft.irfft(estimate, n=ir_true.shape[-1], axis=-1)
    return (np.abs(itd_us(ir_est, d["fs"]) - itd_us(ir_true, d["fs"])),
            np.abs(ild_db(ir_est) - ild_db(ir_true)))
