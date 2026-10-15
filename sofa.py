import h5py

from sphere import EL_MIN

# Pilot subjects pp2-pp6 plus even IDs select parameters; odd IDs from pp7 are held out.
# pp22 and pp88 are the same person and both land in the selection half.
HUTUBS_SELECT = [2, 3, 4, 5, 6] + list(range(8, 96, 2))
HUTUBS_TEST = list(range(7, 96, 2))
SADIE_SUBJECTS = [f"H{i}" for i in range(3, 21)]


def load_sofa(path):
    with h5py.File(path, "r") as f:
        attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f["SourcePosition"].attrs.items()}
        assert attrs["Type"] == "spherical" and attrs["Units"].startswith("degree")
        assert f["ReceiverPosition"][0, 1, 0] > 0  # receiver 0 is the left ear
        assert not f["Data.Delay"][:].any()  # all files used here store the full delay in Data.IR
        position = f["SourcePosition"][:]
        return {"ir": f["Data.IR"][:], "fs": float(f["Data.SamplingRate"][0]),
                "az": position[:, 0], "el": position[:, 1]}


def load_subject(database, subject):
    if database == "hutubs":
        data = load_sofa(f"data/raw/hutubs/pp{subject}_HRIRs_measured.sofa")
    else:
        data = load_sofa(f"data/raw/sadie/{subject}_44K_16bit_256tap_FIR_SOFA.sofa")
    keep = data["el"] >= EL_MIN
    return {key: value if key == "fs" else value[keep] for key, value in data.items()}
