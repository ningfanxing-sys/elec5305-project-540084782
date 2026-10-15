"""Download the HRTF files used so far and write data/manifest.csv.

HUTUBS: measured HRIRs, release doi:10.14279/depositonce-8487 (files from the sofacoustics.org mirror).
SADIE II: version 2-2 from Zenodo (record 12092466). Only the 44.1 kHz HRIR
SOFA file is extracted, so both databases share one sample rate.
"""
import csv
import datetime
import hashlib
import zipfile
from pathlib import Path

import requests

RAW = Path("data/raw")
HUTUBS_URL = "https://sofacoustics.org/data/database/hutubs/pp{}_HRIRs_measured.sofa"
HUTUBS_SUBJECTS = range(2, 96)  # pp1 and pp96 are the FABIAN dummy head

SADIE_RECORD = "https://zenodo.org/api/records/12092466"
SADIE_URL = "https://zenodo.org/records/12092466/files/{}?download=1"
SADIE_SUBJECTS = [f"H{i}" for i in range(3, 21)]  # human subjects; D1 and D2 are dummy heads


def file_hash(path, algo):
    h = hashlib.new(algo)
    h.update(path.read_bytes())
    return h.hexdigest()


def download(url, path):
    if not path.exists():
        print("downloading", url)
        path.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path


rows = []
today = datetime.date.today().isoformat()

for s in HUTUBS_SUBJECTS:
    url = HUTUBS_URL.format(s)
    path = download(url, RAW / "hutubs" / f"pp{s}_HRIRs_measured.sofa")
    rows.append(["HUTUBS measured HRIRs (doi:10.14279/depositonce-8487)", f"pp{s}", url, path,
                 file_hash(path, "sha256"), today])

md5 = {f["key"]: f["checksum"].split(":")[1] for f in requests.get(SADIE_RECORD, timeout=60).json()["files"]}
for subject in SADIE_SUBJECTS:
    name = f"{subject}_HRIR_SOFA.zip"
    url = SADIE_URL.format(name)
    zip_path = download(url, RAW / "sadie" / name)
    assert file_hash(zip_path, "md5") == md5[name], f"md5 mismatch for {name}"
    with zipfile.ZipFile(zip_path) as z:
        sofa = [n for n in z.namelist() if n.endswith(".sofa") and "44" in Path(n).name]
        assert len(sofa) == 1, f"expected one 44.1 kHz SOFA file in {name}, found {sofa}"
        out = RAW / "sadie" / Path(sofa[0]).name
        out.write_bytes(z.read(sofa[0]))
    rows.append(["SADIE II v2-2", subject, url, out, file_hash(out, "sha256"), today])

with open("data/manifest.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["database", "subject", "url", "local_file", "sha256", "downloaded"])
    w.writerows(rows)
print("wrote data/manifest.csv")
