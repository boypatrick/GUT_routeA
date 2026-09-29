#!/usr/bin/env python3
"""Read-only source audit; no thermal fit and no synthetic observations.

Run with --rocit-zip to audit every E2/E3 file in the pinned full release.
Without it, audit the exact first-day excerpt shipped in the repository.
Full and excerpt results use different filenames to preserve provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/clock_comparison_audit"
PAIR = "PTB_Yb1E2_CombYb-PTB_Yb_CombKnoten"
PTB_SHA = "3f38135ebfab1375391195c282aa1980a76657ba9125b59e6e1cb85695a400e6"
ROCIT_SHA = "6168e24a0c29ce0929e9651460f11ea77f151176e1ba3d0fe3428e1e08bd56bd"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def profile(blob):
    text = blob.decode("utf-8-sig")
    a = np.loadtxt(io.StringIO(text), comments="#", ndmin=2)
    assert a.shape[1] == 5, a.shape
    valid = a[:, 2] != 0
    times = a[:, 0]
    return a, dict(
        sha256=sha(blob), bytes=len(blob), rows=len(a), columns=a.shape[1],
        headers=[s for s in text.splitlines() if s.startswith("#")],
        nonfinite_values=int((~np.isfinite(a)).sum()),
        duplicate_timestamps=int(len(a)-len(np.unique(times))),
        nonincreasing_steps=int((np.diff(times) <= 0).sum()),
        time_range_mjd=[float(times.min()), float(times.max())],
        flags={str(int(k)): int(v) for k,v in zip(*np.unique(a[:, 2], return_counts=True))},
        valid_rows=int(valid.sum()),
        uncertainty_A_values=np.unique(a[:, 3]).tolist(),
        uncertainty_B_values=np.unique(a[:, 4]).tolist(),
        nominal_second_steps_outside_rounding=int(
            (abs(np.diff(times)*86400-1) > .1).sum()),
        thermal_columns_present=False,
        note="Five released fields are time, comparator, validity, uA_sys, uB_sys; not raw clock servo plus radiometry.",
    )


def run(full_path=None):
    ptb_bytes=(DATA / "ptb_lpi_2021.zip").read_bytes()
    assert sha(ptb_bytes) == PTB_SHA
    with zipfile.ZipFile(io.BytesIO(ptb_bytes)) as z:
        p = z.read("Data/E3_E2_frequency_ratio_measurement.txt")
        a = np.loadtxt(io.StringIO(p.decode()), comments="#")
        ptb = dict(
            archive_sha256=sha(ptb_bytes), archive_members=z.namelist(),
            ratio_sha256=sha(p), rows=len(a), columns=a.shape[1],
            fields=["MJD", "E3/E2 minus 0.932829404530965 (absolute ratio difference)",
                    "statistical uncertainty", "reproducibility", "combined uncertainty"],
            nonfinite_values=int((~np.isfinite(a)).sum()),
            duplicate_timestamps=int(len(a)-len(np.unique(a[:,0]))),
            time_range_mjd=[float(a[:,0].min()), float(a[:,0].max())],
            uncertainty_quadrature_max_relative_error=float(np.max(
                abs(np.hypot(a[:,2],a[:,3])-a[:,4])/a[:,4])),
            ratio_is_fractional=False,
            eligibility="ineligible: published comparison points, no paired bath/radiometry or correction ledger; paper uses two clock systems",
        )
    profiles = {}
    arrays = []
    metadata = {}
    if full_path is not None:
        blob=Path(full_path).read_bytes()
        assert sha(blob) == ROCIT_SHA, "Full archive differs from reviewed version"
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            # No archive paths are executed or extracted here.
            members=z.namelist()
            names=sorted(n for n in members if n.startswith(PAIR+"/") and n.endswith(".dat"))
            metadata={n:z.read(n).decode() for n in members if n.endswith(".yml")}
            for n in names:
                a, p=profile(z.read(n))
                profiles[n]=p
                arrays.append(a)
        scope="full_selected_pair"
    else:
        names=sorted((DATA / "rocit_2025").glob("*.dat"))
        members=[p.name for p in (DATA / "rocit_2025").iterdir()]
        for p in names:
            a, result=profile(p.read_bytes())
            profiles[p.name]=result
            arrays.append(a)
        metadata={p.name:p.read_text() for p in (DATA / "rocit_2025").glob("*.yml")}
        scope="bundled_first_day_only"
    a=np.concatenate(arrays)
    pair_metadata=[s for n,s in metadata.items() if PAIR in n]
    assert len(pair_metadata)==1
    rocit=dict(
        doi="10.5281/zenodo.17107693", license="CC-BY-4.0",
        full_archive_sha256=ROCIT_SHA, scope=scope,
        all_archive_members=members, all_yaml_metadata=metadata,
        selected_pair=PAIR, selected_files=len(profiles), rows=len(a),
        valid_rows=int((a[:,2]!=0).sum()),
        invalid_rows=int((a[:,2]==0).sum()),
        nonfinite_values=int((~np.isfinite(a)).sum()),
        duplicate_timestamps=int(len(a)-len(np.unique(a[:,0]))),
        time_range_mjd=[float(a[:,0].min()),float(a[:,0].max())],
        metadata_for_pair=pair_metadata[0], files=profiles,
        direction="A=E3, B=E2; released Delta_A_to_B is a scaled transfer beat, not exactly E3/E2 fractional shift",
        eligibility="ineligible for requested thermal test: no two-temperature assignment, independent radiation calibration, or per-run BBR correction ledger in selected files/metadata",
    )
    out=dict(stage="G-R5 data acquisition audit", date="2026-09-24",
             ptb_2021=ptb, rocit_2025=rocit,
             empirical_differential_test_ready=False,
             common_time_effect_measured=False,
             conclusion="Public comparison data obtained; paired uncorrected two-bath-state clock/radiometry data not located in audited releases. No thermal slope fit.")
    path=ROOT / "output" / ("gr5_clock_data_full_audit.json" if full_path else "gr5_clock_data_excerpt_audit.json")
    path.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(dict(scope=scope, ptb_points=ptb['rows'],
                         rocit_files=rocit['selected_files'], rocit_rows=rocit['rows'],
                         valid=rocit['valid_rows'],duplicates=rocit['duplicate_timestamps'],
                         range_mjd=rocit['time_range_mjd'], empirical_test_ready=False)))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--rocit-zip",type=Path)
    run(parser.parse_args().rocit_zip)
