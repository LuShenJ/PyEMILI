from pathlib import Path
import os
import sys
import warnings

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TEST_DIR = Path(__file__).resolve().parent
OUT_DIR = TEST_DIR / "_basic_outputs"

sys.path.insert(0, str(ROOT))


def test_spec_basic():
    from pyemili.Spec import Spec_line_finding

    spectrum = TEST_DIR / "J0608_03may2019MagE.txt"
    Spec_line_finding(
        str(spectrum),
        length=50,
        percentile=30,
        check_continuum=False,
        check_lines=False,
        save_continuum=False,
        snr_threshold=30,
        prominence=8,
    )
    out = OUT_DIR / "J0608_03may2019MagE_linelist.txt"
    if not out.exists():
        raise AssertionError(f"Spec_line_finding did not create {out}")
    print(f"Spec.py basic test passed: {out.name}")


def test_lines_basic():
    from pyemili.Lines import Line_list

    hf22 = np.loadtxt(TEST_DIR / "Hf2-2_linelist.txt", skiprows=1, max_rows=5)
    line_list = Line_list(
        wavelength=hf22[:, 0],
        wavelength_error=10,
        flux=hf22[:, 1],
        flux_error=hf22[:, 1] * hf22[:, 2] * 0.01,
        snr=hf22[:, 3],
        fwhm=hf22[:, 4],
    )
    line_list.identify(
        str(OUT_DIR / "Hf2-2_basic"),
        icf=[0.01, 0.5, 0.4, 0.1, 0.0001],
        v_cor=[0, 0, 0, 0, 0],
        iteration=False,
        abun_type="nebula",
        erc_list=True,
    )
    for suffix in (".out", ".dat"):
        out = OUT_DIR / f"Hf2-2_basic{suffix}"
        if not out.exists():
            raise AssertionError(f"Line_list.identify did not create {out}")
    print("Lines.py basic test passed: Hf2-2_basic.out/dat")


def test_recomb_basic():
    from pyemili.recomb import Recom_Lines

    ic = Recom_Lines(str(TEST_DIR / "Abell46_erc.dat"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = ic.fit(
            ion="O II",
            fit_abundance=True,
            flux_threshold=0.0,
            vel=20.0,
            use_lab_wav=True,
            nwalkers=32,
            max_iter=300,
            progress=False,
            seed=5,
            make_corner=False,
            make_scatter=False,
        )
    if not result.medians:
        raise AssertionError("Recom_Lines.fit returned empty medians")
    print(f"recomb.py basic test passed: {result.medians}")


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    os.chdir(OUT_DIR)

    tests = [
        ("Spec.py", test_spec_basic),
        ("Lines.py", test_lines_basic),
        ("recomb.py", test_recomb_basic),
    ]
    for name, func in tests:
        print(f"\n=== Testing {name} ===")
        func()
