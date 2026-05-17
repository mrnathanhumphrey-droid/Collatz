"""
compute_eps14_from_existing.py
==============================
Compute ε_14 from pre-existing pi_14_truncated.npz via FFT (same approach as ε_13).
"""
from __future__ import annotations

import csv
import sys
import time
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = Path(r"C:\Collatz\probe_epsilon_14")
OUT_DIR.mkdir(exist_ok=True)
PI_PATH = Path(r"C:\Collatz\probe_self_similarity\pi_14_truncated.npz")

EPS_KNOWN = {
    1:  +2.0000000000e-01,
    2:  +9.5238095238e-03,
    3:  -5.0919863259e-03,
    4:  -2.4522582483e-03,
    5:  -1.1517469151e-03,
    6:  -4.9790566522e-04,
    7:  -1.1752368304e-03,
    8:  -7.4554636729e-04,
    9:  -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
    12: +2.2747137206e-03,
    13: +2.9482473172e-03,
}


def main():
    print("=" * 70)
    print("epsilon_14 from existing pi_14_truncated.npz")
    print("=" * 70)

    if not PI_PATH.exists():
        print(f"ERROR: {PI_PATH} not found")
        return
    d = np.load(PI_PATH)
    pi14 = d["pi"]
    coprime14 = d["coprime"]
    v_max = int(d["v_max"])
    k = int(d["k"])
    n = len(pi14)
    print(f"\nLoaded pi_14 from {PI_PATH}")
    print(f"  k = {k}, n = {n:,}, v_max = {v_max}")
    print(f"  sum(pi_14) = {pi14.sum():.15f}")
    print(f"  truncation error bound: 2^-{v_max} = {2.0**-v_max:.2e}")

    assert abs(pi14.sum() - 1.0) < 1e-12, "pi_14 not normalized"
    assert n == 2 * 3 ** 13, f"wrong dim: expected {2*3**13}, got {n}"

    # === S_14 via FFT ===
    print(f"\nFFT cross-check of S_14...")
    t0 = time.time()
    N14 = 3 ** 14
    pi_full = np.zeros(N14, dtype=np.float64)
    pi_full[coprime14] = pi14
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N14)
    mask_nontrivial = xi_arr % 3 != 0
    S14_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps14_fft = S14_fft - 7.0 / 15.0
    t_fft = time.time() - t0
    print(f"  S_14 (FFT) = {S14_fft:.15f}  ({t_fft:.2f}s)")
    print(f"  eps_14 (FFT) = {eps14_fft:+.12e}")

    # Print full eps table
    print(f"\neps_k extended:")
    for kk in sorted(EPS_KNOWN):
        print(f"  eps_{kk} = {EPS_KNOWN[kk]:+.10e}")
    print(f"  eps_14 = {eps14_fft:+.10e}  ← NEW")

    # Envelope check
    env_13 = abs(EPS_KNOWN[13]) * 2**13
    env_14 = abs(eps14_fft) * 2**14
    print(f"\nEnvelope |eps_k| * 2^k:")
    print(f"  k=13: {env_13:.3e}")
    print(f"  k=14: {env_14:.3e}  (ratio {env_14/env_13:.4f})")

    # Hadamard radius update
    had_13 = 1.0 / abs(EPS_KNOWN[13])**(1/13)
    had_14 = 1.0 / abs(eps14_fft)**(1/14)
    print(f"\nHadamard radius |eps_k|^(-1/k):")
    print(f"  k=13: {had_13:.4f}")
    print(f"  k=14: {had_14:.4f}  (inward by {had_13-had_14:.4f})")

    # Ratio
    ratio = eps14_fft / EPS_KNOWN[13]
    print(f"\nRatio eps_14 / eps_13 = {ratio:+.6f}  (|.| = {abs(ratio):.6f})")
    print(f"  sign change: {'YES' if (eps14_fft > 0) != (EPS_KNOWN[13] > 0) else 'NO'}")

    # Save
    out = {
        "k": 14,
        "S_14": S14_fft,
        "epsilon_14": eps14_fft,
        "envelope_2to14": env_14,
        "hadamard_at_14": had_14,
        "ratio_eps14_eps13": ratio,
        "eps_known_through_13": EPS_KNOWN,
        "fft_time_sec": t_fft,
        "pi_path": str(PI_PATH),
        "v_max_truncation": v_max,
    }
    with open(OUT_DIR / "epsilon_14_result.json", "w") as f:
        json.dump(out, f, indent=2)
    with open(OUT_DIR / "S_14_epsilon_14.txt", "w") as f:
        f.write(f"k=14\nS_14 = {S14_fft:.20e}\nepsilon_14 = {eps14_fft:+.20e}\n")
        f.write(f"|eps_14| * 2^14 = {env_14:.6e}\n")
        f.write(f"Hadamard radius at k=14: {had_14:.6f}\n")
    print(f"\nSaved: {OUT_DIR / 'epsilon_14_result.json'}")


if __name__ == "__main__":
    main()
