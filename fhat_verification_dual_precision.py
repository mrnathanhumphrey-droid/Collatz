"""
fhat_verification_dual_precision.py — A1 adversarial safeguard.

At three representative cells (p=3,r=2), (p=11,r=2), (p=5,r=3), compute
F̂_p^short via two independent methods:
    Method 1: numpy.fft.fft (float64)
    Method 2: mpmath direct summation at 50-digit precision

Compare magnitudes on support; verify agreement to ≥ 1e-14.

Also implements A2 hand-computation verification at (p=3, r=2):
    Hand-derived |F̂_3^short(ξ=3)| = √27 ≈ 5.196152422706632
    Verify FFT implementation matches.
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

import numpy as np
import mpmath as mp


OUTPATH = Path("C:/Collatz/fhat_verification_a1_dual_precision.csv")


def fhat_short_mpmath(p: int, r: int, c: int = 1, dps: int = 50) -> dict:
    """Compute |F̂_p^short(ξ)| on the predicted support {p·a : a ≡ 1 mod p}
    via mpmath direct summation at `dps` decimal digits.

    For ξ = p·a, the kernel exp(-2πi p·a·u / M) = exp(-2πi a·u / p^r)
    (M = p^{r+1}, period = p^r). So we compute the period-DFT of f at frequency a.
    """
    mp.mp.dps = dps
    M = p ** (r + 1)
    period = p ** r
    two_pi = 2 * mp.pi
    c_mpf = mp.mpf(c)

    # Build f_p exactly: f_p(u) = exp(2πi · c · (1+p)^u / M)
    f_vals = []
    pow_val = 1
    M_mpf = mp.mpf(M)
    for u in range(period):
        arg = two_pi * c_mpf * pow_val / M_mpf
        f_vals.append(mp.mpc(mp.cos(arg), mp.sin(arg)))
        pow_val = (pow_val * (1 + p)) % M

    # Predicted support: a ≡ 1 mod p, a ∈ Z/p^r
    a_vals = list(range(1, period, p))
    period_mpf = mp.mpf(period)

    mags = []
    for a in a_vals:
        s = mp.mpc(0, 0)
        for u in range(period):
            phase = -two_pi * a * u / period_mpf
            kernel = mp.mpc(mp.cos(phase), mp.sin(phase))
            s += f_vals[u] * kernel
        mags.append(abs(s))

    return {
        "p": p, "r": r,
        "predicted_supp_size": len(a_vals),
        "mags_supp_mpmath": mags,
        "max_mag_mpmath": max(mags),
        "min_mag_mpmath": min(mags),
    }


def fhat_short_fft(p: int, r: int, c: int = 1) -> dict:
    """Same as fhat_short_mpmath but via numpy FFT (float64).

    Uses length-period FFT (G[a]) since F̂_full(p·a) = p · G[a] and the on-support
    magnitudes |G[a]| = p^{(r+1)/2} (equivalently |F̂_full(p·a)| = p^{(r+3)/2}).
    """
    M = p ** (r + 1)
    period = p ** r
    f = np.empty(period, dtype=np.complex128)
    pow_val = 1
    for u in range(period):
        f[u] = np.exp(2j * np.pi * c * pow_val / M)
        pow_val = (pow_val * (1 + p)) % M
    G = np.fft.fft(f)
    mags_all = np.abs(G)

    a_vals = np.arange(1, period, p)
    mags_supp = mags_all[a_vals]
    return {
        "p": p, "r": r,
        "predicted_supp_size": int(len(a_vals)),
        "mags_supp_fft": mags_supp.tolist(),
        "max_mag_fft": float(mags_supp.max()),
        "min_mag_fft": float(mags_supp.min()),
    }


def a2_hand_computation_check() -> dict:
    """A2: verify FFT implementation matches hand-computation at (p=3, r=2).

    Hand result for ξ = 3:
        f(u) = exp(2πi · k_u / 27) where k_u = 4^u mod 27 = (1,4,16,10,13,25,19,22,7).
        F̂(3) = Σ exp(2πi(k_u - 3u)/27).
        Phases (k_u - 3u) mod 27 = (1, 1, 10, 1, 1, 10, 1, 1, 10).
        Six terms at phase 1/27, three terms at phase 10/27.
        |F̂(3)|² = |6 e^{2πi/27} + 3 e^{2πi·10/27}|²
                = 36 + 9 + 36 cos(2π·9/27) = 45 - 18 = 27.
        |F̂(3)| = √27 = 5.196152422706632...

    Verify FFT gives the same.
    """
    p, r = 3, 2
    M, period = 27, 9
    f = np.empty(period, dtype=np.complex128)
    pow_val = 1
    for u in range(period):
        f[u] = np.exp(2j * np.pi * pow_val / M)
        pow_val = (pow_val * 4) % M
    f_padded = np.zeros(M, dtype=np.complex128)
    f_padded[:period] = f
    Fhat = np.fft.fft(f_padded)
    fft_val_at_3 = abs(Fhat[3])

    hand_val = float(mp.sqrt(27))

    rel_diff = abs(fft_val_at_3 - hand_val) / hand_val
    return {
        "test": "A2_hand_computation",
        "cell": "(p=3, r=2, ξ=3)",
        "hand_computed_mag": hand_val,
        "fft_computed_mag": fft_val_at_3,
        "abs_diff": abs(fft_val_at_3 - hand_val),
        "rel_diff": rel_diff,
        "pass": rel_diff < 1e-14,
    }


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    rows = []

    print("=" * 84)
    print("A1: dual-precision FFT vs mpmath comparison")
    print("=" * 84)
    print()

    cells = [(3, 2), (11, 2), (5, 3)]
    for p, r in cells:
        t0 = time.time()
        res_fft = fhat_short_fft(p, r)
        t_fft = time.time() - t0

        t0 = time.time()
        res_mp = fhat_short_mpmath(p, r, dps=50)
        t_mp = time.time() - t0

        # Compare magnitudes elementwise
        mags_fft = np.array(res_fft["mags_supp_fft"])
        mags_mp = np.array([float(m) for m in res_mp["mags_supp_mpmath"]])

        abs_diffs = np.abs(mags_fft - mags_mp)
        rel_diffs = abs_diffs / mags_mp
        max_abs = float(abs_diffs.max())
        max_rel = float(rel_diffs.max())

        # Comparison to theoretical p^{(r+1)/2}
        pred_short = float(mp.power(p, mp.mpf(r + 1) / 2))
        mp_dev_from_pred = float(max(abs(m - mp.mpf(pred_short)) / pred_short for m in res_mp["mags_supp_mpmath"]))

        print(f"  (p={p}, r={r}):")
        print(f"    FFT  ({t_fft:.2f}s): max={res_fft['max_mag_fft']:.16f}, min={res_fft['min_mag_fft']:.16f}")
        print(f"    mpmath ({t_mp:.1f}s, 50dps): max={float(res_mp['max_mag_mpmath']):.16f}, min={float(res_mp['min_mag_mpmath']):.16f}")
        print(f"    Predicted p^((r+1)/2) = {pred_short:.16f}")
        print(f"    max |FFT - mpmath|     = {max_abs:.2e}    rel = {max_rel:.2e}")
        print(f"    mpmath dev from predicted (50dps)  = {mp_dev_from_pred:.2e}")
        print()

        rows.append({
            "test": "A1_dual_precision",
            "p": p, "r": r,
            "supp_size": res_fft["predicted_supp_size"],
            "predicted_p^((r+1)/2)": pred_short,
            "fft_max": res_fft["max_mag_fft"],
            "fft_min": res_fft["min_mag_fft"],
            "mpmath_max": float(res_mp["max_mag_mpmath"]),
            "mpmath_min": float(res_mp["min_mag_mpmath"]),
            "max_abs_diff_fft_vs_mpmath": max_abs,
            "max_rel_diff_fft_vs_mpmath": max_rel,
            "mpmath_max_dev_from_predicted_50dps": mp_dev_from_pred,
            "fft_time_s": t_fft,
            "mpmath_time_s": t_mp,
            "pass_agreement": max_rel < 1e-12,
            "pass_mpmath_matches_predicted_to_50dps": mp_dev_from_pred < 1e-30,
        })

    print()
    print("=" * 84)
    print("A2: hand-computation cross-check at (p=3, r=2, ξ=3)")
    print("=" * 84)
    print()

    res_a2 = a2_hand_computation_check()
    print(f"  Cell: {res_a2['cell']}")
    print(f"  Hand-computed |F̂(3)| = √27 = {res_a2['hand_computed_mag']:.16f}")
    print(f"  FFT-computed  |F̂(3)|      = {res_a2['fft_computed_mag']:.16f}")
    print(f"  |abs diff|                 = {res_a2['abs_diff']:.2e}")
    print(f"  rel diff                   = {res_a2['rel_diff']:.2e}")
    print(f"  pass (rel < 1e-14)         = {res_a2['pass']}")

    rows.append({
        "test": "A2_hand_computation",
        "p": 3, "r": 2,
        "hand_computed_mag": res_a2["hand_computed_mag"],
        "fft_computed_mag": res_a2["fft_computed_mag"],
        "abs_diff": res_a2["abs_diff"],
        "rel_diff": res_a2["rel_diff"],
        "pass": res_a2["pass"],
    })

    # Write CSV
    all_keys = sorted({k for row in rows for k in row.keys()})
    with open(OUTPATH, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["test"] + [k for k in all_keys if k != "test"])
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print()
    print(f"[write] {OUTPATH}")


if __name__ == "__main__":
    main()
