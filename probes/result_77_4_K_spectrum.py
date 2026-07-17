"""
result_77_4_K_spectrum.py
=========================
R77.4 follow-up (option 3): compute the direct spectrum of the level-k
Markov transition operator K_k for k = 3..6 and characterize its structure
near lambda = 1/2.

The R77.4 envelope curve-fitting on n = 2..6 gave verdict (M):
  - Jordan ruled out by direction (b < 0)
  - Logarithmic vs power-law indistinguishable on 5 points

This script bypasses curve-fitting and reads the eigenvalue structure of K_k
directly. Three signatures we look for:

  (J) Jordan at lambda = 1/2: K_k has eigenvalue 1/2 with algebraic
      multiplicity > geometric multiplicity at some k. (Very specific
      numerical signature; usually visible as repeated eigenvalues with
      defective Jordan structure.)

  (P) Branch-cut / point spectrum approaching 1/2: K_k has eigenvalues
      forming a curve approaching 1/2 from below as k grows. The envelope
      decays like n^{-alpha}/2^n.

  (L) Continuous spectrum touching 1/2: K_k has a band of eigenvalues
      whose density near 1/2 grows like a power of N_k. Envelope decays
      like log(n)/2^n or 1/(n log n)/2^n.

  (S) Simple eigenvalue at 1/2 with spectral gap below: rate exactly 1/2
      with constant prefactor. (Empirically inconsistent with the n=3 peak
      in the envelope, but worth checking.)

Output:
  - C:/Collatz/experiments_output/result_77_4_K_spectrum_data.csv
      All eigenvalues of K_k for k = 3, 4, 5, 6 (sorted by |lambda| desc).
  - Stdout: leading-spectrum summary, eigenvalues in band |lambda - 1/2| < 0.1,
            count of eigenvalues in concentric bands around 1/2.
"""
from __future__ import annotations

import csv
import os
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
OUT_CSV = os.path.join(OUTDIR, "result_77_4_K_spectrum_data.csv")


def build_markov_float(k: int):
    """Build K_k as a float64 matrix (right-stochastic, row-stochastic-by-target)."""
    N = 3 ** k
    M = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = np.zeros((n, n), dtype=np.float64)
    Z_v = (2 ** M - 1) / (2 ** M)
    for r in coprime_states:
        i = state_idx[r]
        for r_v in range(1, M + 1):
            p = (1.0 / 2 ** r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[i, state_idx[target]] += p
    return K, coprime_states


def main():
    print("=" * 78)
    print("R77.4 follow-up: direct spectrum of K_k for k = 3..6")
    print("=" * 78)
    print()

    rows_csv = []
    summaries = {}

    for k in [3, 4, 5, 6]:
        t0 = time.time()
        K, coprime = build_markov_float(k)
        n = K.shape[0]
        # K is right-stochastic acting on row-vectors via pi * K = pi.
        # Eigenvalues of K and K^T are identical; use K directly.
        lam = np.linalg.eigvals(K)
        # Sort by |lambda| descending
        idx = np.argsort(-np.abs(lam))
        lam = lam[idx]
        elapsed = time.time() - t0
        print(f"  k = {k}  N = {n:>4}  spec computed in {elapsed:.2f}s")
        print(f"    leading 8 eigenvalues (by |lambda|):")
        for j in range(min(8, len(lam))):
            l = lam[j]
            print(f"      lambda_{j+1} = {l.real:+.10f} {l.imag:+.10f}j   |lambda| = {abs(l):.10f}")
        print()

        # Eigenvalues in band |lambda - 1/2| < 0.1 (real-axis neighborhood)
        band = [l for l in lam if abs(l - 0.5) < 0.1]
        print(f"    eigenvalues in band |lambda - 0.5| < 0.1:  {len(band)}")
        for l in band:
            print(f"      {l.real:+.10f} {l.imag:+.10f}j   |l - 0.5| = {abs(l - 0.5):.6e}")
        print()

        # Concentric bands around 1/2 (count)
        bands = [0.001, 0.01, 0.05, 0.1, 0.2]
        print(f"    counts in concentric bands around lambda = 0.5:")
        for r in bands:
            cnt = sum(1 for l in lam if abs(l - 0.5) < r)
            print(f"      |l - 0.5| < {r:>6.3f}:  {cnt}")
        print()

        # Closest eigenvalue to 1/2
        diffs = np.abs(lam - 0.5)
        j_star = int(np.argmin(diffs))
        print(f"    closest to lambda = 0.5:  lambda = {lam[j_star]:+.10f}, "
              f"|l - 0.5| = {diffs[j_star]:.6e}")

        # Rank-2 Jordan signature: look for repeated eigenvalues at 1/2
        # (within numerical tolerance 1e-8). Geometric multiplicity check
        # would need full eigenvector analysis; here just count algebraic.
        repeated = sum(1 for l in lam if abs(l - 0.5) < 1e-8)
        print(f"    eigenvalues within 1e-8 of 0.5 (algebraic multiplicity proxy):  {repeated}")
        print()

        summaries[k] = {
            "N": n,
            "lam_top8": lam[:8].tolist(),
            "n_band_0p1": len(band),
            "lam_closest": complex(lam[j_star]),
            "diff_closest": float(diffs[j_star]),
            "n_within_1eneg8": repeated,
        }

        # Append CSV rows
        for j, l in enumerate(lam):
            rows_csv.append({
                "k": k,
                "rank": j + 1,
                "lambda_real": l.real,
                "lambda_imag": l.imag,
                "abs_lambda": float(abs(l)),
                "diff_to_half": float(abs(l - 0.5)),
            })

    # Write CSV
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["k", "rank", "lambda_real", "lambda_imag",
                                            "abs_lambda", "diff_to_half"])
        w.writeheader()
        w.writerows(rows_csv)
    print(f"[csv: {OUT_CSV}]")
    print()

    # -------- Cross-k scaling of distance to 1/2 -------- #
    print("=" * 78)
    print("Scaling: closest-eigenvalue distance to 1/2 vs k")
    print("=" * 78)
    print(f"  {'k':>3} {'N_k':>5} {'lambda_closest':>26} {'|l - 0.5|':>15}")
    for k in [3, 4, 5, 6]:
        s = summaries[k]
        l = s["lam_closest"]
        print(f"  {k:>3} {s['N']:>5} {l.real:+.10f}{l.imag:+.10f}j  {s['diff_closest']:.6e}")
    print()

    # If the closest eigenvalue's distance to 1/2 SHRINKS with k -> branch-cut/continuous-spec evidence
    # If it stays bounded away from 1/2 -> simple-eigenvalue evidence
    diffs_seq = [summaries[k]["diff_closest"] for k in [3, 4, 5, 6]]
    print(f"  diff sequence over k=3..6: {[f'{d:.4e}' for d in diffs_seq]}")
    if diffs_seq[-1] < 0.5 * diffs_seq[0]:
        print("  -> diff is SHRINKING with k: consistent with branch-cut / continuous-spectrum")
    elif diffs_seq[-1] > 1.5 * diffs_seq[0]:
        print("  -> diff is GROWING with k: unusual; spectrum drifting away from 1/2")
    else:
        print("  -> diff is APPROXIMATELY STABLE with k: consistent with simple-eigenvalue")
    print()

    # -------- Cross-k count of eigenvalues in band around 1/2 -------- #
    print("=" * 78)
    print("Scaling: count of eigenvalues with |l - 0.5| < 0.1 vs k (and N_k)")
    print("=" * 78)
    print(f"  {'k':>3} {'N_k':>5} {'count_band':>12}  {'count/N_k':>10}")
    for k in [3, 4, 5, 6]:
        s = summaries[k]
        cnt = s["n_band_0p1"]
        print(f"  {k:>3} {s['N']:>5} {cnt:>12}  {cnt / s['N']:>10.4f}")
    print()
    print("  If count grows ~N_k -> continuous spectrum filling [a, b] near 1/2")
    print("  If count stays bounded -> finite point spectrum near 1/2")
    print("  If count grows ~sqrt(N_k) or ~log(N_k) -> intermediate (unusual)")


if __name__ == "__main__":
    main()
