"""
result_qspectrum.py
===================
Top 10 eigenvalues of K_k^(q) at q in {3, 5, 7, 11, 13}.

State-count plan (based on actual phi(q^k) = (q-1)*q^(k-1)):
  q=3,  k=4:    54 states   (smallest available; baseline)
  q=5,  k=4:   500 states
  q=7,  k=4:  1029 states
  q=11, k=3:  1210 states   (k=4 would be 13310, ~1.4 GB)
  q=13, k=3:  2028 states   (k=4 would be 26364, ~5.6 GB)

Method: build dense K in float64 via the same machinery used at q=3 k=6/7.
Compute full spectrum via scipy.linalg.eig, sort by magnitude, report top 10.

Outputs:
  result_qspectrum.md   - tables, q-universality verdict
  result_qspectrum.csv  - lambda_i for i=1..10 across q's
"""
from __future__ import annotations
import csv
import os
import sys
import time

import numpy as np
import scipy.linalg as la

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"


def order_of_two(N):
    assert N % 2 == 1
    m = 1; v = 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def build_K_float(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    p = inv2
    for v in range(M):
        powers_inv2[v] = p
        p = (p * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    K = np.zeros((n, n), dtype=np.float64)
    Z_v = 1.0 - 2.0 ** (-M)
    weights = np.array([(2.0 ** (-v)) / Z_v for v in range(1, M + 1)],
                       dtype=np.float64)
    for i_r, r in enumerate(coprime):
        base = (q * int(r) + 1) % N
        for v in range(1, M + 1):
            tgt = (base * int(powers_inv2[v - 1])) % N
            j = int(state_idx[tgt])
            K[i_r, j] += weights[v - 1]
    return K, M, n


def main():
    print("=" * 78)
    print("Top 10 eigenvalues of K_k^(q) across q in {3, 5, 7, 11, 13}")
    print("=" * 78)

    configs = [
        (3, 4),
        (5, 4),
        (7, 4),
        (11, 3),
        (13, 3),
    ]

    results = {}
    for q, k in configs:
        print(f"\n--- q={q}, k={k} ---")
        t0 = time.time()
        K, M, n = build_K_float(q, k)
        t_build = time.time() - t0
        print(f"  states={n}, M=ord_q^k(2)={M}, build {t_build:.2f}s")
        # Sanity: row sums should equal 1 to machine precision
        row_sum_dev = float(np.max(np.abs(K.sum(axis=1) - 1)))
        print(f"  row-sum check: max |row sum - 1| = {row_sum_dev:.2e}")
        t0 = time.time()
        eigs = la.eigvals(K)
        t_eig = time.time() - t0
        print(f"  eigenvalues: {t_eig:.2f}s")
        eigs_sorted = sorted(eigs, key=lambda x: -abs(x))
        top10 = eigs_sorted[:10]
        results[(q, k)] = {
            "n": n, "M": M, "t_build": t_build, "t_eig": t_eig,
            "top10": top10,
        }
        print(f"  Top 10 |lambda_i|:")
        for i, lam in enumerate(top10, 1):
            print(f"    lambda_{i} = {lam.real:+.10f} + {lam.imag:+.10f}i  "
                  f"|.|={abs(lam):.10f}")

    # Comparison table
    print()
    print("=" * 78)
    print("Cross-q spectrum comparison")
    print("=" * 78)
    print()
    print(f"  {'q':>3}  {'k':>3}  {'n':>5}  "
          + "  ".join(f"|λ_{i}|".rjust(10) for i in range(1, 11)))
    for (q, k), info in results.items():
        mods = [abs(l) for l in info["top10"]]
        line = f"  {q:>3}  {k:>3}  {info['n']:>5}  " + "  ".join(
            f"{m:>10.6f}" for m in mods)
        print(line)

    # |lambda_2| comparison highlighting
    print()
    print("|lambda_2|^(q) (the rate-controlling eigenvalue):")
    for (q, k), info in results.items():
        lam2 = info["top10"][1]
        ref_half = abs(abs(lam2) - 0.5)
        ref_inv_q = abs(abs(lam2) - 1.0/q)
        print(f"  q={q:>3}: |λ_2| = {abs(lam2):.6f}  "
              f"vs 1/2 (Δ={ref_half:.4f}), vs 1/q (Δ={ref_inv_q:.4f})")

    # Gap structure
    print()
    print("Spectral gaps |λ_i| - |λ_{i+1}|:")
    print(f"  {'q':>3}  " + "  ".join(f"gap_{i}".rjust(8) for i in range(1, 6)))
    for (q, k), info in results.items():
        mods = [abs(l) for l in info["top10"]]
        gaps = [mods[i] - mods[i+1] for i in range(5)]
        line = f"  {q:>3}  " + "  ".join(f"{g:>8.5f}" for g in gaps)
        print(line)

    # CSV
    out_csv = os.path.join(OUT_DIR, "result_qspectrum.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["q", "k", "n_states", "M_qk",
                    "i", "lambda_real", "lambda_imag", "lambda_abs"])
        for (q, k), info in results.items():
            for i, lam in enumerate(info["top10"], 1):
                w.writerow([q, k, info["n"], info["M"],
                            i, lam.real, lam.imag, abs(lam)])
    print()
    print(f"saved {out_csv}")

    # Markdown
    md = []
    md.append("# Result: top 10 eigenvalues of K_k^(q) across q ∈ {3,5,7,11,13}")
    md.append("")
    md.append("**Date:** 2026-05-05.  Float64 dense eigensolve (scipy.linalg.eig).")
    md.append("")
    md.append("## Configuration")
    md.append("")
    md.append("| q | k | states | M = ord_{q^k}(2) | build (s) | eig (s) |")
    md.append("|---|---|---|---|---|---|")
    for (q, k), info in results.items():
        md.append(f"| {q} | {k} | {info['n']} | {info['M']} | "
                  f"{info['t_build']:.2f} | {info['t_eig']:.2f} |")
    md.append("")
    md.append("## Top 10 |λ_i| across q")
    md.append("")
    md.append("| q | k | " + " | ".join(f"\\|λ_{i}\\|" for i in range(1, 11)) + " |")
    md.append("|---|---|" + "---|" * 10)
    for (q, k), info in results.items():
        mods = [abs(l) for l in info["top10"]]
        md.append(f"| {q} | {k} | " + " | ".join(f"{m:.6f}" for m in mods) + " |")
    md.append("")
    md.append("## λ_2 (rate-controlling) per q")
    md.append("")
    md.append("| q | λ_2 (real) | λ_2 (imag) | |λ_2| | Δ from 1/2 | Δ from 1/q |")
    md.append("|---|---|---|---|---|---|")
    for (q, k), info in results.items():
        lam2 = info["top10"][1]
        ref_half = abs(abs(lam2) - 0.5)
        ref_inv_q = abs(abs(lam2) - 1.0/q)
        md.append(f"| {q} | {lam2.real:+.6f} | {lam2.imag:+.6f} | "
                  f"{abs(lam2):.6f} | {ref_half:.4f} | {ref_inv_q:.4f} |")
    md.append("")
    md.append("## Spectral gaps")
    md.append("")
    md.append("| q | " + " | ".join(f"\\|λ_{i}\\| − \\|λ_{i+1}\\|" for i in range(1, 6)) + " |")
    md.append("|---|" + "---|" * 5)
    for (q, k), info in results.items():
        mods = [abs(l) for l in info["top10"]]
        gaps = [mods[i] - mods[i+1] for i in range(5)]
        md.append(f"| {q} | " + " | ".join(f"{g:.5f}" for g in gaps) + " |")
    md.append("")
    md.append("## Verdict")
    md.append("")
    # Auto-classify
    lam2s = [abs(results[(q, k)]["top10"][1]) for q, k in configs]
    spread_lam2 = max(lam2s) - min(lam2s)
    median_lam2 = sorted(lam2s)[len(lam2s) // 2]
    md.append(f"- |λ_2| values across q: " + ", ".join(
        f"q={q}: {abs(results[(q, k)]['top10'][1]):.4f}" for q, k in configs))
    md.append(f"- Spread (max - min) = {spread_lam2:.4f}")
    md.append(f"- Median = {median_lam2:.4f}")
    md.append("")
    if spread_lam2 < 0.05:
        verdict = (f"|λ_2| is approximately q-UNIVERSAL at ≈ {median_lam2:.4f} "
                   f"(spread {spread_lam2:.4f}). Spectral structure has q-independent "
                   f"rate-controlling mode.")
    elif spread_lam2 < 0.15:
        verdict = (f"|λ_2| is mildly q-DEPENDENT (spread {spread_lam2:.4f}). "
                   f"Spectrum has a q-dependent component but still "
                   f"clustered around {median_lam2:.4f}.")
    else:
        verdict = (f"|λ_2| is strongly q-DEPENDENT (spread {spread_lam2:.4f}). "
                   f"No q-universal rate; each q has its own spectral structure.")
    md.append(f"**{verdict}**")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_qspectrum.py` — script")
    md.append("- `result_qspectrum.csv` — top 10 eigenvalues for each q")
    md.append("- `result_qspectrum.md` — this writeup")

    out_md = os.path.join(OUT_DIR, "result_qspectrum.md")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"saved {out_md}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
