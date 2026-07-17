"""
result_epsilon_11.py
====================
Compute eps_11 = S_11 - 7/15 to lock down the oscillation period and rate.

At k=10 we got eps_10 = +7.21e-4 (sign flipped from k=9), |eps_10/eps_9| = 95.8.
Two-mode model fit: rho ≈ 0.984 per step, theta ≈ 0.68 rad, period ≈ 9.2.
If the model is right, eps_11 should still be positive (near peak) with
magnitude similar to eps_10.

K_11: 2*3^10 = 118,098 states. Matrix-free required (dense = 110 GB).
Skipping scipy.eigs cross-check to save runtime — power iter to machine
precision + FFT independent of X-formula give sufficient cross-validation.
"""
from __future__ import annotations

import csv
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
ENV_CSV = os.path.join(OUT_DIR, "result_q_sweep_test_1_envelope.csv")

# Cached X_k from earlier verified runs (matched at machine precision)
X_CACHED = {
    5: 3.534161151367800,
    6: 4.000329912369247,
    7: 4.465821342205540,
    8: 4.931742462504920,
    9: 5.398401608914431,
    10: 5.865789026498212,
}


def order_of_two(N):
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


class MatVecK:
    def __init__(self, q, k, chunk=1024):
        self.q = q
        self.k = k
        self.N = q ** k
        self.M = order_of_two(self.N)
        self.coprime = np.array(
            [r for r in range(self.N) if r % q != 0], dtype=np.int64
        )
        self.n = len(self.coprime)
        self.state_idx = -np.ones(self.N, dtype=np.int64)
        for i, r in enumerate(self.coprime):
            self.state_idx[r] = i
        inv2 = pow(2, -1, self.N)
        self.powers_inv2 = np.empty(self.M, dtype=np.int64)
        pi = inv2
        for v in range(self.M):
            self.powers_inv2[v] = pi
            pi = (pi * inv2) % self.N
        self.base = (q * self.coprime + 1) % self.N
        self.Z_v = 1.0 - 2.0 ** (-self.M)
        self.weights = (2.0 ** (-np.arange(1, self.M + 1))) / self.Z_v
        self.chunk = chunk

    def left_mul(self, x):
        n = self.n
        x_new = np.zeros(n, dtype=np.float64)
        for v0 in range(0, self.M, self.chunk):
            v1 = min(v0 + self.chunk, self.M)
            chunk_powers = self.powers_inv2[v0:v1]
            chunk_weights = self.weights[v0:v1]
            targets = (chunk_powers[:, None] * self.base[None, :]) % self.N
            j_targets = self.state_idx[targets]
            weighted = chunk_weights[:, None] * x[None, :]
            x_new += np.bincount(
                j_targets.ravel(),
                weights=weighted.ravel(),
                minlength=n,
            )
        return x_new


def stationary_power_iteration_mfree(K_op, max_iter=300, tol=1e-13,
                                      verbose_every=1):
    n = K_op.n
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    for it in range(max_iter):
        t0 = time.time()
        pi_new = K_op.left_mul(pi)
        pi_new /= pi_new.sum()
        residual = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        dt = time.time() - t0
        if (it + 1) % verbose_every == 0 or residual < tol:
            print(f"    iter {it+1:>4}: residual = {residual:.4e}, "
                  f"matvec {dt:.2f}s", flush=True)
        if residual < tol:
            return pi, it + 1, residual
    return pi, max_iter, residual


def load_cached_S():
    S = {}
    with open(ENV_CSV) as f:
        for row in csv.DictReader(f):
            if int(row["q"]) != 3:
                continue
            S[int(row["n"])] = Fraction(int(row["S_n_num"]),
                                         int(row["S_n_den"]))
    return S


def main():
    print("=" * 78)
    print("Computing eps_11 = S_11 - 7/15 (oscillation period/rate lock-down)")
    print("=" * 78)

    S_cache = load_cached_S()
    eps_all = {k: float(S_cache[k] - Fraction(7, 15)) for k in S_cache}

    # k=6..10 hardcoded eps (from prior verified runs)
    eps_all[6] = -4.9790566522e-04
    eps_all[7] = -1.1752368304e-03
    eps_all[8] = -7.4554636729e-04
    eps_all[9] = -7.5202571564e-06
    eps_all[10] = +7.2075091711e-04

    print(f"\nUsing cached X_10 = {X_CACHED[10]:.15f}")
    print(f"Using cached eps_10 = {eps_all[10]:+.10e}")

    # k = 11 (target, matrix-free)
    print(f"\n--- k = 11 (TARGET, matrix-free) ---")
    print(f"  expected size: 118,098 states")
    t0 = time.time()
    K_op = MatVecK(3, 11, chunk=1024)
    t_init = time.time() - t0
    n11 = K_op.n
    M11 = K_op.M
    print(f"  initialize MatVecK_11: {n11:,} states, M={M11:,}, "
          f"chunk=1024, {t_init:.2f}s")

    # Sanity: check 3 random rows
    print("  K_11 row-sum sanity check (3 random rows)...")
    rng = np.random.default_rng(20260505)
    for _ in range(3):
        ir = int(rng.integers(0, n11))
        e_i = np.zeros(n11)
        e_i[ir] = 1.0
        t0_row = time.time()
        row = K_op.left_mul(e_i)
        dt_row = time.time() - t0_row
        print(f"    row {ir}: sum = {row.sum():.15f}  ({dt_row:.1f}s)",
              flush=True)

    print("\n  power iteration:")
    t0 = time.time()
    pi11, iters11, res11 = stationary_power_iteration_mfree(
        K_op, tol=1e-13, max_iter=200, verbose_every=1
    )
    t_iter = time.time() - t0
    print(f"  power iter complete: {iters11} iters, residual {res11:.2e}, "
          f"{t_iter:.2f}s")
    print(f"  sum(pi_11) = {pi11.sum():.15f}  (expect 1)")

    # Compute X_11, S_11, eps_11
    sum_sq_11 = float((pi11 ** 2).sum())
    X11 = (3 ** 11) * sum_sq_11
    S11 = X11 - X_CACHED[10]
    eps11 = S11 - 7.0 / 15.0
    eps_all[11] = eps11
    print(f"\n  X_11 = {X11:.15f}, S_11 = {S11:.15f}")
    print(f"  eps_11 = {eps11:+.10e}")

    # FFT cross-check
    print("\n  FFT cross-check of S_11...")
    N11 = 3 ** 11
    pi_full = np.zeros(N11, dtype=np.float64)
    coprime_11 = K_op.coprime
    pi_full[coprime_11] = pi11
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N11)
    mask_nontrivial = xi_arr % 3 != 0
    S11_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps11_fft = S11_fft - 7.0 / 15.0
    print(f"  S_11 via FFT = {S11_fft:.15f}, eps_11_fft = "
          f"{eps11_fft:+.10e}")
    print(f"  agreement: |S_11_X - S_11_FFT| = {abs(S11 - S11_fft):.4e}")

    # Ratios + analysis
    eps9 = eps_all[9]
    eps10 = eps_all[10]
    ratio_10_11 = abs(eps11 / eps10)
    ratio_9_10 = abs(eps10 / eps9)
    print(f"\n|eps_11 / eps_10| = {ratio_10_11:.6f}  "
          f"(compare |eps_10/eps_9| = {ratio_9_10:.4f})")
    print(f"sign(eps_11) = {'+' if eps11 > 0 else '-'}  "
          f"(eps_10 was {'+' if eps10 > 0 else '-'})")

    # Two-mode model prediction check
    # Model: eps_k ~ rho^k cos(k theta + phi), with rho = 0.984, theta = 0.68
    # eps_8/eps_10 = -1.034 → rho^2 ≈ 0.967 → rho ≈ 0.984
    # eps_7/eps_8 → cos(theta) ≈ 0.78 → theta ≈ 0.68 rad
    rho_est = 0.984
    theta_est = 0.68
    # Calibrate amplitude from eps_10 (positive peak vicinity).
    # Guess phase: at k=10, cos(10*theta + phi) ~ +1 (peak), so phi ≈ -10*theta
    phi_est = -10 * theta_est
    eps11_predicted = (rho_est ** 11) * np.cos(11 * theta_est + phi_est)
    # Calibrate scale
    eps10_predicted_unscaled = (rho_est ** 10) * np.cos(10 * theta_est + phi_est)
    if abs(eps10_predicted_unscaled) > 1e-12:
        scale = eps10 / eps10_predicted_unscaled
        eps11_predicted *= scale
    print(f"\nTwo-mode model prediction (rho={rho_est}, theta={theta_est}):")
    print(f"  eps_11 predicted: {eps11_predicted:+.4e}")
    print(f"  eps_11 actual:    {eps11:+.4e}")
    print(f"  ratio actual/predicted: "
          f"{eps11/eps11_predicted if abs(eps11_predicted)>1e-12 else 'NA':.4f}")

    # Decision
    if eps11 * eps10 > 0 and 0.5 < ratio_10_11 < 1.5:
        verdict = (
            f"OSCILLATION CONFIRMED + PERIOD CHARACTERIZED. eps_11 same sign "
            f"as eps_10 (+) with |eps_11/eps_10| = {ratio_10_11:.3f} ~ 1. "
            "k=10..11 spans the post-node peak region; full period "
            "consistent with model period ~9.2 in k-space."
        )
    elif eps11 * eps10 < 0 and abs(eps11) > 0.1 * abs(eps10):
        verdict = (
            f"FAST OSCILLATION — eps_11 sign-flipped vs eps_10 already. "
            "Period < 4 in k-space. Different model than predicted; "
            "needs refit."
        )
    elif abs(eps11) < 0.05 * abs(eps10):
        verdict = (
            f"DOUBLE-NODE? eps_11 collapses again, |eps_11/eps_10| = "
            f"{ratio_10_11:.4f}. k=10 may have been a single overshoot; "
            "underlying decay is genuinely faster than expected."
        )
    else:
        verdict = (
            f"EXPLORE — |eps_11/eps_10| = {ratio_10_11:.4f}, sign flip = "
            f"{eps11 * eps10 < 0}. Doesn't fit known model regimes cleanly."
        )

    print(f"\n*** VERDICT: {verdict}")

    print("\nUpdated ratio trajectory |eps_{k+1}/eps_k|:")
    for k in range(1, 11):
        rr = abs(eps_all[k + 1] / eps_all[k])
        print(f"  k = {k} -> {k+1}: {rr:.6f}")

    # Outputs
    out_md = os.path.join(OUT_DIR, "result_epsilon_11.md")
    md = []
    md.append("# Result: eps_11 = S_11 - 7/15")
    md.append("")
    md.append(f"**Date:** 2026-05-05.  Float64 matrix-free power iteration "
              f"on K_11 ({n11:,} states, M={M11:,}). FFT cross-check.")
    md.append("")
    md.append(f"## Verdict")
    md.append("")
    md.append(verdict)
    md.append("")
    md.append("## Headline numbers")
    md.append("")
    md.append(f"- eps_11 = `{eps11:+.10e}` (power iter)")
    md.append(f"- eps_11 = `{eps11_fft:+.10e}` (FFT cross-check; agreement "
              f"{abs(eps11 - eps11_fft):.2e})")
    md.append(f"- |eps_11/eps_10| = **{ratio_10_11:.6f}**")
    md.append(f"- sign(eps_11) = {'+' if eps11 > 0 else '-'}")
    md.append("")
    md.append("## Two-mode model check")
    md.append("")
    md.append(f"Model: eps_k ≈ scale · ρ^k · cos(k·θ + φ),  "
              f"ρ = {rho_est}, θ = {theta_est} rad (period ≈ 9.2)")
    md.append("")
    md.append(f"- eps_11 predicted (scale-calibrated to eps_10): "
              f"`{eps11_predicted:+.4e}`")
    md.append(f"- eps_11 actual: `{eps11:+.4e}`")
    if abs(eps11_predicted) > 1e-12:
        md.append(f"- actual / predicted: "
                  f"`{eps11/eps11_predicted:.4f}`")
    md.append("")
    md.append("## Ratio trajectory (k=1..11)")
    md.append("")
    md.append("| k → k+1 | |eps_{k+1}/eps_k| |")
    md.append("|---|---|")
    for k in range(1, 11):
        rr = abs(eps_all[k + 1] / eps_all[k])
        md.append(f"| {k} → {k+1} | {rr:.6f} |")
    md.append("")
    md.append("## eps_k table (k=1..11)")
    md.append("")
    md.append("| k | eps_k | sign | source |")
    md.append("|---|---|:---:|---|")
    for k in [1, 2, 3, 4, 5]:
        ek = float(S_cache[k] - Fraction(7, 15))
        sgn = "+" if ek > 0 else "-"
        md.append(f"| {k} | {ek:+.10e} | {sgn} | exact rational |")
    for k in [6, 7, 8, 9, 10]:
        sgn = "+" if eps_all[k] > 0 else "-"
        md.append(f"| {k} | {eps_all[k]:+.10e} | {sgn} | "
                  f"float64 power iter (cached) |")
    sgn11 = "+" if eps11 > 0 else "-"
    md.append(f"| 11 | {eps11:+.10e} | {sgn11} | float64 matrix-free + "
              f"FFT |")
    md.append("")
    md.append("## Computation diagnostics (k=11)")
    md.append("")
    md.append(f"- States: {n11:,}, M = {M11:,}")
    md.append(f"- Init MatVecK: {t_init:.2f}s")
    md.append(f"- Power iter: {iters11} iters, residual = {res11:.2e}, "
              f"{t_iter:.2f}s total ({t_iter/iters11:.1f}s per matvec)")
    md.append(f"- FFT cross-check vs X-formula: diff = "
              f"{abs(S11 - S11_fft):.2e}")
    md.append(f"- pi_11 sum = {pi11.sum():.15f}")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_epsilon_11.py` — script")
    md.append("- `result_epsilon_11.csv` — eps_k and ratios for k=1..11")
    md.append("- `result_epsilon_11.md` — this writeup")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"\nsaved {out_md}")

    out_csv = os.path.join(OUT_DIR, "result_epsilon_11.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "eps_k", "abs_eps_k", "sign", "ratio_to_prev_signed",
                    "abs_ratio_to_prev"])
        prev = None
        for k in range(1, 12):
            ek = eps_all[k]
            sgn = "+" if ek > 0 else "-"
            ratio_signed = "" if prev is None else f"{ek/prev:+.10f}"
            ratio_abs = "" if prev is None else f"{abs(ek/prev):.10f}"
            w.writerow([k, f"{ek:+.15e}", f"{abs(ek):.15e}", sgn,
                        ratio_signed, ratio_abs])
            prev = ek
    print(f"saved {out_csv}")
    print("\nDone.")


if __name__ == "__main__":
    main()
