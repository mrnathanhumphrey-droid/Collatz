"""
result_epsilon_10.py
====================
Compute eps_10 = S_10 - 7/15 to test the oscillation hypothesis.

At k=9 we got eps_9 = -7.52e-6, a 100x collapse from eps_8 (|ratio| = 0.010).
This is too sharp for single-mode decay (pre-spike rates were 0.43-0.55), so
either (a) oscillation with a near-zero node at k=9, or (b) one-time plunge
to a much faster underlying decay.

If oscillation: eps_10 should rebound to substantial magnitude (next half-
  cycle peak). |eps_10/eps_9| > 10 likely.
If one-time plunge: eps_10 stays tiny. |eps_10/eps_9| < 1.
If level: |eps_10/eps_9| ~ 1, slow envelope at small magnitude.

K_10 has 2*3^9 = 39,366 states. Dense float64 = 12.4 GB (won't fit). Using
matrix-free left-multiply: precompute base[i] = (3i+1) mod N and the orbit
powers_inv2[v]; for each v compute targets vectorized and accumulate via
bincount. Chunked to keep working memory bounded.
"""
from __future__ import annotations

import csv
import gc
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
ENV_CSV = os.path.join(OUT_DIR, "result_q_sweep_test_1_envelope.csv")


def order_of_two(N):
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


class MatVecK:
    """Matrix-free Syracuse Markov K_k. Provides left-multiply x @ K
    (which equals K.T @ x). For power iteration on stationary measure."""

    def __init__(self, q, k, chunk=256):
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
        """Compute x_new = x @ K = K.T @ x (left action / stationary iter)."""
        n = self.n
        x_new = np.zeros(n, dtype=np.float64)
        for v0 in range(0, self.M, self.chunk):
            v1 = min(v0 + self.chunk, self.M)
            chunk_powers = self.powers_inv2[v0:v1]               # (c,)
            chunk_weights = self.weights[v0:v1]                  # (c,)
            # targets[v_offset, i] = (base[i] * chunk_powers[v_offset]) % N
            targets = (chunk_powers[:, None] * self.base[None, :]) % self.N
            j_targets = self.state_idx[targets]                  # (c, n)
            weighted = chunk_weights[:, None] * x[None, :]       # (c, n)
            x_new += np.bincount(
                j_targets.ravel(),
                weights=weighted.ravel(),
                minlength=n,
            )
        return x_new


def stationary_power_iteration_mfree(K_op, max_iter=10000, tol=1e-13,
                                      verbose_every=2):
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


def stationary_eigs_mfree(K_op):
    n = K_op.n
    op = LinearOperator(
        (n, n),
        matvec=lambda v: K_op.left_mul(v),
        dtype=np.float64,
    )
    vals, vecs = eigs(op, k=1, which="LM", maxiter=10000, tol=1e-12)
    val = float(vals[0].real)
    vec = vecs[:, 0].real
    s = vec.sum()
    if s < 0:
        vec = -vec
        s = -s
    vec = vec / s
    return val, vec


# Dense versions for k <= 9 (sanity)
def build_K_float_dense(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    pi = inv2
    for v in range(M):
        powers_inv2[v] = pi
        pi = (pi * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    Z_v = 1.0 - 2.0 ** (-M)
    weights = (2.0 ** (-np.arange(1, M + 1))) / Z_v
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (q * r + 1) % N
        targets = (base * powers_inv2) % N
        js = state_idx[targets]
        K[i_r] = np.bincount(js, weights=weights, minlength=n)
    return K, M


def stationary_power_dense(K, max_iter=10000, tol=1e-14):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    for it in range(max_iter):
        pi_new = pi @ K
        pi_new /= pi_new.sum()
        residual = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
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
    print("Computing eps_10 = S_10 - 7/15 (oscillation vs one-time plunge)")
    print("=" * 78)

    S_cache = load_cached_S()
    X = {0: 1.0}
    for k in [1, 2, 3, 4, 5]:
        X[k] = X[k - 1] + float(S_cache[k])
    print(f"\nX_5 (exact, cached) = {X[5]:.15f}")

    eps_all = {k: float(S_cache[k] - Fraction(7, 15)) for k in S_cache}

    # k=6..9 sanity recompute (dense, fast)
    for k in [6, 7, 8, 9]:
        print(f"\n--- k = {k} (sanity recompute) ---")
        t0 = time.time()
        K, M = build_K_float_dense(3, k)
        t_build = time.time() - t0
        print(f"  build K_{k} (states={K.shape[0]}, M={M}): "
              f"{t_build:.2f}s")
        t0 = time.time()
        pi, iters, res = stationary_power_dense(K, tol=1e-14)
        t_iter = time.time() - t0
        print(f"  power iter: {iters} iters, residual {res:.2e}, "
              f"{t_iter:.2f}s")
        sum_sq = float((pi ** 2).sum())
        X[k] = (3 ** k) * sum_sq
        S_k = X[k] - X[k - 1]
        eps_k = S_k - 7.0 / 15.0
        eps_all[k] = eps_k
        print(f"  X_{k} = {X[k]:.15f}, eps_{k} = {eps_k:+.10e}")
        del K, pi
        gc.collect()

    # k=10 (target, matrix-free)
    print(f"\n--- k = 10 (TARGET, matrix-free) ---")
    t0 = time.time()
    K_op = MatVecK(3, 10, chunk=512)
    t_build_10 = time.time() - t0
    n10 = K_op.n
    M10 = K_op.M
    print(f"  initialize MatVecK_10: {n10} states, M={M10}, chunk=512, "
          f"{t_build_10:.2f}s")

    # Sanity: row sums of K should be 1. Verify on a few rows.
    print("  K_10 row-sum sanity check (5 random rows)...")
    rng = np.random.default_rng(20260505)
    sample_rows = rng.choice(n10, size=5, replace=False)
    for ir in sample_rows:
        e_i = np.zeros(n10)
        e_i[ir] = 1.0
        # row i of K: K[i, :] = e_i^T @ K? No; e_i @ K = K[i, :].
        # Wait: e_i @ K[j] = sum_k e_i[k] K[k,j] = K[i, j]. So e_i @ K = K[i, :].
        row = K_op.left_mul(e_i)
        print(f"    row {ir}: sum = {row.sum():.15f}")

    t0 = time.time()
    pi10, iters10, res10 = stationary_power_iteration_mfree(
        K_op, tol=1e-13, max_iter=200, verbose_every=1
    )
    t_iter_10 = time.time() - t0
    print(f"  power iteration: {iters10} iters, residual {res10:.2e}, "
          f"{t_iter_10:.2f}s")
    print(f"  sum(pi_10) = {pi10.sum():.15f}  (expect 1)")

    print("\n  cross-check via scipy.sparse.linalg.eigs (Arnoldi)...")
    t0 = time.time()
    val_eigs, pi10_eigs = stationary_eigs_mfree(K_op)
    t_eigs = time.time() - t0
    print(f"  eigs: leading eigenvalue = {val_eigs:.12f}, {t_eigs:.2f}s")
    diff_l1 = float(np.linalg.norm(pi10 - pi10_eigs, ord=1))
    diff_max = float(np.max(np.abs(pi10 - pi10_eigs)))
    print(f"  |pi10_power - pi10_eigs|_1 = {diff_l1:.4e}, "
          f"max-abs = {diff_max:.4e}")

    sum_sq_10 = float((pi10 ** 2).sum())
    X[10] = (3 ** 10) * sum_sq_10
    S10 = X[10] - X[9]
    eps10 = S10 - 7.0 / 15.0
    eps_all[10] = eps10
    print(f"\n  X_10 = {X[10]:.15f}, S_10 = {S10:.15f}")
    print(f"  eps_10 = {eps10:+.10e}")

    sum_sq_10_eigs = float((pi10_eigs ** 2).sum())
    X10_eigs = (3 ** 10) * sum_sq_10_eigs
    S10_eigs = X10_eigs - X[9]
    eps10_eigs = S10_eigs - 7.0 / 15.0
    print(f"  cross-check eps_10 via eigs vector: {eps10_eigs:+.10e}")
    print(f"  agreement: |eps10_power - eps10_eigs| = "
          f"{abs(eps10 - eps10_eigs):.4e}")

    print("\n  FFT cross-check of S_10...")
    N10 = 3 ** 10
    pi_full = np.zeros(N10, dtype=np.float64)
    coprime_10 = np.array(
        [r for r in range(N10) if r % 3 != 0], dtype=np.int64
    )
    pi_full[coprime_10] = pi10
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N10)
    mask_nontrivial = xi_arr % 3 != 0
    S10_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps10_fft = S10_fft - 7.0 / 15.0
    print(f"  S_10 via FFT = {S10_fft:.15f}, eps_10_fft = "
          f"{eps10_fft:+.10e}")
    print(f"  agreement: |S_10_X - S_10_FFT| = {abs(S10 - S10_fft):.4e}")

    eps8 = eps_all[8]
    eps9 = eps_all[9]
    ratio_910 = abs(eps10 / eps9)
    ratio_89 = abs(eps9 / eps8)
    print(f"\n|eps_10 / eps_9| = {ratio_910:.6f}  "
          f"(compare |eps_9/eps_8| = {ratio_89:.6f})")
    print(f"sign(eps_10) = {'+' if eps10 > 0 else '-'}")

    # Decision
    if ratio_910 > 5.0:
        verdict = (
            f"OSCILLATION CONFIRMED — |eps_10/eps_9| = {ratio_910:.2f} >> 1. "
            "The k=9 collapse was a near-zero node, and eps_10 rebounds to a "
            "magnitude consistent with the next half-cycle peak. The k-space "
            "trajectory of eps_k carries a non-trivial oscillating component "
            "(complex eigenpair in K)."
        )
    elif 1.5 < ratio_910 <= 5.0:
        verdict = (
            f"PARTIAL REBOUND — |eps_10/eps_9| = {ratio_910:.2f} above 1.5 "
            "but well below an idealized half-cycle peak. Consistent with "
            "damped oscillation: amplitude decays between cycles, k=10 is "
            "post-node but with reduced peak."
        )
    elif 0.5 < ratio_910 <= 1.5:
        verdict = (
            f"LEVEL-OUT — |eps_10/eps_9| ~ 1 ({ratio_910:.3f}). Neither "
            "rebound nor continued plunge. The k=9 collapse may have been "
            "the leading edge of a slow envelope at small magnitude; "
            "structure is non-oscillating but at lower amplitude than "
            "expected from pre-spike trajectory."
        )
    else:
        verdict = (
            f"CONTINUED PLUNGE — |eps_10/eps_9| = {ratio_910:.4f} < 0.5. "
            "The k=9 collapse was not a node; underlying decay is genuinely "
            "much faster than pre-spike (rate < 0.5). Asymptotic rate "
            "implied is below 1/2."
        )

    print(f"\n*** VERDICT: {verdict}")

    print("\nUpdated ratio trajectory |eps_{k+1}/eps_k|:")
    for k in range(1, 10):
        rr = abs(eps_all[k + 1] / eps_all[k])
        print(f"  k = {k} -> {k+1}: {rr:.6f}")

    # Outputs
    out_md = os.path.join(OUT_DIR, "result_epsilon_10.md")
    md = []
    md.append("# Result: eps_10 = S_10 - 7/15")
    md.append("")
    md.append(f"**Date:** 2026-05-05.  Float64 matrix-free power iteration "
              f"on K_10 ({n10} states, M={M10}). Dense storage (12.4 GB) "
              f"avoided via chunked bincount per orbit-power. scipy.eigs "
              f"and FFT cross-checks.")
    md.append("")
    md.append(f"## Verdict")
    md.append("")
    md.append(verdict)
    md.append("")
    md.append("## Headline numbers")
    md.append("")
    md.append(f"- eps_10 = `{eps10:+.10e}` (power iter)")
    md.append(f"- eps_10 = `{eps10_eigs:+.10e}` (eigs cross-check; "
              f"agreement {abs(eps10 - eps10_eigs):.2e})")
    md.append(f"- eps_10 = `{eps10_fft:+.10e}` (FFT cross-check; agreement "
              f"{abs(eps10 - eps10_fft):.2e})")
    md.append(f"- |eps_10/eps_9| = **{ratio_910:.6f}**")
    md.append(f"- |eps_9/eps_8| = {ratio_89:.6f} (prior)")
    md.append("")
    md.append("## Ratio trajectory (k=1..10)")
    md.append("")
    md.append("| k → k+1 | |eps_{k+1}/eps_k| |")
    md.append("|---|---|")
    for k in range(1, 10):
        rr = abs(eps_all[k + 1] / eps_all[k])
        md.append(f"| {k} → {k+1} | {rr:.6f} |")
    md.append("")
    md.append("## eps_k table (k=1..10)")
    md.append("")
    md.append("| k | eps_k | source |")
    md.append("|---|---|---|")
    for k in [1, 2, 3, 4, 5]:
        md.append(f"| {k} | {float(S_cache[k] - Fraction(7,15)):+.10e} | "
                  f"exact rational (cached) |")
    for k in [6, 7, 8, 9]:
        md.append(f"| {k} | {eps_all[k]:+.10e} | float64 dense power iter |")
    md.append(f"| 10 | {eps10:+.10e} | float64 matrix-free + eigs + FFT |")
    md.append("")
    md.append("## Computation diagnostics (k=10)")
    md.append("")
    md.append(f"- States: {n10:,}, M = {M10:,}")
    md.append(f"- Init MatVecK: {t_build_10:.2f}s")
    md.append(f"- Power iter: {iters10} iterations, residual = {res10:.2e}, "
              f"{t_iter_10:.2f}s total")
    md.append(f"- eigs Arnoldi cross-check: {t_eigs:.2f}s, leading eval = "
              f"{val_eigs:.12f}")
    md.append(f"- |pi10_power - pi10_eigs|_1 = {diff_l1:.2e}, max-abs = "
              f"{diff_max:.2e}")
    md.append(f"- S_10 via X_10 - X_9 vs S_10 via FFT: diff = "
              f"{abs(S10 - S10_fft):.2e}")
    md.append(f"- pi_10 sum = {pi10.sum():.15f}")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_epsilon_10.py` — script (with matrix-free K_10)")
    md.append("- `result_epsilon_10.csv` — eps_k and ratios for k=1..10")
    md.append("- `result_epsilon_10.md` — this writeup")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"\nsaved {out_md}")

    out_csv = os.path.join(OUT_DIR, "result_epsilon_10.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "eps_k", "abs_eps_k", "ratio_to_prev_signed",
                    "abs_ratio_to_prev"])
        prev = None
        for k in range(1, 11):
            ek = eps_all[k]
            ratio_signed = "" if prev is None else f"{ek/prev:+.10f}"
            ratio_abs = "" if prev is None else f"{abs(ek/prev):.10f}"
            w.writerow([k, f"{ek:+.15e}", f"{abs(ek):.15e}",
                        ratio_signed, ratio_abs])
            prev = ek
    print(f"saved {out_csv}")
    print("\nDone.")


if __name__ == "__main__":
    main()
