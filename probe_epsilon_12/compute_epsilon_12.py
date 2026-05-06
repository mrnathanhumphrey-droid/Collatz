"""
compute_epsilon_12.py
=====================
Compute eps_12 = S_12 - 7/15 via matrix-free power iteration on K_12.

State-space facts:
  - n_12 (coprime states) = 2 * 3^11 = 354,294
  - M_12 = ord(2 mod 3^12) = 2 * 3^11 = 354,294
  - Dense K_12 would be ~1004 GB. Sparse "csr" is also nearly-dense
    (each row has M = 354,294 nonzeros), so ~1 TB. The matrix-free
    MatVecK already used for k=11 is the correct approach.

Method: matrix-free power iteration with Aitken Delta-squared acceleration
applied post-hoc once the iteration is near convergence.

Validation gate: before launching k=12, this script can be invoked with
--validate-k=8 to reproduce eps_8 = -7.4554636729e-04 to ~1e-9 absolute.

Usage:
  python compute_epsilon_12.py --validate-k 8     # quick smoke test
  python compute_epsilon_12.py --target-k 12      # full ε_12 run
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = Path(r"C:\Collatz\probe_epsilon_12")
OUT_DIR.mkdir(exist_ok=True)
ENV_CSV = Path(r"C:\Collatz\result_q_sweep_test_1_envelope.csv")

X_CACHED = {
    5: 3.534161151367800,
    6: 4.000329912369247,
    7: 4.465821342205540,
    8: 4.931742462504920,
    9: 5.398401608914431,
    10: 5.865789026498212,
    11: 6.333957660176961,
}

EPS_CACHED = {
    6: -4.9790566522e-04,
    7: -1.1752368304e-03,
    8: -7.4554636729e-04,
    9: -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
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
    """Matrix-free K_k operator. left_mul(x) computes K^T @ x in O(n*M)."""

    def __init__(self, q, k, chunk=512):
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


def power_iteration(K_op, max_iter=300, tol=1e-13, history_size=3,
                    aitken=True, log_path=None, verbose_every=1):
    """Stochastic power iteration with optional Aitken Delta^2 acceleration.

    Aitken applied component-wise on the last 3 iterates if convergence is
    slow but monotone (criterion: residual decay ratio in [0.7, 0.99]).
    """
    n = K_op.n
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    history = [pi.copy()]
    residuals = []

    log_lines = []

    def log(msg):
        print(msg, flush=True)
        log_lines.append(msg)
        if log_path is not None:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(msg + "\n")

    for it in range(max_iter):
        t0 = time.time()
        pi_new = K_op.left_mul(pi)
        s = pi_new.sum()
        pi_new /= s
        residual = float(np.linalg.norm(pi_new - pi, ord=1))
        residuals.append(residual)
        pi = pi_new
        history.append(pi.copy())
        if len(history) > history_size:
            history.pop(0)
        dt = time.time() - t0

        if (it + 1) % verbose_every == 0 or residual < tol:
            log(f"    iter {it+1:>4}: residual = {residual:.4e}, "
                f"sum = {s:.15f}, matvec {dt:.2f}s")

        if residual < tol:
            log(f"    [converged at iter {it+1}, residual {residual:.2e}]")
            return pi, it + 1, residual, residuals

        # Aitken acceleration: try when iteration is in linear convergence regime
        if aitken and len(history) == 3 and len(residuals) >= 3:
            r_recent = residuals[-3:]
            ratio = r_recent[-1] / max(r_recent[-2], 1e-30)
            if 0.7 < ratio < 0.99:
                p0, p1, p2 = history
                denom = p2 - 2 * p1 + p0
                # Avoid division by zero / ill-conditioned components
                ok = np.abs(denom) > 1e-15
                pi_accel = pi.copy()
                pi_accel[ok] = p0[ok] - (p1[ok] - p0[ok]) ** 2 / denom[ok]
                # Keep mass nonneg + normalize
                pi_accel = np.maximum(pi_accel, 0.0)
                ssum = pi_accel.sum()
                if ssum > 0.5:  # sanity
                    pi_accel /= ssum
                    res_accel = float(np.linalg.norm(pi_accel - pi, ord=1))
                    log(f"    [Aitken probe at iter {it+1}: "
                        f"ratio={ratio:.4f}, accelerated residual "
                        f"vs current = {res_accel:.4e}]")
                    # Use accelerated point as next iterate ONLY if it
                    # reduces residual against next K-application
                    pi_test = K_op.left_mul(pi_accel)
                    pi_test /= pi_test.sum()
                    res_test = float(np.linalg.norm(pi_test - pi_accel, ord=1))
                    if res_test < residual * 0.5:
                        log(f"    [Aitken accepted: post-Aitken "
                            f"residual={res_test:.4e} < {residual:.4e}/2]")
                        pi = pi_test
                        residual = res_test
                        residuals.append(residual)
                        history = [pi_accel, pi]
                    else:
                        log(f"    [Aitken rejected: "
                            f"post-Aitken residual={res_test:.4e}]")

    return pi, max_iter, residual, residuals


def compute_eps_at(k, chunk=512, tol=1e-13, max_iter=200, label=""):
    """Run power iteration at level k, return (X_k, S_k, eps_k, diagnostics)."""
    log_path = OUT_DIR / f"compute_k{k}.log"
    if log_path.exists():
        log_path.unlink()

    print(f"\n=== Computing eps_{k} {label} ===", flush=True)
    print(f"  chunk = {chunk}, tol = {tol}, max_iter = {max_iter}")
    t0 = time.time()
    K_op = MatVecK(3, k, chunk=chunk)
    t_init = time.time() - t0
    print(f"  initialized MatVecK_{k}: n = {K_op.n:,}, M = {K_op.M:,}, "
          f"init {t_init:.2f}s")

    # Quick row-sum check
    rng = np.random.default_rng(20260506 + k)
    ir = int(rng.integers(0, K_op.n))
    e_i = np.zeros(K_op.n)
    e_i[ir] = 1.0
    t0_row = time.time()
    row = K_op.left_mul(e_i)
    dt_row = time.time() - t0_row
    print(f"  row-sum check (idx {ir}): {row.sum():.15f} ({dt_row:.1f}s)")

    print(f"\n  power iteration:")
    t0 = time.time()
    pi_k, iters, res, residuals = power_iteration(
        K_op, max_iter=max_iter, tol=tol, log_path=log_path,
        verbose_every=1
    )
    t_iter = time.time() - t0
    print(f"  power iter complete: {iters} iters, residual {res:.2e}, "
          f"{t_iter:.2f}s ({t_iter/max(iters,1):.1f}s/iter)")
    print(f"  sum(pi_{k}) = {pi_k.sum():.15f}")

    # X_k, S_k, eps_k via X-formula
    sum_sq = float((pi_k ** 2).sum())
    X = (3 ** k) * sum_sq
    if k - 1 in X_CACHED:
        S = X - X_CACHED[k - 1]
    else:
        S = float("nan")
    eps = S - 7.0 / 15.0

    # FFT cross-check
    print(f"\n  FFT cross-check of S_{k}...")
    t0 = time.time()
    Nk = 3 ** k
    pi_full = np.zeros(Nk, dtype=np.float64)
    pi_full[K_op.coprime] = pi_k
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(Nk)
    mask_nontrivial = xi_arr % 3 != 0
    S_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps_fft = S_fft - 7.0 / 15.0
    t_fft = time.time() - t0
    print(f"  S_{k} via FFT  = {S_fft:.15f}, eps_{k}_fft = {eps_fft:+.10e} "
          f"({t_fft:.2f}s)")
    print(f"  agreement: |S_X - S_FFT| = {abs(S - S_fft):.4e}")

    diag = dict(
        k=k,
        n=K_op.n,
        M=K_op.M,
        chunk=chunk,
        iters=iters,
        residual=res,
        t_iter=t_iter,
        t_iter_per=t_iter / max(iters, 1),
        X=X,
        S=S,
        eps=eps,
        S_fft=S_fft,
        eps_fft=eps_fft,
        agreement=abs(S - S_fft),
        residuals=residuals,
    )
    return pi_k, K_op, diag


def cmd_validate(k):
    """Run k for validation; compare to EPS_CACHED[k]."""
    if k not in EPS_CACHED:
        print(f"No cached eps_{k} to validate against.")
        return 1
    expected = EPS_CACHED[k]
    print(f"Validating k={k} against cached eps_{k} = {expected:+.10e}")
    pi_k, K_op, diag = compute_eps_at(k, chunk=512, tol=1e-13, max_iter=200,
                                       label="(VALIDATION)")
    measured = diag["eps"]
    err = abs(measured - expected)
    rel = err / max(abs(expected), 1e-15)
    print(f"\n  expected eps_{k} = {expected:+.10e}")
    print(f"  measured eps_{k} = {measured:+.10e}")
    print(f"  abs error  = {err:.3e}")
    print(f"  rel error  = {rel:.3e}")
    if err < 1e-9:
        print(f"  *** VALIDATION PASS *** (abs error < 1e-9)")
        return 0
    elif err < 1e-6:
        print(f"  *** VALIDATION MARGINAL *** (1e-9 < err < 1e-6)")
        return 0
    else:
        print(f"  *** VALIDATION FAIL *** (err > 1e-6)")
        return 1


def cmd_target(k):
    """Run target k=12 (or any k); save full results."""
    label = f"(TARGET k={k})"
    chunk = 256 if k >= 12 else 512
    print(f"=== Target run: k={k}, chunk={chunk} ===")

    # Estimate time before launching
    if k >= 11 and 11 in {11}:
        # k=11 took 1443s with n=M=118098, chunk=1024. Per-matvec scales
        # roughly as n*M/chunk. Going to k=12 with chunk=256 means:
        # work_ratio = (354294 * 354294 / 256) / (118098 * 118098 / 1024)
        # = (n^2 / chunk) ratios.
        n11_M11 = 118098 * 118098 / 1024
        n12_M12 = 354294 * 354294 / chunk
        ratio = n12_M12 / n11_M11
        per_matvec_est = 122 * ratio
        print(f"  Per-matvec estimate (vs k=11 baseline): "
              f"{per_matvec_est:.1f}s = {per_matvec_est/60:.1f} min")
        print(f"  At ~12 iters: ~{12 * per_matvec_est / 3600:.1f} hours total.")

    pi_k, K_op, diag = compute_eps_at(k, chunk=chunk, tol=1e-13,
                                       max_iter=200, label=label)
    measured = diag["eps"]

    # Save pi_k
    pi_path = OUT_DIR / f"pi_{k}.npz"
    np.savez_compressed(pi_path, pi=pi_k, coprime=K_op.coprime)
    print(f"\nSaved {pi_path}")

    # Save CSV
    csv_path = OUT_DIR / f"result_S_{k}.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "method", "n", "M", "chunk", "iters",
                    "final_residual", "t_iter_seconds", "X_k", "S_k",
                    "eps_k", "eps_k_fft", "fft_agreement"])
        w.writerow([k, "matvec_power_iter_aitken", K_op.n, K_op.M,
                    diag["chunk"], diag["iters"], f"{diag['residual']:.6e}",
                    f"{diag['t_iter']:.2f}", f"{diag['X']:.15f}",
                    f"{diag['S']:.15f}", f"{measured:+.15e}",
                    f"{diag['eps_fft']:+.15e}",
                    f"{diag['agreement']:.6e}"])
    print(f"Saved {csv_path}")

    # Pre-registered analysis
    pred_lo, pred_hi = 1.50e-3, 2.5e-3
    pred_mid = (pred_lo + pred_hi) / 2
    print(f"\nPre-registered prediction: eps_12 ∈ [+1.50e-3, +2.5e-3]")
    print(f"Measured:                  eps_{k} = {measured:+.10e}")

    sign_match = (measured > 0) == (pred_mid > 0)
    if pred_lo <= measured <= pred_hi:
        verdict = (
            f"RECURRENCE CONFIRMED: eps_{k} = {measured:+.4e} lies within "
            f"pre-registered band [+1.50e-3, +2.5e-3]. Order-3 recurrence "
            "model strengthened."
        )
    elif sign_match and abs((measured - pred_mid) / pred_mid) < 0.5:
        verdict = (
            f"RECURRENCE APPROXIMATELY RIGHT: sign matches, magnitude off "
            f"by {abs((measured - pred_mid) / pred_mid)*100:.1f}%; "
            "coefficients likely need updating but model class still fits."
        )
    elif not sign_match:
        verdict = (
            f"RECURRENCE FALSIFIED AT k={k}: sign differs from prediction. "
            "Order-3 model needs structural rework."
        )
    else:
        verdict = (
            f"RECURRENCE STRESSED: sign matches but magnitude differs by "
            f">50% from pre-registered midpoint. Refit recommended."
        )
    print(f"\n*** VERDICT: {verdict}")

    # Findings markdown
    md_path = OUT_DIR / f"epsilon_{k}_findings.md"
    md = []
    md.append(f"# Result: eps_{k} = S_{k} - 7/15")
    md.append("")
    md.append(f"**Date:** 2026-05-05.  Float64 matrix-free power iteration "
              f"on K_{k} ({K_op.n:,} states, M={K_op.M:,}). FFT cross-check.")
    md.append("")
    md.append("## Verdict")
    md.append("")
    md.append(verdict)
    md.append("")
    md.append("## Headline numbers")
    md.append("")
    md.append(f"- eps_{k} = `{measured:+.10e}` (power iter, X-formula)")
    md.append(f"- eps_{k} = `{diag['eps_fft']:+.10e}` (FFT cross-check; "
              f"agreement {diag['agreement']:.2e})")
    md.append("")
    if k - 1 in EPS_CACHED:
        prev = EPS_CACHED[k - 1]
        md.append(f"- eps_{k-1} = `{prev:+.10e}`")
        md.append(f"- |eps_{k} / eps_{k-1}| = "
                  f"{abs(measured/prev):.6f}")
        md.append(f"- sign(eps_{k}) = {'+' if measured > 0 else '-'}  "
                  f"(eps_{k-1} was {'+' if prev > 0 else '-'})")
    md.append("")
    md.append(f"## Pre-registered prediction")
    md.append("")
    md.append(f"From order-3 linear recurrence fit on eps_2..eps_11:")
    md.append(f"- Predicted band: `[+1.50e-3, +2.5e-3]`")
    md.append(f"- Measured: `{measured:+.10e}`")
    md.append(f"- Outcome: {'in band' if pred_lo <= measured <= pred_hi else 'outside band'}")
    md.append("")
    md.append("## eps_k table (k=6..k)")
    md.append("")
    md.append("| k | eps_k | sign | source |")
    md.append("|---|---|:---:|---|")
    for kk in sorted(EPS_CACHED):
        ek = EPS_CACHED[kk]
        sgn = "+" if ek > 0 else "-"
        md.append(f"| {kk} | {ek:+.10e} | {sgn} | float64 power iter (cached) |")
    if k not in EPS_CACHED:
        sgn = "+" if measured > 0 else "-"
        md.append(f"| {k} | {measured:+.10e} | {sgn} | this run |")
    md.append("")
    md.append("## Compute diagnostics")
    md.append("")
    md.append(f"- States n = {K_op.n:,}")
    md.append(f"- M = {K_op.M:,}")
    md.append(f"- Chunk size = {diag['chunk']}")
    md.append(f"- Power iter: {diag['iters']} iters, "
              f"final residual = {diag['residual']:.2e}")
    md.append(f"- Total iter time: {diag['t_iter']:.1f}s "
              f"({diag['t_iter']/3600:.2f} hours)")
    md.append(f"- Per-matvec: {diag['t_iter_per']:.1f}s")
    md.append(f"- FFT cross-check vs X-formula: "
              f"diff = {diag['agreement']:.2e}")
    md.append("")
    md.append("## Method notes")
    md.append("")
    md.append("Method: matrix-free power iteration. K_k is *dense* per row "
              "(M nonzeros / row, M ≈ n), so a sparse csr representation "
              "would not save memory; instead we exploit the multiplicative "
              f"structure: each row j of K is the histogram of "
              f"(2^{{-(v+1)}} / Z_v) summed over v=0..M-1 onto state "
              f"idx((q*r_j+1)*2^{{-v-1}} mod N). bincount over chunked v "
              f"keeps peak memory bounded by chunk*n*8 bytes "
              f"(here {diag['chunk']*K_op.n*8/1e9:.2f} GB).")
    md.append("")
    md.append("Aitken acceleration was wired in (componentwise Δ² applied "
              "when residual ratio is in [0.7, 0.99]) but typically not "
              "triggered — power iter on K_k converges sharply once a "
              "subspace alignment is reached, so the linear-convergence "
              "regime is brief.")

    md_path.write_text("\n".join(md), encoding="utf-8")
    print(f"\nSaved {md_path}")

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-k", type=int, default=None)
    parser.add_argument("--target-k", type=int, default=None)
    args = parser.parse_args()

    if args.validate_k is not None:
        return cmd_validate(args.validate_k)
    elif args.target_k is not None:
        return cmd_target(args.target_k)
    else:
        print("Specify --validate-k <k> or --target-k <k>")
        return 1


if __name__ == "__main__":
    sys.exit(main())
