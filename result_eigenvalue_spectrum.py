"""
result_eigenvalue_spectrum.py
=============================
Top-10 eigenvalue spectrum of K_k for k = 5, 6, 7. Pure data acquisition,
no fitting.

Reuses the float Markov kernel construction from result_77_4_K_spectrum.py.
For k=5 (n=162) and k=6 (n=486): dense np.linalg.eig.
For k=7 (n=1458): sparse scipy.sparse.linalg.eigs(k=10, which='LM') first,
fall back to dense if sparse solver fails to converge.

Reports for each k:
  - lambda_1, ..., lambda_10 (magnitude and argument)
  - spectral gap = 1 - |lambda_2|
  - |lambda_2| - |lambda_3| (multi-mode candidate flag)
  - |lambda_2| vs 1/2
  - Perron-Frobenius sanity (lambda_1 = 1, eigenvector real positive)
  - sparsity statistics (n_nonzero per row, max/mean)
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz"
OUT_MD = os.path.join(OUTDIR, "result_eigenvalue_spectrum.md")
OUT_CSV = os.path.join(OUTDIR, "result_eigenvalue_spectrum.csv")
OUT_DIAG = os.path.join(OUTDIR, "result_eigenvalue_spectrum_diagnostic.md")


def build_markov_float(k: int):
    """Right-stochastic K_k as float64 dense matrix on coprime states.

    For large M (= 2*3^{k-1}) probabilities 2^{-r_v} underflow float64 once
    r_v > 1074. Z = 1 - 2^{-M} -> 1 for any M >= 60 (relative error < 2^-60).
    """
    N = 3 ** k
    M = 2 * 3 ** (k - 1)
    M_eff = min(M, 1074)  # below this, 2^{-r_v} underflows in float64
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M_eff + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = np.zeros((n, n), dtype=np.float64)
    # Z = 1 - 2^{-M}; for M >= 60 this is 1.0 in float64 to relative err < 2^-60.
    Z_v = 1.0 if M >= 60 else (2 ** M - 1) / (2 ** M)
    p = 0.5 / Z_v  # p for r_v = 1
    for r in coprime_states:
        i = state_idx[r]
        p_cur = p
        for r_v in range(1, M_eff + 1):
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[i, state_idx[target]] += p_cur
            p_cur *= 0.5
    return K, coprime_states


def sparsity_stats(K: np.ndarray, tol: float = 1e-15):
    nz_per_row = np.sum(np.abs(K) > tol, axis=1)
    return {
        "n_rows": K.shape[0],
        "n_total": int(K.size),
        "n_nonzero": int(np.sum(nz_per_row)),
        "density": float(np.sum(nz_per_row) / K.size),
        "nz_per_row_min": int(nz_per_row.min()),
        "nz_per_row_max": int(nz_per_row.max()),
        "nz_per_row_mean": float(nz_per_row.mean()),
    }


def perron_check(K: np.ndarray, lam: complex, vec: np.ndarray):
    """Check lambda_1 ≈ 1 and corresponding eigenvector is real-positive (after sign norm)."""
    err_one = abs(lam - 1.0)
    v = vec.copy()
    if abs(v.real).sum() > abs(v.imag).sum():
        v = v.real
    else:
        v = v.imag
    if v.sum() < 0:
        v = -v
    v_normed = v / np.abs(v).max()
    n_neg = int(np.sum(v_normed < -1e-8))
    return err_one, n_neg, float(v_normed.min()), float(v_normed.max())


def top_eigs_dense(K: np.ndarray, n_top: int = 10):
    """Full dense eigendecomp; return top n_top eigenpairs by |lambda|."""
    t0 = time.time()
    lam, V = np.linalg.eig(K)
    elapsed = time.time() - t0
    idx = np.argsort(-np.abs(lam))
    lam = lam[idx]
    V = V[:, idx]
    return lam[:n_top], V[:, :n_top], elapsed, "dense (np.linalg.eig)"


def top_eigs_sparse(K: np.ndarray, n_top: int = 10):
    """Sparse Arnoldi for top n_top by |lambda|. Falls back to dense on failure."""
    t0 = time.time()
    Ksp = sp.csr_matrix(K)
    try:
        lam, V = spla.eigs(Ksp, k=n_top, which="LM", maxiter=5000, tol=1e-10)
        elapsed = time.time() - t0
        idx = np.argsort(-np.abs(lam))
        lam = lam[idx]
        V = V[:, idx]
        return lam, V, elapsed, "sparse (scipy.sparse.linalg.eigs, LM)"
    except spla.ArpackNoConvergence as e:
        print(f"    [sparse failed: {e}; falling back to dense]")
        return top_eigs_dense(K, n_top)


def main():
    print("=" * 78)
    print("Eigenvalue spectrum of K_k for k = 5, 6, 7")
    print("=" * 78)
    print()

    csv_rows = []
    md_blocks = []
    diag_blocks = []
    summaries = {}

    for k in [5, 6, 7]:
        print(f"--- k = {k} ---")
        t_build = time.time()
        K, coprime = build_markov_float(k)
        t_build = time.time() - t_build
        n = K.shape[0]
        spar = sparsity_stats(K)
        print(f"  build: {t_build:.2f}s, n = {n}, density = {spar['density']:.4f}")
        print(f"  nz/row: min={spar['nz_per_row_min']}, "
              f"max={spar['nz_per_row_max']}, mean={spar['nz_per_row_mean']:.1f}")

        if k == 5:
            lam, V, t_eig, solver = top_eigs_dense(K, 10)
        else:
            lam, V, t_eig, solver = top_eigs_sparse(K, 10)
        print(f"  spectrum: {t_eig:.2f}s ({solver})")

        err_one, n_neg, vmin, vmax = perron_check(K, lam[0], V[:, 0])
        pf_ok = err_one < 1e-8 and n_neg == 0
        print(f"  Perron-Frobenius: |lambda_1 - 1| = {err_one:.2e}, "
              f"n_neg_components = {n_neg}, v_range = [{vmin:.4f}, {vmax:.4f}]  "
              f"=> {'OK' if pf_ok else 'FAIL'}")

        print(f"  top 10 eigenvalues (sorted by |lambda|):")
        print(f"    {'i':>3}  {'Re(lambda)':>14}  {'Im(lambda)':>14}  "
              f"{'|lambda|':>12}  {'arg(lambda)':>14}")
        for i, l in enumerate(lam):
            arg = math.atan2(l.imag, l.real)
            print(f"    {i+1:>3}  {l.real:+14.10f}  {l.imag:+14.10f}  "
                  f"{abs(l):>12.10f}  {arg:+14.10f}")
            csv_rows.append({
                "k": k,
                "rank": i + 1,
                "lambda_real": float(l.real),
                "lambda_imag": float(l.imag),
                "abs_lambda": float(abs(l)),
                "arg_lambda": float(arg),
            })

        gap = 1.0 - abs(lam[1])
        l2 = abs(lam[1])
        l3 = abs(lam[2])
        diff_23 = abs(l2 - l3)
        diff_l2_half = abs(l2 - 0.5)
        print()
        print(f"  spectral gap        = 1 - |lambda_2| = {gap:.6f}")
        print(f"  |lambda_2|          = {l2:.6f}")
        print(f"  |lambda_3|          = {l3:.6f}")
        print(f"  |lambda_2|-|lambda_3| = {l2 - l3:.6f}")
        print(f"  |lambda_2| - 1/2    = {l2 - 0.5:+.6f}")
        print()

        summaries[k] = {
            "n": n,
            "lam": lam.copy(),
            "gap": gap,
            "l2": l2,
            "l3": l3,
            "diff_23": diff_23,
            "diff_l2_half": diff_l2_half,
            "pf_ok": pf_ok,
            "pf_err_one": err_one,
            "pf_n_neg": n_neg,
            "build_s": t_build,
            "eig_s": t_eig,
            "solver": solver,
            "sparsity": spar,
        }

    # ---------- write CSV ----------
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "k", "rank", "lambda_real", "lambda_imag", "abs_lambda", "arg_lambda"
        ])
        w.writeheader()
        w.writerows(csv_rows)
    print(f"[csv: {OUT_CSV}]")

    # ---------- write main MD ----------
    with open(OUT_MD, "w", encoding="utf-8") as fh:
        fh.write("# Eigenvalue spectrum of K_k (k = 5, 6, 7)\n\n")
        fh.write("Top-10 eigenvalues by |lambda| of the Tao-Syracuse Markov kernel K_k "
                 "on coprime states of Z/3^k. Pure data acquisition, no fitting.\n\n")
        fh.write("## Summary table\n\n")
        fh.write("| k  | n    | |lambda_2| | |lambda_3| | gap = 1-|lambda_2| | "
                 "|lambda_2|-1/2 | |lambda_2|-|lambda_3| | PF OK |\n")
        fh.write("|----|------|------------|------------|--------------------|"
                 "------------------|-----------------------|-------|\n")
        for k in [5, 6, 7]:
            s = summaries[k]
            fh.write(f"| {k} | {s['n']} | {s['l2']:.8f} | {s['l3']:.8f} | "
                     f"{s['gap']:.8f} | {s['l2']-0.5:+.8f} | {s['l2']-s['l3']:+.8f} | "
                     f"{'yes' if s['pf_ok'] else 'NO'} |\n")
        fh.write("\n")

        for k in [5, 6, 7]:
            s = summaries[k]
            fh.write(f"## k = {k} (n = {s['n']})\n\n")
            fh.write("```\n")
            fh.write(f"{'i':>3}  {'Re(lambda)':>16}  {'Im(lambda)':>16}  "
                     f"{'|lambda|':>14}  {'arg(lambda)':>14}\n")
            for i, l in enumerate(s["lam"]):
                arg = math.atan2(l.imag, l.real)
                fh.write(f"{i+1:>3}  {l.real:+16.12f}  {l.imag:+16.12f}  "
                         f"{abs(l):>14.12f}  {arg:+14.10f}\n")
            fh.write("```\n\n")

        fh.write("## Decision rubric (from brief)\n\n")
        fh.write("- |lambda_2| ~ 0.5 with clean gap to |lambda_3| < 0.4: "
                 "single-mode rate-1/2.\n")
        fh.write("- |lambda_2| ~ 0.5 but |lambda_3| also ~ 0.5: multi-mode, "
                 "comparable rates.\n")
        fh.write("- |lambda_2| significantly off 1/2: rate-1/2 conjecture wrong "
                 "at algebraic level.\n")
        fh.write("- Spectrum character changes with k: K_k structure is "
                 "k-dependent.\n\n")

        fh.write("## Empirical verdict\n\n")
        for k in [5, 6, 7]:
            s = summaries[k]
            fh.write(f"- **k={k}**: |lambda_2| = {s['l2']:.6f} "
                     f"(diff to 1/2 = {s['l2']-0.5:+.6f}), "
                     f"|lambda_3| = {s['l3']:.6f} "
                     f"(gap |lambda_2|-|lambda_3| = {s['l2']-s['l3']:+.6f}).\n")
        fh.write("\n")

        # tag verdict by simple rule from numbers
        l2s = [summaries[k]["l2"] for k in [5, 6, 7]]
        l3s = [summaries[k]["l3"] for k in [5, 6, 7]]
        any_l2_off_half = any(abs(x - 0.5) > 0.05 for x in l2s)
        any_l3_close_to_l2 = any((l2s[i] - l3s[i]) < 0.05 for i in range(3))
        l2_drift = abs(l2s[2] - l2s[0])

        fh.write(f"Numeric checks:\n\n")
        fh.write(f"- max |lambda_2 - 1/2| over k=5..7: {max(abs(x-0.5) for x in l2s):.6f}\n")
        fh.write(f"- min (|lambda_2| - |lambda_3|) over k=5..7: "
                 f"{min(l2s[i]-l3s[i] for i in range(3)):.6f}\n")
        fh.write(f"- |lambda_2| drift across k=5..7: {l2_drift:.6f}\n")
        fh.write(f"- any |lambda_2| off 1/2 by >0.05?  {'yes' if any_l2_off_half else 'no'}\n")
        fh.write(f"- any |lambda_2|-|lambda_3| < 0.05? {'yes' if any_l3_close_to_l2 else 'no'}\n")
    print(f"[md:  {OUT_MD}]")

    # ---------- write diagnostic MD ----------
    with open(OUT_DIAG, "w", encoding="utf-8") as fh:
        fh.write("# Eigenvalue spectrum — diagnostic\n\n")
        fh.write("## Solver, sparsity, timing\n\n")
        fh.write("| k | n    | density | nz/row min | nz/row max | nz/row mean | "
                 "build (s) | eig (s) | solver |\n")
        fh.write("|---|------|---------|------------|------------|-------------|"
                 "-----------|---------|--------|\n")
        for k in [5, 6, 7]:
            s = summaries[k]
            sp_ = s["sparsity"]
            fh.write(f"| {k} | {s['n']} | {sp_['density']:.4f} | "
                     f"{sp_['nz_per_row_min']} | {sp_['nz_per_row_max']} | "
                     f"{sp_['nz_per_row_mean']:.1f} | "
                     f"{s['build_s']:.2f} | {s['eig_s']:.2f} | "
                     f"{s['solver']} |\n")
        fh.write("\n## Perron-Frobenius checks\n\n")
        fh.write("| k | |lambda_1 - 1| | n_neg_components (eigvec_1) | PF OK |\n")
        fh.write("|---|----------------|------------------------------|-------|\n")
        for k in [5, 6, 7]:
            s = summaries[k]
            fh.write(f"| {k} | {s['pf_err_one']:.2e} | {s['pf_n_neg']} | "
                     f"{'yes' if s['pf_ok'] else 'NO'} |\n")
        fh.write("\nKernel construction: float64 K_k from Tao-Syracuse build_markov_float "
                 "(coprime states r in Z/3^k with r mod 3 != 0; transition r -> "
                 "((3r+1)*2^{-v}) mod 3^k with prob (1/2^v)/Z_M).\n")
    print(f"[diag: {OUT_DIAG}]")


if __name__ == "__main__":
    main()
