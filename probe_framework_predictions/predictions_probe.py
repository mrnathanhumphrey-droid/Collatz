"""
Probe B: quantitative-prediction tests of the AS / DG / DS-C framework on K_k.

Companion to probe_ayyer_singla_test/ ("Probe A"), which concluded Outcome C
(predicted B_k spectrum lives on circle |λ − 2/3| = 1/3, |λ| in [1/3, 1];
K_k's actual non-trivial spectrum at |λ| ~ 10^(-3); +1 affine shift breaks
group-walk structure).

Probe B tests three quantitative predictions:
  Phase 1 — per-character iteration: |Q̂_n(χ)|^2 ~ |F(χ)|^(2n).
            Compare extracted |F(χ)| to K_k's measured |λ| at the same rank.
            <5% match = confirm; >20% = decisive failure.
  Phase 2 — Jordan-block: rank((K-λI)^j) drops linearly in j up to block size m,
            then plateaus. Identify block sizes per top eigenvalue.
  Phase 3 — DS-C bound: t_mix ≤ C·γ_k^2·log(1/ε). Compare to empirical t_mix
            of K_k. Cayley walk natural for B_k; bound likely misaligned with K_k.

Run k = 5, 6, 7. Outputs:
  per_character_rates_k{k}.csv
  jordan_blocks_k{k}.csv
  mixing_time_k{k}.csv
"""

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import eigs as sparse_eigs

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path("C:/Collatz/probe_framework_predictions")
OUTDIR.mkdir(exist_ok=True)


def order_of_2_mod(N):
    o, p = 1, 2 % N
    while p != 1:
        p = (p * 2) % N
        o += 1
    return o


def build_K(k):
    """Standard K_k builder (full Syracuse with +1)."""
    N = 3 ** k
    M = order_of_2_mod(N)
    coprime = [r for r in range(N) if r % 3 != 0]
    n = len(coprime)
    idx = {r: i for i, r in enumerate(coprime)}
    inv2 = pow(2, -1, N)
    inv2_pow = [1, inv2]
    for v in range(2, M + 1):
        inv2_pow.append((inv2_pow[-1] * inv2) % N)
    Z = 1.0 - 2.0 ** (-M)
    weights = [None] + [(2.0 ** (-v)) / Z for v in range(1, M + 1)]
    K = np.zeros((n, n))
    for i, r in enumerate(coprime):
        base_mod = (3 * r + 1) % N
        for v in range(1, M + 1):
            target = (base_mod * inv2_pow[v]) % N
            j = idx[target]
            K[i, j] += weights[v]
    return K, coprime, idx, M


def build_B(k):
    """Multiplicative-only chain B_k: x -> x * 2^(-V) on (Z/3^k)*."""
    N = 3 ** k
    M = order_of_2_mod(N)
    coprime = [r for r in range(N) if r % 3 != 0]
    n = len(coprime)
    idx = {r: i for i, r in enumerate(coprime)}
    inv2 = pow(2, -1, N)
    inv2_pow = [1, inv2]
    for v in range(2, M + 1):
        inv2_pow.append((inv2_pow[-1] * inv2) % N)
    Z = 1.0 - 2.0 ** (-M)
    weights = [None] + [(2.0 ** (-v)) / Z for v in range(1, M + 1)]
    B = np.zeros((n, n))
    for i, r in enumerate(coprime):
        for v in range(1, M + 1):
            target = (r * inv2_pow[v]) % N
            j = idx[target]
            B[i, j] += weights[v]
    return B, coprime, idx, M


def power_iter(K, max_iter=50000, tol=1e-15):
    n = K.shape[0]
    pi = np.ones(n) / n
    KT = K.T.copy()
    for it in range(max_iter):
        pi_new = KT @ pi
        s = pi_new.sum()
        pi_new /= s
        diff = float(np.linalg.norm(pi_new - pi, ord=np.inf))
        if diff < tol:
            return pi_new, it + 1, diff
        pi = pi_new
    return pi, max_iter, diff


# =========================================================================
# Phase 1: per-character iteration
# =========================================================================

def discrete_logs(coprime, k):
    """For (Z/3^k)* cyclic, discrete log of each coprime under generator g.
    Returns dict r -> a such that r = g^a mod 3^k.
    Generator: -1 has order 2, so we use the generator of the cyclic group.
    For (Z/3^k)*, 2 has order 2*3^(k-1) = full order, so 2 generates.
    """
    N = 3 ** k
    n_k = 2 * 3 ** (k - 1)
    powers = {1: 0}
    g = 2
    cur = 1
    for a in range(1, n_k):
        cur = (cur * g) % N
        powers[cur] = a
    # All coprime residues should be in powers
    return powers


def phase1_per_character(K, pi, coprime, idx, k, n_steps=20):
    """Iterate K_k from delta_{r_0=1}. Compute Fourier coefficients at each
    character chi_j (j=0..n_k-1). Track |Q_hat_n(chi_j)|^2 over n.

    For each character, fit log|Q_hat_n|^2 vs n -> slope = 2 log|F(chi)|.
    Compare extracted |F(chi)| to top eigenvalues |lambda| of K.
    """
    N = 3 ** k
    n_k = len(coprime)
    n = K.shape[0]

    # Discrete log table: r -> a such that r = 2^a mod 3^k
    dlog = discrete_logs(coprime, k)

    # Centered measure: Q_n - pi (so mean drops out, characters track decay)
    KT = K.T.copy()

    # Use uniform initial - pi (mean-zero starting deviation)
    # Or use delta_1 - pi
    delta = np.zeros(n)
    delta[idx[1]] = 1.0
    Qn = delta - pi  # initial deviation

    # Precompute character matrix: chi[j, r_idx] = exp(2 pi i j a_r / n_k)
    a = np.array([dlog[r] for r in coprime])
    omega = np.exp(-2j * np.pi / n_k)
    # We'll just compute selected chi values. n_k can be up to 1458 — need ~20 chi.
    # Selected: j = 1..n_k-1, but want top-20 by initial magnitude
    # First do FFT-style: chi[j, i] = omega^(j * a[i])
    # Compute Q_hat_n[j] = sum_i chi_j(coprime[i]) Qn[i] for n=0..n_steps
    # We'll select top-N j's after first iteration

    # Iterate, recording all chi
    Qn_hist = [Qn.copy()]
    for _ in range(n_steps):
        Qn = KT @ Qn
        Qn_hist.append(Qn.copy())

    # Compute character coefficients efficiently via DFT-like weighting
    Qn_hat = np.zeros((n_steps + 1, n_k), dtype=complex)
    for n in range(n_steps + 1):
        # Q_hat_n[j] = sum_i exp(-2 pi i j a[i] / n_k) * Qn_hist[n][i]
        # For each j, dot product. Use vectorized.
        for j in range(n_k):
            chi_vals = omega ** (j * a)
            Qn_hat[n, j] = np.dot(chi_vals, Qn_hist[n])

    # |Q_hat_n[j]|^2
    Qn_hat_sq = np.abs(Qn_hat) ** 2

    # Fit log|Q_hat_n|^2 vs n (skip n=0 which is initial)
    # Slope = 2 log|F(chi)|
    # |F(chi)| = exp(slope/2)
    F_extracted = np.zeros(n_k)
    for j in range(n_k):
        ys = Qn_hat_sq[1:, j]
        # Use only points where |Q_hat_n|^2 > floor (else underflow / float noise)
        mask = ys > 1e-30
        if mask.sum() < 3:
            F_extracted[j] = float('nan')
            continue
        valid_ys = np.log(ys[mask])
        valid_xs = np.arange(1, n_steps + 1)[mask]
        slope, intercept = np.polyfit(valid_xs, valid_ys, 1)
        F_extracted[j] = np.exp(slope / 2)

    return Qn_hat, Qn_hat_sq, F_extracted


# =========================================================================
# Phase 2: Jordan-block structure
# =========================================================================

def phase2_jordan(K, eigvals, k, top_N=20, rtol=1e-10):
    """For each top-N eigenvalue λ, compute rank((K - λI)^j) for j=1..k.
    If rank drops linearly until j = m then plateaus at n - m, then Jordan
    block of size m exists at λ.

    Returns: list of dicts with rank_drops per λ, inferred block_size.
    """
    n = K.shape[0]
    # Sort eigvals by |λ| descending; pick top-N unique magnitudes (cluster within rtol)
    sorted_idx = np.argsort(-np.abs(eigvals))
    selected = []
    seen_mags = []
    for ii in sorted_idx:
        if len(selected) >= top_N:
            break
        lam = eigvals[ii]
        mag = abs(lam)
        if mag < 1e-12:
            continue
        # cluster: skip if mag close to a seen mag (within 1e-8 abs)
        is_new = all(abs(mag - m) > 1e-8 * max(mag, 1e-10) for m in seen_mags)
        if is_new:
            selected.append(lam)
            seen_mags.append(mag)

    results = []
    for lam in selected:
        n_minus_lam_I = K - lam * np.eye(n, dtype=K.dtype)
        # Computing (K - λI)^j: matrix multiplication. For large j, condition explodes.
        # Use SVD-based rank.
        ranks = []
        Mj = n_minus_lam_I.copy()
        for j in range(1, k + 1):
            # rank of Mj
            sv = np.linalg.svd(Mj, compute_uv=False)
            tol = max(Mj.shape) * sv.max() * rtol if sv.size > 0 else 0.0
            r_j = int((sv > tol).sum())
            ranks.append(r_j)
            if j < k:
                Mj = Mj @ n_minus_lam_I

        # null spaces grow: nullity[j] = n - rank[j]
        nullities = [n - r for r in ranks]
        # block sizes: difference of nullities until plateau
        # nullity[j] - nullity[j-1] = number of blocks of size >= j
        block_size_at_least = [nullities[0]]
        for j in range(1, len(nullities)):
            block_size_at_least.append(nullities[j] - nullities[j - 1])
        # block size m exists if nullity[m] > nullity[m-1] AND nullity[m+1] = nullity[m]
        # Estimated max block size: largest j where nullity[j] > nullity[j-1]
        max_block_size = 1
        for j in range(1, len(nullities)):
            if nullities[j] > nullities[j - 1] + 0:  # stricter would use rtol
                max_block_size = j + 1

        # Geometric multiplicity (= nullity at j=1)
        geom_mult = nullities[0]
        # Algebraic multiplicity from Jordan: sum of all block sizes for this λ
        # = max nullity (= dim of generalized eigenspace)
        alg_mult = max(nullities)

        results.append({
            "lambda_re": lam.real,
            "lambda_im": lam.imag,
            "abs_lambda": abs(lam),
            "geom_mult": geom_mult,
            "alg_mult": alg_mult,
            "max_block_size": max_block_size,
            "nullities": nullities,
            "rank_drops": [n - r for r in ranks],
            "k": k,
        })

    return results


# =========================================================================
# Phase 3: DS-C mixing time vs Cayley diameter bound
# =========================================================================

def cayley_diameter_and_balls(k, gen_set):
    """Compute Cayley diameter of (Z/3^k)* under generating set gen_set.
    gen_set is a list of group elements (units mod 3^k) that generates with their inverses.
    Returns (diameter, ball_sizes) where ball_sizes[r] = |B(r)|.
    """
    N = 3 ** k
    coprime = [r for r in range(N) if r % 3 != 0]
    coprime_set = set(coprime)

    # Symmetrize: include inverses
    gen_with_inv = set()
    for g in gen_set:
        gen_with_inv.add(g % N)
        gen_with_inv.add(pow(g, -1, N))

    visited = {1: 0}
    frontier = [1]
    diameter = 0
    while frontier:
        next_frontier = []
        for x in frontier:
            for g in gen_with_inv:
                y = (x * g) % N
                if y not in visited and y in coprime_set:
                    visited[y] = visited[x] + 1
                    diameter = max(diameter, visited[y])
                    next_frontier.append(y)
        frontier = next_frontier

    ball_sizes = [0] * (diameter + 1)
    for x, d in visited.items():
        ball_sizes[d] += 1
    # cumulative
    cum_balls = []
    cum = 0
    for s in ball_sizes:
        cum += s
        cum_balls.append(cum)

    return diameter, cum_balls, len(visited)


def empirical_mixing_time(K, pi, idx, coprime, max_n=500, eps=0.25):
    """Smallest n such that ||K^n delta_r - pi||_TV < eps for r in sample of starts."""
    n_states = K.shape[0]
    KT = K.T.copy()
    sample_starts = coprime[: min(20, len(coprime))]  # sample first 20 starts
    times = []
    for r0 in sample_starts:
        Q = np.zeros(n_states)
        Q[idx[r0]] = 1.0
        t = -1
        for n in range(1, max_n + 1):
            Q = KT @ Q
            tv = 0.5 * np.sum(np.abs(Q - pi))
            if tv < eps:
                t = n
                break
        times.append(t)
    return times, sample_starts


# =========================================================================
# Main loop
# =========================================================================

def main():
    summary = []
    k_list = [5, 6, 7]

    for k in k_list:
        print(f"\n{'='*72}")
        print(f"k = {k}")
        print(f"{'='*72}")
        N = 3 ** k

        # Build K_k and B_k
        t0 = time.time()
        K, coprime, idx, M = build_K(k)
        n = K.shape[0]
        pi, iters, residual = power_iter(K)
        print(f"K_k built (n={n}, M={M}); pi iters={iters} residual={residual:.1e} ({time.time()-t0:.1f}s)")

        t0 = time.time()
        B, _, _, _ = build_B(k)
        pi_B, iters_B, _ = power_iter(B)
        print(f"B_k built; pi_B iters={iters_B} ({time.time()-t0:.1f}s)")

        # Compute K_k full eigvals via dense eig
        t0 = time.time()
        eigvals_K, eigvecs_K = la.eig(K)
        print(f"K_k full eig done ({time.time()-t0:.1f}s)")
        # Sort descending |λ|
        order = np.argsort(-np.abs(eigvals_K))
        eigvals_K_sorted = eigvals_K[order]

        # =========================================================
        # Phase 1: per-character iteration
        # =========================================================
        print("\n--- Phase 1: per-character iteration ---")
        t0 = time.time()
        Qn_hat, Qn_hat_sq, F_extracted = phase1_per_character(
            K, pi, coprime, idx, k, n_steps=20
        )
        print(f"Phase 1 done ({time.time()-t0:.1f}s)")

        # Top-20 characters by |Q_hat_1|
        top_chi_ranks = np.argsort(-Qn_hat_sq[1, :])[:20]
        print(f"\n{'rank':>5} {'chi_j':>6} {'|F(chi)|_extr':>14} {'|lambda|_meas':>14} {'ratio':>8}")
        per_char_rows = []
        for rank, j in enumerate(top_chi_ranks):
            F_val = F_extracted[j]
            lam_mag = abs(eigvals_K_sorted[rank + 1]) if rank + 1 < len(eigvals_K_sorted) else float('nan')
            ratio = F_val / lam_mag if lam_mag > 0 else float('nan')
            print(f"{rank+1:>5} {j:>6} {F_val:>14.6e} {lam_mag:>14.6e} {ratio:>8.3f}")
            per_char_rows.append({
                "rank": rank + 1,
                "chi_j": int(j),
                "F_extracted": F_val,
                "lambda_meas_abs": lam_mag,
                "ratio": ratio,
            })

        # Save per-character CSV
        with open(OUTDIR / f"per_character_rates_k{k}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "chi_j", "F_extracted", "lambda_meas_abs", "ratio"])
            for r in per_char_rows:
                w.writerow([r["rank"], r["chi_j"],
                            f"{r['F_extracted']:.6e}",
                            f"{r['lambda_meas_abs']:.6e}",
                            f"{r['ratio']:.6f}"])

        # =========================================================
        # Phase 2: Jordan blocks
        # =========================================================
        print("\n--- Phase 2: Jordan-block analysis ---")
        t0 = time.time()
        jordan_results = phase2_jordan(K, eigvals_K, k, top_N=10)
        print(f"Phase 2 done ({time.time()-t0:.1f}s)")

        print(f"\n{'|λ|':>10} {'geom':>5} {'alg':>5} {'max_block':>10} {'nullities (j=1..k)':>40}")
        with open(OUTDIR / f"jordan_blocks_k{k}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["lambda_re", "lambda_im", "abs_lambda",
                        "geom_mult", "alg_mult", "max_block_size", "nullities"])
            for jr in jordan_results:
                null_str = " ".join(str(n) for n in jr["nullities"])
                print(f"{jr['abs_lambda']:>10.4e} {jr['geom_mult']:>5} {jr['alg_mult']:>5} "
                      f"{jr['max_block_size']:>10} {null_str:>40}")
                w.writerow([f"{jr['lambda_re']:.6e}", f"{jr['lambda_im']:.6e}",
                            f"{jr['abs_lambda']:.6e}", jr['geom_mult'],
                            jr['alg_mult'], jr['max_block_size'], null_str])

        # =========================================================
        # Phase 3: DS-C mixing time
        # =========================================================
        print("\n--- Phase 3: DS-C mixing time ---")
        t0 = time.time()

        # Cayley diameter under minimal symmetric generating set {2, 2^{-1}}
        gen_set_2 = [2]
        diam_2, balls_2, n_visited = cayley_diameter_and_balls(k, gen_set_2)
        print(f"Cayley diameter under <2>: {diam_2}, |G| visited = {n_visited}, |coprime| = {n}")

        # Moderate growth: |B(r)| / |G| vs (r/diam)^d
        # Fit log(|B(r)|/|G|) = log(A) + d * log(r/γ) for r in [1..γ]
        log_r = np.log(np.arange(1, len(balls_2)) / diam_2)
        log_b = np.log(np.array(balls_2[1:]) / n_visited)
        slope_d, intercept_logA = np.polyfit(log_r, log_b, 1)
        A_const = np.exp(intercept_logA)
        d_growth = slope_d
        print(f"Moderate growth fit: |B(r)|/|G| ≈ {A_const:.3f} · (r/γ)^{d_growth:.3f}")

        # DS-C bound: t_mix(ε) ≤ C(A,d) γ^2 log(1/ε)
        # Standard: C(A,d) ≈ (A^2/d^d) up to constants. We'll report bound at C=1 for a baseline.
        eps = 0.25
        ds_c_bound_C1 = (diam_2 ** 2) * np.log(1.0 / eps)
        print(f"DS-C bound (C=1, ε=1/4): t_mix ≤ γ² log(1/ε) = {diam_2}² · {np.log(1/eps):.3f} = {ds_c_bound_C1:.1f}")

        # Empirical t_mix of K_k
        t_mix_K, starts_K = empirical_mixing_time(K, pi, idx, coprime, max_n=200, eps=eps)
        t_mix_K_arr = np.array([t for t in t_mix_K if t > 0])
        print(f"K_k empirical t_mix(1/4) over 20 starts: median={np.median(t_mix_K_arr):.1f}, "
              f"mean={t_mix_K_arr.mean():.1f}, max={t_mix_K_arr.max()}")

        # Empirical t_mix of B_k
        t_mix_B, starts_B = empirical_mixing_time(B, pi_B, idx, coprime, max_n=2000, eps=eps)
        t_mix_B_arr = np.array([t for t in t_mix_B if t > 0])
        if len(t_mix_B_arr) > 0:
            print(f"B_k empirical t_mix(1/4) over 20 starts: median={np.median(t_mix_B_arr):.1f}, "
                  f"mean={t_mix_B_arr.mean():.1f}, max={t_mix_B_arr.max()}")
        else:
            print(f"B_k empirical t_mix did not converge within 2000 steps")
        print(f"Phase 3 done ({time.time()-t0:.1f}s)")

        # Save mixing CSV
        with open(OUTDIR / f"mixing_time_k{k}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["k", "n_states", "cayley_diameter_gen2", "moderate_growth_A",
                        "moderate_growth_d", "DS_C_bound_C1_eps025",
                        "K_k_t_mix_median", "K_k_t_mix_mean", "K_k_t_mix_max",
                        "B_k_t_mix_median", "B_k_t_mix_mean", "B_k_t_mix_max",
                        "ratio_bound_to_K_k_median", "ratio_bound_to_B_k_median"])
            B_med = float(np.median(t_mix_B_arr)) if len(t_mix_B_arr) > 0 else float('nan')
            B_mean = float(t_mix_B_arr.mean()) if len(t_mix_B_arr) > 0 else float('nan')
            B_max = int(t_mix_B_arr.max()) if len(t_mix_B_arr) > 0 else -1
            K_med = float(np.median(t_mix_K_arr))
            ratio_K = ds_c_bound_C1 / K_med if K_med > 0 else float('nan')
            ratio_B = ds_c_bound_C1 / B_med if B_med > 0 else float('nan')
            w.writerow([k, n, diam_2, f"{A_const:.4f}", f"{d_growth:.4f}",
                        f"{ds_c_bound_C1:.4f}",
                        f"{K_med:.4f}", f"{t_mix_K_arr.mean():.4f}", int(t_mix_K_arr.max()),
                        f"{B_med:.4f}", f"{B_mean:.4f}", B_max,
                        f"{ratio_K:.4f}", f"{ratio_B:.4f}"])

        summary.append({
            "k": k, "n": n,
            "max_F_extracted": float(np.nanmax(F_extracted)),
            "max_meas_lambda": float(abs(eigvals_K_sorted[1])),  # |λ_2|
            "phase1_top_ratio": per_char_rows[0]["ratio"] if per_char_rows else None,
            "phase2_top_block_size": jordan_results[0]["max_block_size"] if jordan_results else None,
            "cayley_diameter": diam_2,
            "DS_C_bound": ds_c_bound_C1,
            "K_k_t_mix_median": float(np.median(t_mix_K_arr)),
            "B_k_t_mix_median": B_med,
        })

    print("\n\n========================================")
    print("SUMMARY")
    print("========================================")
    print(f"{'k':>3} {'n':>5} {'|λ_2|':>10} {'F_top':>10} {'P1 ratio':>10} {'P2 block':>10} "
          f"{'γ_2':>5} {'DS-C':>8} {'K_t_mix':>8} {'B_t_mix':>8}")
    for s in summary:
        print(f"{s['k']:>3} {s['n']:>5} {s['max_meas_lambda']:>10.3e} "
              f"{s['max_F_extracted']:>10.3e} {s['phase1_top_ratio']:>10.3f} "
              f"{s['phase2_top_block_size']:>10} {s['cayley_diameter']:>5} "
              f"{s['DS_C_bound']:>8.1f} {s['K_k_t_mix_median']:>8.1f} {s['B_k_t_mix_median']:>8.1f}")

    print(f"\nOutputs in {OUTDIR}")


if __name__ == "__main__":
    main()
