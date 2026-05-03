"""
Re-walk orbits at N=2^36, track v_t for t=1..T at each Syracuse step.
Stratify by σ-quartile (bottom 25% of σ_resid). For bottom-quartile orbits,
compute E[v_t | q125] and P(v_t = k | q125) at each t.

If E[v_t | q125] is flat ~2.22 for all t: orbits keep their fast-descent
character → multi-step contribution structurally derivable as t-fold extension
of single-step density amplification.

If E[v_t | q125] decays toward unconditional Geom(1/2) E[v]=2: orbits "forget"
their initial fast-descent advantage at some mixing rate λ. Multi-step piece
is a relaxation calculation.

Either result gives closure structure for the residual ~20% of E[v]_q125.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_syracuse_with_v_seq(starts, max_value, max_syr_steps, T_track):
    """Walk Syracuse iteration. Record:
      - sigma_arr[i]: total Collatz step count (odd steps + halvings)
      - n_odd_arr[i]: total Syracuse step count
      - v_seq[i, t]: v_2(3m+1) at Syracuse step t (or -1 if t ≥ n_odd)
      - ok_arr[i]: True iff orbit completed within max_syr_steps
    """
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    n_odd_arr = np.zeros(n, dtype=np.int32)
    v_seq = np.full((n, T_track), -1, dtype=np.int8)
    ok_arr = np.zeros(n, dtype=np.bool_)

    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0
        syr_steps = 0
        failed = False
        while m != 1 and syr_steps < max_syr_steps:
            # m must be odd here (Syracuse iterates odd → odd)
            # If m even (only possible if start is even, which we don't use), halve
            if (m & 1) == 0:
                m = m >> 1
                sigma_total += 1
                continue
            if m > max_value // 3:
                failed = True; break
            threex_p1 = 3 * m + 1
            v = 0
            tmp = threex_p1
            while (tmp & 1) == 0:
                tmp >>= 1
                v += 1
            if syr_steps < T_track:
                v_seq[i, syr_steps] = v
            m = tmp
            sigma_total += 1 + v
            syr_steps += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            n_odd_arr[i] = syr_steps
            ok_arr[i] = True
    return sigma_arr, n_odd_arr, v_seq, ok_arr


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N = 36
    N = 1 << log2N
    n_starts = 500_000
    T_track = 50    # track v_t for t = 0..49

    print(f"# v_t conditional distribution along bottom-σ-quartile orbits", flush=True)
    print(f"# N = 2^{log2N} = {N:,}, sampling {n_starts:,} orbits, tracking T={T_track} steps", flush=True)
    print(f"# Reference: unconditional Geom(1/2) has E[v] = 2.0", flush=True)
    print(f"# Asymptote target (Result 14 follow-up 2): E[v]_q125,∞ = 2.216", flush=True)

    seeds = [42, 137, 271, 314, 1729]
    n_per_seed = 100_000  # 5 seeds × 100K = 500K total

    # Aggregate across seeds
    P_v_uncond = np.zeros((T_track, 30))  # P(v_t = k | uncond) for k=0..29
    P_v_cond = np.zeros((T_track, 30))
    counts_uncond = np.zeros(T_track, dtype=np.int64)
    counts_cond = np.zeros(T_track, dtype=np.int64)
    sum_v_uncond = np.zeros(T_track)
    sum_v_cond = np.zeros(T_track)
    sumsq_v_uncond = np.zeros(T_track)
    sumsq_v_cond = np.zeros(T_track)

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_per_seed, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, n_odd, v_seq, ok = walk_syracuse_with_v_seq(
            starts, MAX_VAL, 1_000_000, T_track)
        elapsed = time.perf_counter() - t0
        print(f"  seed={seed}: walked in {elapsed:.1f}s, ok={int(ok.sum()):,}, failed={int((~ok).sum())}", flush=True)

        starts_ok = starts[ok]
        sigma_ok = sigma[ok].astype(np.float64)
        n_odd_ok = n_odd[ok]
        v_seq_ok = v_seq[ok]

        # σ_resid
        log_n = np.log(starts_ok.astype(np.float64))
        log_n_c = log_n - log_n.mean()
        sigma_c = sigma_ok - sigma_ok.mean()
        beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
        alpha = float(sigma_ok.mean() - beta * log_n.mean())
        sigma_resid = sigma_ok - (alpha + beta * log_n)
        q25 = float(np.percentile(sigma_resid, 25))
        mask_low = sigma_resid <= q25

        # For each t, aggregate v_t over orbits with n_odd > t (orbit hasn't terminated)
        for t in range(T_track):
            valid = n_odd_ok > t  # orbits still running at step t
            v_at_t = v_seq_ok[valid, t].astype(np.int32)

            # Unconditional
            counts_uncond[t] += len(v_at_t)
            sum_v_uncond[t] += float(v_at_t.sum())
            sumsq_v_uncond[t] += float((v_at_t.astype(np.int64) ** 2).sum())
            for k in range(30):
                P_v_uncond[t, k] += float((v_at_t == k).sum())

            # Conditional on bottom-σ-quartile
            valid_low = valid & mask_low
            v_at_t_low = v_seq_ok[valid_low, t].astype(np.int32)
            counts_cond[t] += len(v_at_t_low)
            sum_v_cond[t] += float(v_at_t_low.sum())
            sumsq_v_cond[t] += float((v_at_t_low.astype(np.int64) ** 2).sum())
            for k in range(30):
                P_v_cond[t, k] += float((v_at_t_low == k).sum())

    # Compute final means and probabilities
    E_v_uncond = sum_v_uncond / np.maximum(counts_uncond, 1)
    E_v_cond = sum_v_cond / np.maximum(counts_cond, 1)
    Var_v_uncond = sumsq_v_uncond / np.maximum(counts_uncond, 1) - E_v_uncond**2
    Var_v_cond = sumsq_v_cond / np.maximum(counts_cond, 1) - E_v_cond**2
    P_v_uncond = P_v_uncond / np.maximum(counts_uncond[:, None], 1)
    P_v_cond = P_v_cond / np.maximum(counts_cond[:, None], 1)

    # Print E[v_t] trajectory
    print(f"\n# E[v_t] vs t trajectory (T = {T_track}):", flush=True)
    print(f"  {'t':>3}  {'E[v_t] uncond':>14}  {'E[v_t] cond':>13}  {'shift':>9}  "
          f"{'Var cond':>10}  {'#cond':>10}  {'#uncond':>10}", flush=True)
    for t in range(T_track):
        if counts_cond[t] < 100:
            break
        print(f"  {t:>3}  {E_v_uncond[t]:>14.5f}  {E_v_cond[t]:>13.5f}  "
              f"{E_v_cond[t]-E_v_uncond[t]:>+9.5f}  {Var_v_cond[t]:>10.4f}  "
              f"{counts_cond[t]:>10,}  {counts_uncond[t]:>10,}", flush=True)

    # Decay-rate fit if E[v_t | cond] is decreasing toward unconditional
    valid_t = counts_cond >= 1000
    t_arr = np.arange(T_track)[valid_t]
    E_v_cond_arr = E_v_cond[valid_t]
    E_v_uncond_arr = E_v_uncond[valid_t]
    # Detrend by E[v]_uncond — model: (E_cond - E_uncond) = A · exp(-λ·t)
    diff = E_v_cond_arr - E_v_uncond_arr
    # Only fit where diff > 0
    pos = diff > 1e-3
    if pos.sum() >= 5:
        log_diff = np.log(diff[pos])
        t_fit = t_arr[pos]
        # Linear regression log_diff = log_A - λ · t
        xc = t_fit - t_fit.mean()
        yc = log_diff - log_diff.mean()
        slope = float((xc * yc).sum() / (xc * xc).sum())
        intercept = float(log_diff.mean() - slope * t_fit.mean())
        A = float(np.exp(intercept))
        lam = -slope
        print(f"\n# Exponential decay fit: E[v_t | cond] - E[v_t | uncond] ≈ A · exp(-λt)", flush=True)
        print(f"  A = {A:.4f}", flush=True)
        print(f"  λ = {lam:.5f}", flush=True)
        print(f"  Half-life (steps): {np.log(2) / max(lam, 1e-9):.2f}", flush=True)
        # If decays to 0: orbits forget initial fast-descent
        # If A is small but persistent: stays elevated structurally

    # Average over orbit lifetime: E[V_per_orbit | q125] should = mean(E[v_t | q125]) over t weighted by orbit-length distribution
    # Estimate empirical E[V] over conditional orbits
    print(f"\n# Cross-check: avg E[v_t | q125] over t = 0..{T_track-1}, weighted by counts:", flush=True)
    avg_uncond = float((sum_v_uncond.sum()) / max(counts_uncond.sum(), 1))
    avg_cond = float((sum_v_cond.sum()) / max(counts_cond.sum(), 1))
    print(f"  Total-step weighted E[v] uncond = {avg_uncond:.5f}", flush=True)
    print(f"  Total-step weighted E[v] cond   = {avg_cond:.5f}", flush=True)
    print(f"  Shift: {avg_cond - avg_uncond:+.5f}", flush=True)
    print(f"  Reference (Result 14 ff.2): empirical E[V]_q125 at 2^36 = ~2.32, shift = +0.27", flush=True)
    print(f"                              asymptote target = 2.216, shift = +0.221", flush=True)

    # Show P(v_t = k) for selected t values
    print(f"\n# P(v_t = k | bottom q) at selected t (k=1..6):", flush=True)
    print(f"  {'t':>3}", end="", flush=True)
    for k in range(1, 7):
        print(f"  {'P(v='+str(k)+')':>9}", end="", flush=True)
    print("", flush=True)
    for t in [0, 1, 2, 3, 5, 10, 20, 30, 49]:
        if t >= T_track or counts_cond[t] < 100: continue
        print(f"  {t:>3}", end="", flush=True)
        for k in range(1, 7):
            print(f"  {P_v_cond[t, k]:>9.4f}", end="", flush=True)
        print("", flush=True)

    # Save
    rows = []
    for t in range(T_track):
        if counts_cond[t] < 100: continue
        rows.append({
            't': t,
            'E_v_uncond': float(E_v_uncond[t]),
            'E_v_cond': float(E_v_cond[t]),
            'shift': float(E_v_cond[t] - E_v_uncond[t]),
            'Var_v_cond': float(Var_v_cond[t]),
            'counts_uncond': int(counts_uncond[t]),
            'counts_cond': int(counts_cond[t]),
        })
    out_csv = out_dir / "58_v_t_along_q125.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}", flush=True)

    # Plot if matplotlib available
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        T_show = min(T_track, len([t for t in range(T_track) if counts_cond[t] >= 100]))
        ts = np.arange(T_show)
        axes[0].plot(ts, E_v_uncond[:T_show], 'b-', label='E[v_t] uncond', alpha=0.7)
        axes[0].plot(ts, E_v_cond[:T_show], 'r-', label='E[v_t] | bottom-σ-quartile', alpha=0.9)
        axes[0].axhline(2.0, color='gray', linestyle=':', alpha=0.5, label='Geom(1/2) E[v]=2.0')
        axes[0].axhline(2.216, color='purple', linestyle=':', alpha=0.5, label='asymptote 2.216')
        axes[0].set_xlabel('t (Syracuse step)')
        axes[0].set_ylabel('E[v_t]')
        axes[0].legend()
        axes[0].grid(alpha=0.3)
        axes[0].set_title(f'E[v_t] vs t at N=2^{log2N}')

        axes[1].semilogy(ts, np.maximum(E_v_cond[:T_show] - E_v_uncond[:T_show], 1e-4), 'r-')
        axes[1].set_xlabel('t')
        axes[1].set_ylabel('E[v_t | cond] - E[v_t | uncond] (log scale)')
        axes[1].grid(alpha=0.3, which='both')
        axes[1].set_title('Conditional shift vs t (log scale; linear = exponential decay)')

        plt.tight_layout()
        plot_path = out_dir / "58_v_t_along_q125.png"
        plt.savefig(plot_path, dpi=120)
        print(f"[save] {plot_path}", flush=True)
    except Exception as e:
        print(f"  (plot skipped: {e})", flush=True)


if __name__ == "__main__":
    main()
