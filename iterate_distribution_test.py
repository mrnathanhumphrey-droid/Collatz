"""
iterate_distribution_test.py — probe the trajectory measure via fixed-iterate
distributions.

Three tests:
  (1) P(v_t = k | step t, σ_S > t) — per-step v-distribution at iterate t.
      Does it converge to Geom(1/2), stay biased, or evolve?
  (2) P(m_t mod 32 | t, σ_S > t) — residue distribution at iterate t.
      Connects to Result 23's forward-chain λ_max = 1.264.
  (3) Mellin M[P(m_t)](s) — multiplicative Mellin at fixed t.
      Look for ζ-like structure that σ-marginal washed out.

Setup: Walk 1M orbits at N=2^32. At each step t, record (m_t, v_t).
Bucket by t, condition on "still alive" (σ_S > t).

Compute time: ~30s.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_record_steps(starts, T_max=80, max_total_steps=400000):
    """Walk Syracuse, record m_t at fixed steps t = 0, 1, ..., T_max - 1.
    For each t, also record v_t (the v at step t→t+1).
    If orbit absorbed before step t: m_t = -1 (sentinel).
    """
    n = len(starts)
    # m_at_step[i, t] = m at step t for orbit i (or -1 if absorbed)
    m_at_step = np.full((n, T_max), -1, dtype=np.int64)
    # v_at_step[i, t] = v applied at step t→t+1 (or -1 if no step)
    v_at_step = np.full((n, T_max), -1, dtype=np.int8)
    sigma_arr = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        m = starts[i]
        m_at_step[i, 0] = m
        steps = 0
        while m != 1 and steps < max_total_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            v = 0
            tmp = three_m
            while (tmp & 1) == 0:
                tmp >>= 1
                v += 1
            if steps < T_max:
                v_at_step[i, steps] = min(v, 30)  # cap v at 30
            m = tmp
            steps += 1
            if steps < T_max:
                m_at_step[i, steps] = m
        sigma_arr[i] = steps
    return m_at_step, v_at_step, sigma_arr


def main():
    log2N = 32
    N = 1 << log2N
    n_orbits = 1_000_000
    seed = 42
    T_max = 80

    print(f"# Iterate-distribution test at N=2^{log2N}, {n_orbits:,} orbits")
    print(f"# Track m_t and v_t up to step T_max = {T_max}\n")

    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    t0 = time.perf_counter()
    m_at_step, v_at_step, sigma = walk_record_steps(starts, T_max=T_max)
    elapsed = time.perf_counter() - t0
    print(f"# Walk completed in {elapsed:.1f}s\n")

    # === Test 1: P(v_t = k | step t, alive) ===
    print(f"# Test 1: P(v_t = k | step t, orbit still alive)")
    print(f"# Geom(1/2) baseline: P(v=k) = 2^(-k) for k = 1, 2, ...")
    print(f"# Stage 1 found v=4 spike (P(v=4) ≈ 1.37 × 2^(-4) under uniform-m sampling)")
    print()
    print(f"#   {'t':>3}  {'n_alive':>9}  {'P(v=1)':>8}  {'P(v=2)':>8}  {'P(v=3)':>8}  "
          f"{'P(v=4)':>8}  {'P(v=5)':>8}  {'P(v=6)':>8}  {'⟨v⟩':>6}")
    test_t_values = [0, 1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70]
    v_dist_per_t = {}
    for t in test_t_values:
        if t >= T_max - 1:
            continue
        v_t = v_at_step[:, t]
        valid = v_t > 0  # has a step
        v_t = v_t[valid]
        n_alive = len(v_t)
        if n_alive < 100:
            continue
        # P(v=k)
        max_v = int(v_t.max())
        counts = np.bincount(v_t, minlength=max_v + 1)
        P = counts / counts.sum()
        v_dist_per_t[t] = P
        Ev = (np.arange(len(P)) * P).sum()
        line = f"#   {t:>3}  {n_alive:>9,}  "
        for k in [1, 2, 3, 4, 5, 6]:
            if k < len(P):
                line += f"{P[k]:>8.5f}  "
            else:
                line += f"{0:>8.5f}  "
        line += f"{Ev:>6.4f}"
        print(line)

    # Compare to Geom(1/2) and check v=4 spike per t
    print(f"\n# v=4 spike strength per t (P(v=4) / Geom(1/2) = P(v=4) · 16):")
    print(f"#   {'t':>3}  {'P(v=4)/Geom':>15}")
    for t in test_t_values:
        if t in v_dist_per_t:
            P = v_dist_per_t[t]
            if 4 < len(P):
                spike = P[4] * 16
                print(f"#   {t:>3}  {spike:>15.4f}")

    # === Test 2: P(m_t mod 32 | t, alive) — residue distribution evolution ===
    print(f"\n# Test 2: P(m_t mod 32 = r | step t, alive)")
    print(f"# Tracks how residue distribution evolves under iterated Syracuse.")
    print(f"# Result 23 (forward chain) found stationary distribution; check empirical convergence.\n")

    print(f"#   {'t':>3}  {'n_alive':>9}  CV(P(r))  range(P(r))  argmax(P)  argmin(P)")
    for t in test_t_values:
        if t >= T_max:
            continue
        m_t = m_at_step[:, t]
        valid = m_t > 0
        m_t = m_t[valid]
        n_alive = len(m_t)
        if n_alive < 100:
            continue
        residues = (m_t & 31).astype(np.int64)  # mod 32
        # Only odd residues (m is always odd in Syracuse)
        odd_res = residues[(residues & 1).astype(bool)]
        odd_res_idx = (odd_res - 1) // 2  # map {1,3,...,31} to {0..15}
        counts = np.bincount(odd_res_idx, minlength=16)
        P = counts / counts.sum()
        cv = P.std() / P.mean()
        argmax = (np.argmax(P) * 2 + 1)
        argmin = (np.argmin(P) * 2 + 1)
        rng_p = P.max() - P.min()
        print(f"#   {t:>3}  {n_alive:>9,}  {cv:>8.4f}  [{P.min():.4f}, {P.max():.4f}]  "
              f"r={argmax:>2}    r={argmin:>2}")

    # === Test 3: Mellin of P(m_t) at large t ===
    print(f"\n# Test 3: Mellin of P(m_t) at fixed step t (probe multiplicative structure)")
    print(f"#   For each t, compute M[P(m_t)](s) = Σ_m P(m_t) · m^(s-1)")
    print(f"#   Look for poles, zeros on critical line s = 1/2 + it_imag")
    print()
    for t in [5, 10, 20, 40]:
        if t >= T_max:
            continue
        m_t = m_at_step[:, t]
        valid = m_t > 0
        m_t = m_t[valid].astype(np.float64)
        n_alive = len(m_t)
        if n_alive < 1000:
            continue
        log_m = np.log(m_t)
        print(f"#   t={t}: n_alive = {n_alive:,}, ⟨log m_t⟩ = {log_m.mean():.4f} "
              f"(vs log m_0 ≈ {(log2N - 1) * LOG2:.4f}), sd(log m_t) = {log_m.std(ddof=1):.4f}")
        # Mellin via empirical sum
        # Need histogram of m_t. Since m_t can be large integers, use empirical sum
        # M(s) = (1/n_alive) Σ_i m_i^(s-1)
        # Real s
        for s_real in [0.5, 1.5, 2.0]:
            M = float(np.mean(m_t ** (s_real - 1)))
            print(f"#     M({s_real}) = (1/n) Σ m^({s_real-1}) = {M:.4e}")
        # Imaginary axis
        t_grid = np.array([1.0, 2.0, 5.0, 10.0, 20.0])
        for t_im in t_grid:
            s = 1j * t_im
            # m^(s-1) = m^(-1) · m^(it_im) = m^(-1) · exp(it_im · log m)
            M_complex = np.mean(m_t ** (-1.0) * np.exp(1j * t_im * log_m))
            print(f"#     |M(it={t_im})| = {abs(M_complex):.4e}, arg = {np.angle(M_complex):.3f}")
        # Critical line
        for t_im in [5.0, 10.0, 20.0]:
            s = 0.5 + 1j * t_im
            M_complex = np.mean(m_t ** (-0.5) * np.exp(1j * t_im * log_m))
            print(f"#     |M(1/2+i{t_im})| = {abs(M_complex):.4e}")

    # === Test 4: Fourier of P(v=k) sequence — look for periodic structure ===
    print(f"\n# Test 4: Fourier of P(v=·) sequence at t=0 and t=20")
    print(f"#   F(P)(ω) = Σ_k P(v=k) · exp(2πi·k·ω)")
    print(f"#   Geom(1/2) gives F(ω) = (1/2)/(1 - exp(2πiω)/2) = closed form\n")
    for t in [0, 20]:
        if t in v_dist_per_t:
            P = v_dist_per_t[t]
            print(f"#   t={t}: P = [P(v=1), ..., P(v=10)] = {[f'{P[k]:.4f}' for k in range(1, min(11, len(P)))]}")
            # FFT
            P_padded = np.zeros(32)
            P_padded[:len(P)] = P
            F = np.fft.fft(P_padded)
            print(f"#     |F(ω)| at ω = 0, 1/8, 1/4, 1/2:")
            for omega_idx in [0, 4, 8, 16]:
                if omega_idx < len(F):
                    print(f"#       ω={omega_idx/32:.3f}: |F| = {abs(F[omega_idx]):.4f}, "
                          f"arg = {np.angle(F[omega_idx]):.3f}")

    # Save summary CSV
    rows = []
    for t in test_t_values:
        if t in v_dist_per_t:
            P = v_dist_per_t[t]
            rows.append({
                't': t,
                'P_v_1': float(P[1]) if len(P) > 1 else 0,
                'P_v_2': float(P[2]) if len(P) > 2 else 0,
                'P_v_3': float(P[3]) if len(P) > 3 else 0,
                'P_v_4': float(P[4]) if len(P) > 4 else 0,
                'P_v_5': float(P[5]) if len(P) > 5 else 0,
                'P_v_4_spike': float(P[4] * 16) if len(P) > 4 else 0,
                'mean_v': float((np.arange(len(P)) * P).sum()),
            })
    pl.DataFrame(rows).write_csv("experiments_output/iterate_v_distribution.csv")
    print(f"\n[save] experiments_output/iterate_v_distribution.csv")


if __name__ == "__main__":
    main()
