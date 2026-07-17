"""
qsd_extended_horizon.py — Extend empirical D(r, t) to T=200 with denser
time snapshots, to test whether the trajectory measure asymptotes (QSD)
or genuinely drifts (no QSD).

Walks 10M orbits at N=2^32, records m mod 32 at every 5 steps to T=200.
Critical diagnostic: per-step survival rate vs time.
"""
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")
MAX_VAL = np.int64(2**62)

OUT = Path(r"C:\Collatz")


@njit(parallel=True, cache=True)
def walk_record_residues(starts, T_MAX=200, max_steps=400000):
    n = len(starts)
    residues = np.full((n, T_MAX), -1, dtype=np.int8)
    for i in prange(n):
        m = starts[i]
        residues[i, 0] = m & 31
        steps = 0
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            while (three_m & 1) == 0:
                three_m >>= 1
            m = three_m
            steps += 1
            if steps < T_MAX:
                residues[i, steps] = m & 31
    return residues


def main():
    log2N = 32
    N = 1 << log2N
    n_orbits = 10_000_000
    seed = 137
    T_MAX = 200

    print(f"# Extended horizon: N=2^{log2N}, {n_orbits:,} orbits, T_MAX={T_MAX}")
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    t0 = time.perf_counter()
    residues = walk_record_residues(starts, T_MAX=T_MAX)
    print(f"# Walk: {time.perf_counter() - t0:.1f}s")

    # Compute survival rate at every t
    odd_r32 = list(range(1, 32, 2))
    n_alive_at_t = []
    for t in range(T_MAX):
        n_alive = (residues[:, t] > 0).sum()
        n_alive_at_t.append(n_alive)

    # Survival rate per-step, smoothed over 10-step windows
    print("\n# Per-step survival rate (smoothed over 10-step windows)")
    print("# t-range  n_alive_start  n_alive_end  per_step_survival")
    for t_start in range(0, T_MAX - 10, 10):
        t_end = t_start + 10
        ns = n_alive_at_t[t_start]
        ne = n_alive_at_t[t_end]
        if ns == 0:
            break
        survival = (ne / ns) ** (1 / 10)
        print(f"#  {t_start:>3}-{t_end:>3}  {ns:>13,}  {ne:>11,}  {survival:.6f}")

    # D(r, t) at t = 30, 50, 70, 90, 110, 130, 150, 170, 190
    # Use Chang's pi as denominator (computed elsewhere); for self-contained,
    # just compute rho(r, t) and let the analysis script handle it
    print("\n# rho(r, t) for r in odd mod 32, at sparse t")
    test_t = list(range(0, T_MAX, 10))

    rho_data = []
    for t in test_t:
        r_t = residues[:, t]
        valid = r_t > 0
        r_alive = r_t[valid]
        n_alive = len(r_alive)
        if n_alive < 100:
            continue
        odd_idx = (r_alive & 1).astype(bool)
        r_alive_odd = r_alive[odd_idx]
        counts = np.bincount(r_alive_odd, minlength=32)
        total = counts[1::2].sum()
        rho = {r: counts[r] / total for r in odd_r32}
        rho_data.append((t, n_alive, rho))

    # Print rho table
    print(f"\n# {'t':>3}  {'n_alive':>10}  " +
          "  ".join(f"r={r:>2}" for r in odd_r32))
    for t, n_alive, rho in rho_data:
        line = f"# {t:>3}  {n_alive:>10,}  "
        line += "  ".join(f"{rho[r]:.4f}" for r in odd_r32)
        print(line)

    # Save extended data
    with open(OUT / "qsd_extended.csv", "w") as f:
        f.write("t,n_alive," + ",".join(f"rho_r{r}" for r in odd_r32) + "\n")
        for t, n_alive, rho in rho_data:
            f.write(f"{t},{n_alive}," +
                    ",".join(f"{rho[r]:.8f}" for r in odd_r32) + "\n")
    print(f"\n[wrote] qsd_extended.csv ({len(rho_data)} time snapshots)")

    # Save survival trajectory
    with open(OUT / "qsd_survival_trajectory.csv", "w") as f:
        f.write("t,n_alive,survival_from_0,survival_per_step_smoothed\n")
        for t in range(0, T_MAX, 5):
            ns = n_alive_at_t[t]
            surv = ns / n_orbits
            # Per-step survival smoothed over [t-5, t+5]
            t_lo = max(0, t - 5)
            t_hi = min(T_MAX - 1, t + 5)
            ns_lo = n_alive_at_t[t_lo]
            ns_hi = n_alive_at_t[t_hi]
            if ns_lo > 0 and t_hi - t_lo > 0:
                step_surv = (ns_hi / ns_lo) ** (1 / (t_hi - t_lo))
            else:
                step_surv = float("nan")
            f.write(f"{t},{ns},{surv:.10f},{step_surv:.10f}\n")
    print("[wrote] qsd_survival_trajectory.csv")


if __name__ == "__main__":
    main()
