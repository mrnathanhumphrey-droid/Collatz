"""
generate_descent_data.py - 3x+1 descent data for two stacked visualizations.

Outputs to C:/Collatz/visualization/viz_outputs/:
  descent_a.csv  per-step trajectory points (~5M rows)
                 columns: x=log(value at step t), y=step_idx t, z=a_star_idx,
                          a_star_idx (== z for color), status
  descent_b.csv  per-orbit summary (50K rows)
                 columns: x=log(n_start), y=sigma, z=log(peak/n_start),
                          a_star_idx (color), status

q=3 only. All orbits expected to converge (3x+1 verified through n ~ 2^68).
"""
import csv
import importlib.util
import time
from pathlib import Path

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
EXP29 = ROOT / "experiments" / "29_qx1_cycle_classification.py"
spec = importlib.util.spec_from_file_location("exp29", EXP29)
exp29 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp29)
qx1_prefix = exp29.qx1_prefix

# Config
N_SAMPLES = 50_000
N_MAX = 1_000_000
K_MOD = 6
MAX_STEPS = 5_000
SUBSAMPLE_VIS1 = 100
SEED = 42

OUT_DIR = Path(__file__).parent / "viz_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DESCENT_A = OUT_DIR / "descent_a.csv"
DESCENT_B = OUT_DIR / "descent_b.csv"

S_CONVERGED = 0
S_TIMEOUT = 1


@njit(cache=True)
def trajectory_3x1(n_start, max_steps, log_buf):
    """Walk T_3 from n_start. Fill log_buf with log(value) at each step.
    Returns (n_logged, status, sigma, peak_value).
    """
    n = n_start
    log_buf[0] = np.log(np.float64(n))
    n_logged = 1
    peak = n
    for step in range(1, max_steps + 1):
        if n == 1:
            return n_logged, S_CONVERGED, step - 1, peak
        if n % 2 == 0:
            n = n >> 1
        else:
            n = 3 * n + 1
        if n > peak:
            peak = n
        log_buf[step] = np.log(np.float64(n))
        n_logged += 1
    return n_logged, S_TIMEOUT, max_steps, peak


def a_star_idx_for_n(n, k=K_MOD):
    """Return j where a_star = 3^j (a_star_idx in {0..k})."""
    r = n % (1 << k)
    if r == 0:
        return 0
    _, a_star, _ = qx1_prefix(r, k, 3)
    j = 0
    a = a_star
    while a > 1 and a % 3 == 0:
        a //= 3
        j += 1
    if a != 1:
        j = -1
    return j


def main():
    rng = np.random.default_rng(SEED)
    odd_pool = np.arange(3, N_MAX + 1, 2, dtype=np.int64)
    starts = rng.choice(odd_pool, size=N_SAMPLES, replace=False)
    starts.sort()
    log_buf = np.zeros(MAX_STEPS + 1, dtype=np.float64)

    # JIT warmup
    _ = trajectory_3x1(np.int64(7), MAX_STEPS, log_buf)

    a_count = 0
    b_count = 0
    timeout_count = 0
    sigma_max = 0
    peak_log_max = 0.0
    peak_ratio_max = 0.0

    t0 = time.perf_counter()

    with open(DESCENT_A, "w", newline="") as fa, open(DESCENT_B, "w", newline="") as fb:
        wa = csv.writer(fa)
        wb = csv.writer(fb)
        wa.writerow(["x", "y", "z", "a_star_idx", "status"])
        wb.writerow(["x", "y", "z", "a_star_idx", "status"])

        for n in starts:
            n_int = int(n)
            n_logged, status, sigma, peak = trajectory_3x1(np.int64(n), MAX_STEPS, log_buf)
            if status == S_TIMEOUT:
                timeout_count += 1
                continue  # skip non-convergent (expect 0 for n <= 10^6 under 3x+1)
            j = a_star_idx_for_n(n_int)
            log_n = float(np.log(n_int))
            log_peak = float(np.log(peak))
            log_ratio = log_peak - log_n  # log(peak / n_start)

            sigma_max = max(sigma_max, sigma)
            peak_log_max = max(peak_log_max, log_peak)
            peak_ratio_max = max(peak_ratio_max, log_ratio)

            # Vis 1: per-step trajectory points (subsample to SUBSAMPLE_VIS1 if longer)
            if n_logged > SUBSAMPLE_VIS1:
                idxs = np.linspace(0, n_logged - 1, SUBSAMPLE_VIS1).astype(np.int64)
            else:
                idxs = np.arange(n_logged, dtype=np.int64)
            for t_idx in idxs:
                wa.writerow((
                    f"{log_buf[t_idx]:.6f}",
                    int(t_idx),
                    j,
                    j,
                    "converged",
                ))
                a_count += 1

            # Vis 2: per-orbit (log_n_start, sigma, log(peak/n_start))
            wb.writerow((
                f"{log_n:.6f}",
                int(sigma),
                f"{log_ratio:.6f}",
                j,
                "converged",
            ))
            b_count += 1

    elapsed = time.perf_counter() - t0
    print(f"\n[gen] DONE in {elapsed:.1f}s")
    print(f"[gen] descent_a.csv rows: {a_count:,}  ({DESCENT_A})")
    print(f"[gen] descent_b.csv rows: {b_count:,}  ({DESCENT_B})")
    print(f"[gen] timeouts (skipped): {timeout_count}")
    print(f"[gen] max sigma: {sigma_max}  max log(peak): {peak_log_max:.2f}  max log(peak/start): {peak_ratio_max:.2f}")


if __name__ == "__main__":
    main()
