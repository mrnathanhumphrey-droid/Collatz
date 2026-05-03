"""
generate_cumv_data.py - Vis 1 with cumulative-v coloring (Cramer race tracker).

For each step of each trajectory, track:
  cum_v = sum of v_2(3n+1) over all odd steps so far
  cum_v_heuristic = 2 * n_odd_steps (heuristic E[v_2] = 2)
  cumv_residual = cum_v - cum_v_heuristic  (positive = orbit "luckier" than heuristic)

Color trajectory points by cumv_residual to visualize where each orbit
beats vs underperforms the Crandall heuristic.

Output: viz_outputs/cumv_a.csv with same column shape as descent_a.csv plus a
'cumv_resid' column for coloring.
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

N_SAMPLES = 50_000
N_MAX = 1_000_000
K_MOD = 6
MAX_STEPS = 5_000
SUBSAMPLE = 100
SEED = 42

OUT_DIR = Path(__file__).parent / "viz_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CUMV_CSV = OUT_DIR / "cumv_a.csv"


@njit(cache=True)
def trajectory_with_cumv(n_start, max_steps, log_buf, cumv_buf, n_odd_buf):
    """Walk T_3, recording log(value), cumulative v, and odd-step count at each step."""
    n = n_start
    log_buf[0] = np.log(np.float64(n))
    cumv_buf[0] = 0
    n_odd_buf[0] = 0
    n_logged = 1
    cum_v = 0
    n_odd = 0
    for step in range(1, max_steps + 1):
        if n == 1:
            return n_logged
        if n % 2 == 0:
            n = n >> 1
        else:
            # Odd step: m = 3n+1, then strip factors of 2 ourselves to count v
            n_odd += 1
            m = 3 * n + 1
            v = 0
            while m % 2 == 0:
                m >>= 1
                v += 1
            cum_v += v
            # The current step is the multiplication (n -> 3n+1). We then take v halves
            # implicitly. Record the value AFTER the multiplication, then advance.
            # To keep step granularity comparable to the original walker, record state
            # at each individual halve too.
            # First, record the (3n+1) state:
            log_buf[n_logged] = np.log(np.float64(3 * n + 1))
            cumv_buf[n_logged] = cum_v
            n_odd_buf[n_logged] = n_odd
            n_logged += 1
            # Then record each halve:
            interim = 3 * n + 1
            for _ in range(v):
                interim >>= 1
                if n_logged >= log_buf.shape[0]:
                    return n_logged
                log_buf[n_logged] = np.log(np.float64(interim))
                cumv_buf[n_logged] = cum_v
                n_odd_buf[n_logged] = n_odd
                n_logged += 1
            n = interim
            continue
        log_buf[n_logged] = np.log(np.float64(n))
        cumv_buf[n_logged] = cum_v
        n_odd_buf[n_logged] = n_odd
        n_logged += 1
    return n_logged


def a_star_idx_for_n(n, k=K_MOD):
    r = n % (1 << k)
    if r == 0:
        return 0
    _, a_star, _ = qx1_prefix(r, k, 3)
    j = 0
    a = a_star
    while a > 1 and a % 3 == 0:
        a //= 3
        j += 1
    return j if a == 1 else -1


def main():
    rng = np.random.default_rng(SEED)
    odd_pool = np.arange(3, N_MAX + 1, 2, dtype=np.int64)
    starts = rng.choice(odd_pool, size=N_SAMPLES, replace=False)
    starts.sort()

    BUF_SIZE = MAX_STEPS + 200  # allow extra for halve-cascades on rare big v
    log_buf = np.zeros(BUF_SIZE, dtype=np.float64)
    cumv_buf = np.zeros(BUF_SIZE, dtype=np.int64)
    n_odd_buf = np.zeros(BUF_SIZE, dtype=np.int64)

    # JIT warmup
    _ = trajectory_with_cumv(np.int64(7), MAX_STEPS, log_buf, cumv_buf, n_odd_buf)

    a_count = 0
    t0 = time.perf_counter()
    skip = 0

    with open(CUMV_CSV, "w", newline="") as fa:
        wa = csv.writer(fa)
        wa.writerow(["x", "y", "z", "cumv_resid", "a_star_idx", "status"])

        for n in starts:
            n_int = int(n)
            n_logged = trajectory_with_cumv(np.int64(n), MAX_STEPS, log_buf, cumv_buf, n_odd_buf)
            if n_logged >= BUF_SIZE - 1:
                skip += 1
                continue
            j = a_star_idx_for_n(n_int)

            # Subsample to SUBSAMPLE points
            if n_logged > SUBSAMPLE:
                idxs = np.linspace(0, n_logged - 1, SUBSAMPLE).astype(np.int64)
            else:
                idxs = np.arange(n_logged, dtype=np.int64)
            for t_idx in idxs:
                cumv = int(cumv_buf[t_idx])
                n_odd = int(n_odd_buf[t_idx])
                cumv_pred = 2 * n_odd  # heuristic E[v]=2 per odd step
                cumv_resid = cumv - cumv_pred
                wa.writerow((
                    f"{log_buf[t_idx]:.6f}",
                    int(t_idx),
                    j,
                    cumv_resid,
                    j,
                    "converged",
                ))
                a_count += 1

    elapsed = time.perf_counter() - t0
    print(f"\n[gen] DONE in {elapsed:.1f}s")
    print(f"[gen] cumv_a.csv rows: {a_count:,}  ({CUMV_CSV})")
    print(f"[gen] skipped (overflowed step buffer): {skip}")


if __name__ == "__main__":
    main()
