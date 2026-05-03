"""
generate_viz_data.py - qx+1 phase-space data for 3D visualization.

Outputs two CSVs to C:/Collatz/visualization/:
  vis_a.csv  trajectory point cloud (one row per visible step)
              columns: x=log(value), y=step_idx, z=q, a_star, a_star_power, status
  vis_b.csv  stopping-time summary (one row per starting n)
              columns: x=log(n), y=sigma (or -1 sentinel), z=q, a_star, a_star_power, status
  vis_b_div.csv  divergence-aware summary (one row per starting n, all 60K plottable)
              columns: x=log(n), y=log(final_value_at_termination), z=q, a_star, a_star_power, status
              For converged orbits y=0 (since final value = 1). For divergent/timeout, y = log of
              last value before the orbit hit the int64 safe-cap or step cap.

Reuses qx1_prefix from experiments/29_qx1_cycle_classification.py for a* lookup.
The a_star_power column = j where a_star = q^j; this is the categorical invariant
for coloring across q (raw a_star = q^j varies wildly across q).
"""
import csv
import importlib.util
import time
from pathlib import Path

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
EXP29 = ROOT / "experiments" / "29_qx1_cycle_classification.py"

# Import qx1_prefix from the digit-prefixed source file
spec = importlib.util.spec_from_file_location("exp29", EXP29)
exp29 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp29)
qx1_prefix = exp29.qx1_prefix  # (r, k, q) -> (prefix_steps, a_star, c_star)

# Config
Q_VALUES = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
N_PER_Q = 5000
N_MAX = 1_000_000
K_MOD = 6
MAX_STEPS = 10_000
MAX_VALUE_INT = 1 << 62  # ~4.6e18, leaves headroom for q*x+1 in int64
VIS_SUBSAMPLE = 200
SEED = 42

OUT_DIR = Path(__file__).parent
VIS_A_CSV = OUT_DIR / "vis_a.csv"
VIS_B_CSV = OUT_DIR / "vis_b.csv"
VIS_B_DIV_CSV = OUT_DIR / "vis_b_div.csv"

S_CONVERGED = 0
S_DIVERGENT = 1
S_TIMEOUT   = 2
STATUS_NAMES = {0: "converged", 1: "divergent", 2: "timeout"}


@njit(cache=True)
def trajectory_qx1(n_start, q, max_steps, max_value_int, log_buf):
    """Walk T_q from n_start, recording log(value) at each step into log_buf.
    Returns (n_logged, status, sigma).
    """
    safe_cap = max_value_int // q
    n = n_start
    log_buf[0] = np.log(np.float64(n))
    n_logged = 1
    for step in range(1, max_steps + 1):
        if n == 1:
            return n_logged, S_CONVERGED, step - 1
        if n % 2 == 0:
            n = n >> 1
        else:
            if n > safe_cap:
                return n_logged, S_DIVERGENT, step - 1
            n = q * n + 1
        if n > max_value_int:
            return n_logged, S_DIVERGENT, step
        log_buf[step] = np.log(np.float64(n))
        n_logged += 1
    return n_logged, S_TIMEOUT, max_steps


def a_star_for_n(n, k, q):
    """Return (a_star, j) where a_star = q^j is the qx+1 prefix's terminal odd a-coef."""
    r = n % (1 << k)
    if r == 0:
        return 1, 0
    _, a_star, _ = qx1_prefix(r, k, q)
    j = 0
    a = a_star
    while a > 1 and a % q == 0:
        a //= q
        j += 1
    if a != 1:
        j = -1  # shouldn't happen for qx+1 prefix; defensive
    return a_star, j


def main():
    rng = np.random.default_rng(SEED)
    odd_pool = np.arange(3, N_MAX + 1, 2, dtype=np.int64)
    log_buf = np.zeros(MAX_STEPS + 1, dtype=np.float64)

    # Warm up JIT
    _ = trajectory_qx1(np.int64(7), np.int64(3), MAX_STEPS, MAX_VALUE_INT, log_buf)

    vis_a_count = 0
    vis_b_count = 0
    status_counts = {q: {"converged": 0, "divergent": 0, "timeout": 0} for q in Q_VALUES}
    sigma_max_seen = 0

    t0 = time.perf_counter()

    with open(VIS_A_CSV, "w", newline="") as fa, \
         open(VIS_B_CSV, "w", newline="") as fb, \
         open(VIS_B_DIV_CSV, "w", newline="") as fbd:
        wa = csv.writer(fa)
        wb = csv.writer(fb)
        wbd = csv.writer(fbd)
        wa.writerow(["x", "y", "z", "a_star", "a_star_power", "status"])
        wb.writerow(["x", "y", "z", "a_star", "a_star_power", "status"])
        wbd.writerow(["x", "y", "z", "a_star", "a_star_power", "status"])

        for q in Q_VALUES:
            starts = rng.choice(odd_pool, size=N_PER_Q, replace=False)
            starts.sort()
            q64 = np.int64(q)
            for n in starts:
                n_int = int(n)
                n_logged, status, sigma = trajectory_qx1(
                    np.int64(n), q64, MAX_STEPS, MAX_VALUE_INT, log_buf
                )
                a_star, a_pow = a_star_for_n(n_int, K_MOD, q)
                status_str = STATUS_NAMES[status]
                status_counts[q][status_str] += 1
                if status == S_CONVERGED and sigma > sigma_max_seen:
                    sigma_max_seen = sigma

                # Vis A: trajectory points (subsample to VIS_SUBSAMPLE if longer)
                if n_logged > VIS_SUBSAMPLE:
                    idxs = np.linspace(0, n_logged - 1, VIS_SUBSAMPLE).astype(np.int64)
                else:
                    idxs = np.arange(n_logged, dtype=np.int64)
                for t_idx in idxs:
                    wa.writerow((
                        f"{log_buf[t_idx]:.6f}",
                        int(t_idx),
                        q,
                        a_star,
                        a_pow,
                        status_str,
                    ))
                    vis_a_count += 1

                # Vis B: per-orbit summary (sigma; sentinel -1 for divergent/timeout)
                log_n = float(np.log(n_int))
                wb.writerow((
                    f"{log_n:.6f}",
                    int(sigma) if status == S_CONVERGED else -1,
                    q,
                    a_star,
                    a_pow,
                    status_str,
                ))
                # Vis B': divergence-aware (log of final value before termination; 0 for converged)
                log_final = float(log_buf[n_logged - 1])
                wbd.writerow((
                    f"{log_n:.6f}",
                    f"{log_final:.6f}",
                    q,
                    a_star,
                    a_pow,
                    status_str,
                ))
                vis_b_count += 1

            elapsed = time.perf_counter() - t0
            sc = status_counts[q]
            print(
                f"[gen] q={q:>2}  conv={sc['converged']:>4}  div={sc['divergent']:>4}  "
                f"to={sc['timeout']:>4}  vis_a={vis_a_count:>9,}  vis_b={vis_b_count:>6,}  "
                f"t={elapsed:6.1f}s",
                flush=True,
            )

    print(f"\n[gen] DONE")
    print(f"[gen] vis_a.csv     rows: {vis_a_count:,}  ({VIS_A_CSV})")
    print(f"[gen] vis_b.csv     rows: {vis_b_count:,}  ({VIS_B_CSV})")
    print(f"[gen] vis_b_div.csv rows: {vis_b_count:,}  ({VIS_B_DIV_CSV})")
    print(f"[gen] max sigma observed (converged orbits): {sigma_max_seen}")
    print(f"[gen] status breakdown by q:")
    for q in Q_VALUES:
        sc = status_counts[q]
        total = sum(sc.values())
        print(
            f"       q={q:>2}: conv={sc['converged']:>4} ({sc['converged']/total:.1%})  "
            f"div={sc['divergent']:>4} ({sc['divergent']/total:.1%})  "
            f"to={sc['timeout']:>4} ({sc['timeout']/total:.1%})"
        )


if __name__ == "__main__":
    main()
