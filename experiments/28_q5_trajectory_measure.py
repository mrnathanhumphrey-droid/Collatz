"""
Experiment 28 — q=5 trajectory v-distribution on convergent orbits.

The qx+1 trajectory measure on v_2(qm+1) is i.i.d. Geom(1/2) at q ∈ {5, 7, 9, 11}
when sampled UNCONDITIONALLY across all orbits (mostly divergent at q≥5)
— verified in findings.md 2026-05-02.

This experiment asks the *conditional* question: do convergent orbits at q=5
also have Geom(1/2) v's, or is the v-distribution tilted on the rare-event
subset? Cramér large deviations theory predicts that convergent orbits sit on
the upper-v-biased tilted measure (more v=high steps drag iterates down faster).

Method:
  - Load convergent q=5 orbits from data/q_main_q5_N100000000.parquet
  - For each, simulate Syracuse trajectory m → (5m+1)/2^v until m == 1
  - Pool all v's, histogram, compare to Geom(1/2) and to 3x+1 traj from exp 25

Usage:
    python 28_q5_trajectory_measure.py
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")

MAX_V = 35


def collect_v_q5(starts):
    """Pure-Python (int) Syracuse walk at q=5 until m == 1 or overflow.
    starts: list of int (odd, convergent at q=5).
    Returns: histogram array, total v-samples, max trajectory length seen.
    """
    hist = np.zeros(MAX_V + 2, dtype=np.int64)
    total = 0
    max_len = 0
    not_converged = 0
    for n in starts:
        m = int(n)
        steps = 0
        # safety cap: convergent orbits at q=5 N=10^8 have max odd_steps = 265
        while m != 1 and steps < 5000:
            tmp = 5 * m + 1
            v = 0
            while tmp % 2 == 0:
                tmp //= 2
                v += 1
            if v <= MAX_V:
                hist[v] += 1
            else:
                hist[MAX_V + 1] += 1
            total += 1
            m = tmp
            steps += 1
        if m != 1:
            not_converged += 1
        max_len = max(max_len, steps)
    return hist, total, max_len, not_converged


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[load] q=5 N=10^8 parquet", flush=True)
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / "q_main_q5_N100000000.parquet")
    conv = df.filter(pl.col("converged"))
    starts = conv["n"].to_list()
    print(f"        {len(starts):,} convergent orbits  ({time.perf_counter()-t0:.1f}s)",
          flush=True)

    print("[run] q=5 trajectory v collection (Python int Syracuse) ...", flush=True)
    t0 = time.perf_counter()
    hist, total, max_len, not_conv = collect_v_q5(starts)
    print(f"        done in {time.perf_counter()-t0:.1f}s   "
          f"v-samples: {total:,}   max traj len: {max_len}   "
          f"didn't converge: {not_conv}",
          flush=True)
    print()

    # Load q=3 traj for comparison
    q3 = pl.read_csv(out_dir / "25_trajectory_measure_Nstart100000000_T200.csv")

    # Compare
    print(f"{'v':>3} {'2^-v Geom':>11} {'q5 count':>11} {'q5 ratio':>10} "
          f"{'q3 ratio (exp25)':>17} {'q5 SE':>9}")
    for v in range(1, MAX_V + 1):
        cq5 = int(hist[v])
        if total == 0:
            break
        emp_q5 = cq5 / total
        pred = 2.0 ** (-v)
        ratio_q5 = emp_q5 / pred if pred > 0 else float("nan")
        # Get q=3 ratio for same v
        q3_row = q3.filter(pl.col("v") == v)
        ratio_q3 = float(q3_row["traj_ratio"][0]) if len(q3_row) else float("nan")
        # SE: ~ 1/sqrt(expected_count) = 1/sqrt(total * 2^-v)
        expected = total * pred
        se = 1.0 / np.sqrt(expected) if expected > 0 else float("nan")
        print(f"{v:>3} {pred:>11.6f} {cq5:>11,} {ratio_q5:>10.4f} "
              f"{ratio_q3:>17.4f} {se:>9.4f}")

    # Moments
    v_arr = np.arange(1, MAX_V + 1, dtype=np.float64)
    p_q5 = hist[1:MAX_V + 1].astype(np.float64) / total
    E_v = (v_arr * p_q5).sum()
    E_v2 = (v_arr ** 2 * p_q5).sum()
    Var_v = E_v2 - E_v ** 2
    print()
    print(f"q=5 convergent-trajectory v-distribution moments:")
    print(f"  E[v]   = {E_v:.6f}    (Geom(1/2): 2.000000)   delta = {E_v-2:+.6f}  ({100*(E_v-2)/2:+.4f}%)")
    print(f"  Var[v] = {Var_v:.6f}  (Geom(1/2): 2.000000)   delta = {Var_v-2:+.6f}  ({100*(Var_v-2)/2:+.4f}%)")

    # MGF probe: Cramer rate uses E[(2/q)^v]; under Geom(1/2) this is 1/(q-1)
    print()
    print(f"q=5 MGF probe E[(2/q)^v] (under Geom(1/2): 1/(q-1)):")
    for q in [3, 5, 7, 9, 11]:
        mgf = (((2.0 / q) ** v_arr) * p_q5).sum()
        pred = 1.0 / (q - 1)
        print(f"  q={q:2d}: E[(2/q)^v] = {mgf:.6f}   1/(q-1) = {pred:.6f}   "
              f"delta = {mgf-pred:+.6f}  ({100*(mgf-pred)/pred:+.4f}%)")

    # Save CSV
    rows = []
    for v in range(1, MAX_V + 1):
        cq5 = int(hist[v])
        emp_q5 = cq5 / total if total > 0 else 0
        pred = 2.0 ** (-v)
        ratio_q5 = emp_q5 / pred if pred > 0 else float("nan")
        rows.append({"v": v, "geom_pred": pred, "q5_count": cq5,
                     "q5_ratio": ratio_q5})
    out = pl.DataFrame(rows)
    out_csv = out_dir / "28_q5_trajectory_v_convergent.csv"
    out.write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
