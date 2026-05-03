"""
Experiment 15 — Trajectory step-variance measurement for qx+1.

Random-walk derivation predicts log(conv_rate(j)) slope = -2*log(q/4)*log(q)/V(q)
where V(q) is the variance of log(per-Syracuse-step change) under the trajectory
measure. Each Syracuse step takes odd x to (q*x+1)/2^v with v = nu_2(q*x+1) >= 1.
Step log-change: Delta = log(q) - v*log(2) (approximately). Variance:
    V(q) = Var(v) * log(2)^2

Geom(1/2) prediction (i.i.d. v_2 ~ Geom(1/2)): Var(v) = 2, V_geom = 2*log(2)^2 ~= 0.961.

If the empirical j-slope is -C*log(q/4), gambler's ruin gives C = 2*log(q)/V(q),
so V(q) = 2*log(q)/C.

Empirical slope at q=5 is -0.5619 over j increments. C_empirical = -slope/log(q/4) =
2.518. Predicted V(5) under C=5/2: V(5) = 2*log(5)/(5/2) = 1.288.
Predicted under Geom(1/2) baseline: V(5) = 0.961.

Measure empirical V(5) by replaying convergent qx+1 trajectories at q=5 and
recording v_2 per Syracuse step.

Sample size: 1000 convergent q=5 orbits, replayed step-by-step.

Usage:
    python 15_step_variance.py --q 5 --N 100000000 --sample 1000
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")


def replay_trajectory(n0, q, max_steps=100000):
    """Replay qx+1 trajectory from n0; record v_2 per Syracuse step.
    Returns list of v_2 values, one per Syracuse step."""
    x = int(n0)
    v_seq = []
    steps = 0
    while x != 1 and steps < max_steps:
        # Standard map: if even, halve; if odd, q*x+1.
        # Syracuse step: starting at odd x, do q*x+1 then strip all 2's.
        if x % 2 == 0:
            x = x // 2
            steps += 1
            continue
        # x is odd
        x = q * x + 1
        steps += 1
        v = 0
        while x % 2 == 0:
            x //= 2
            v += 1
            steps += 1
        if v == 0:
            # shouldn't happen; q*x+1 with x odd is even
            continue
        v_seq.append(v)
    return v_seq, x == 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--N", type=int, default=100_000_000)
    ap.add_argument("--sample", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    print(f"[load] q={args.q}, N={args.N:,}, sample={args.sample}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet").filter(pl.col("converged"))
    n_conv = len(df)
    print(f"        {n_conv:,} convergent orbits in data", flush=True)
    if n_conv == 0:
        print("[error] no convergent orbits to sample")
        return

    rng = np.random.default_rng(args.seed)
    sample_size = min(args.sample, n_conv)
    idx = rng.choice(n_conv, sample_size, replace=False)
    ns = df["n"].to_numpy()[idx]

    print(f"[replay] sampling {sample_size} convergent orbits ...", flush=True)
    all_vs = []
    n_orbits_done = 0
    for n in ns:
        vs, ok = replay_trajectory(int(n), args.q)
        if ok and len(vs) > 0:
            all_vs.extend(vs)
            n_orbits_done += 1
    print(f"        {n_orbits_done} orbits replayed; {len(all_vs)} Syracuse steps total", flush=True)

    arr = np.array(all_vs, dtype=np.float64)
    v_mean = arr.mean()
    v_var = arr.var(ddof=1)
    v_min, v_max = arr.min(), arr.max()

    log2 = np.log(2.0)
    V_emp = v_var * log2 ** 2
    V_geom = 2.0 * log2 ** 2  # Geom(1/2) variance is 2; scale by log(2)^2

    # Predicted V(q) under empirical-slope C = 2.518 (for q=5)
    # C_emp = 2*log(q)/V(q)  =>  V(q) = 2*log(q)/C_emp
    if args.q == 5:
        C_emp = 2.518
    elif args.q == 7:
        C_emp = 2.445
    elif args.q == 9:
        C_emp = 2.531
    elif args.q == 11:
        C_emp = 1.628
    else:
        C_emp = 2.5
    V_pred_from_C = 2.0 * np.log(args.q) / C_emp

    # Cross-check: under C = 5/2 exactly
    V_pred_C25 = 2.0 * np.log(args.q) / 2.5

    # Drift mu = log(q) - E[v]*log(2) = log(q/4) under E[v]=2
    drift_emp = np.log(args.q) - v_mean * log2
    drift_geom = np.log(args.q / 4.0)

    print()
    print(f"=== Trajectory v_2 distribution (q={args.q}, n_orbits={n_orbits_done}, n_steps={len(all_vs):,}) ===")
    print(f"  Empirical v_2 distribution:")
    print(f"    mean  = {v_mean:.4f}    (Geom(1/2) prediction: 2.0)")
    print(f"    var   = {v_var:.4f}    (Geom(1/2) prediction: 2.0)")
    print(f"    range = [{int(v_min)}, {int(v_max)}]")
    print(f"  Per-step v_2 histogram (top 10):")
    counts = np.bincount(arr.astype(int))[:11]
    total = counts.sum()
    geom_pmf = np.array([0.5 ** v for v in range(1, 11)])  # P(v=k) = 0.5^k for Geom(1/2) on {1,2,...}
    geom_pmf = geom_pmf / geom_pmf.sum()
    for v in range(1, min(11, len(counts))):
        emp_p = counts[v] / total if total > 0 else 0
        geom_p = geom_pmf[v - 1]
        print(f"    v={v:>2}: emp={emp_p:.4f}  geom={geom_p:.4f}  ratio={emp_p / geom_p:.3f}")

    print()
    print(f"=== V(q) = Var(v)*log(2)^2 — the trajectory step variance ===")
    print(f"  V_empirical                       = {V_emp:.4f}")
    print(f"  V_Geom(1/2) (i.i.d. baseline)     = {V_geom:.4f}")
    print(f"  V_predicted from empirical C={C_emp:.3f} = {V_pred_from_C:.4f}")
    print(f"  V_predicted under C=5/2 exactly    = {V_pred_C25:.4f}")

    print()
    print(f"=== Drift check: log(q) - E[v]*log(2) ===")
    print(f"  Empirical drift = {drift_emp:.4f}")
    print(f"  Predicted drift = log(q/4) = {drift_geom:.4f}")
    print(f"  Ratio emp/pred  = {drift_emp / drift_geom:.4f}")

    print()
    print(f"=== Verdict ===")
    if abs(V_emp - V_pred_from_C) / V_pred_from_C < 0.10:
        print(f"  V_emp matches V_predicted_from_C within 10%. Trajectory step variance is")
        print(f"  the right mechanism. The random-walk derivation works once we use the")
        print(f"  empirical step variance (rather than Geom(1/2) baseline).")
    elif abs(V_emp - V_geom) / V_geom < 0.10:
        print(f"  V_emp matches Geom(1/2) baseline. Trajectory step variance is NOT the")
        print(f"  source of the empirical correction. The C != log(q)/log(2)^2 deviation")
        print(f"  has a different mechanism (correlations, finite-N selection, cycle mixing).")
    else:
        print(f"  V_emp = {V_emp:.4f} matches neither prediction cleanly.")
        print(f"  Either the random-walk framework is too simplified for qx+1, or there's")
        print(f"  a different correction factor we haven't identified.")


if __name__ == "__main__":
    main()
