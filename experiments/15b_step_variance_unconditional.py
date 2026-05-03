"""
Experiment 15b — UNCONDITIONAL trajectory step-variance measurement.

Sample random odd n in [1, N] (no convergence filter), run qx+1 trajectory
for up to max_steps Syracuse steps, record v_2 per step. Pool across orbits.

This gives the unconditional v_2 distribution, which is what gambler's ruin
uses to predict conv_rate. Earlier (15) sampled CONVERGENT orbits only —
that's biased toward large v (the survival bias).

Comparison targets at q=5:
  Geom(1/2) baseline:        Var(v)=2,  V=2*log(2)^2 = 0.961
  Under empirical C=2.518:    V = 2*log(5)/2.518 = 1.278
  Under C=5/2 exactly:        V = 2*log(5)/2.5  = 1.288

Usage:
    python 15b_step_variance_unconditional.py --q 5 --sample 5000 --max_steps 200
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


def replay_to_max_syracuse_steps(n0, q, max_syr_steps, max_value):
    """Run qx+1 trajectory from odd n0, recording v_2 per Syracuse step.
    Stop after max_syr_steps Syracuse steps OR if x exceeds max_value (overflow guard)
    OR if x reaches 1.
    Returns list of v_2 values, one per Syracuse step taken."""
    x = int(n0)
    v_seq = []
    syr_steps = 0
    while syr_steps < max_syr_steps and x != 1 and x < max_value:
        # x is odd at start of each Syracuse step (we only enter with odd x; pre-strip if needed)
        if x % 2 == 0:
            # strip evens (shouldn't happen with odd seed unless mid-trajectory)
            while x % 2 == 0:
                x //= 2
            continue
        x = q * x + 1
        v = 0
        while x % 2 == 0:
            x //= 2
            v += 1
        v_seq.append(v)
        syr_steps += 1
    return v_seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--N_max", type=int, default=10_000_000,
                    help="Sample odd n uniformly in [1, N_max].")
    ap.add_argument("--sample", type=int, default=5000)
    ap.add_argument("--max_steps", type=int, default=200,
                    help="Syracuse steps per orbit (truncate divergent ones).")
    ap.add_argument("--max_value", type=float, default=1e80,
                    help="Overflow guard; if x exceeds this, stop.")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    # Sample random odd n
    odd_candidates = rng.integers(1, args.N_max // 2, size=args.sample) * 2 + 1
    odd_candidates = np.unique(odd_candidates)
    if len(odd_candidates) < args.sample:
        # top up
        more = rng.integers(1, args.N_max // 2, size=args.sample - len(odd_candidates)) * 2 + 1
        odd_candidates = np.unique(np.concatenate([odd_candidates, more]))
    odd_candidates = odd_candidates[:args.sample]

    print(f"[setup] q={args.q}, sample={len(odd_candidates)} random odd n in [1, {args.N_max:,}]", flush=True)
    print(f"        max_syr_steps={args.max_steps}, max_value={args.max_value:.2e}", flush=True)

    all_vs = []
    n_truncated = 0
    n_converged = 0
    n_overflow = 0
    max_value = int(args.max_value)
    for n in odd_candidates:
        vs = replay_to_max_syracuse_steps(int(n), args.q, args.max_steps, max_value)
        all_vs.extend(vs)
        # Classify ending
        if len(vs) == args.max_steps:
            n_truncated += 1
        else:
            # ended either by x==1 or x>max_value, replay to find which
            x = int(n)
            for v in vs:
                if x % 2 == 0:
                    while x % 2 == 0: x //= 2
                x = args.q * x + 1
                while x % 2 == 0: x //= 2
            if x == 1:
                n_converged += 1
            else:
                n_overflow += 1
    print(f"[replay] orbits done: converged={n_converged}, truncated_at_max_steps={n_truncated}, overflow={n_overflow}", flush=True)
    print(f"        total Syracuse steps recorded: {len(all_vs):,}", flush=True)

    arr = np.array(all_vs, dtype=np.float64)
    if len(arr) == 0:
        print("[error] no steps recorded")
        return

    v_mean = arr.mean()
    v_var = arr.var(ddof=1)
    v_min, v_max = arr.min(), arr.max()

    log2 = np.log(2.0)
    V_emp = v_var * log2 ** 2
    V_geom = 2.0 * log2 ** 2

    if args.q == 5: C_emp = 2.518
    elif args.q == 7: C_emp = 2.445
    elif args.q == 9: C_emp = 2.531
    elif args.q == 11: C_emp = 1.628
    else: C_emp = 2.5
    V_pred_from_C = 2.0 * np.log(args.q) / C_emp
    V_pred_C25 = 2.0 * np.log(args.q) / 2.5

    drift_emp = np.log(args.q) - v_mean * log2
    drift_pred = np.log(args.q / 4.0)

    print()
    print(f"=== UNCONDITIONAL v_2 distribution (q={args.q}, n_steps={len(all_vs):,}) ===")
    print(f"  mean  = {v_mean:.4f}    (Geom(1/2) prediction: 2.0)")
    print(f"  var   = {v_var:.4f}    (Geom(1/2) prediction: 2.0)")
    print(f"  range = [{int(v_min)}, {int(v_max)}]")
    print(f"  v_2 PMF (first 10):")
    counts = np.bincount(arr.astype(int))[:11]
    total = counts.sum()
    geom_pmf = np.array([0.5 ** v for v in range(1, 11)])
    geom_pmf = geom_pmf / geom_pmf.sum()
    for v in range(1, min(11, len(counts))):
        emp_p = counts[v] / total if total > 0 else 0
        geom_p = geom_pmf[v - 1]
        ratio = emp_p / geom_p if geom_p > 0 else float("inf")
        print(f"    v={v:>2}: emp={emp_p:.4f}  geom={geom_p:.4f}  ratio={ratio:.3f}")

    print()
    print(f"=== V(q) = Var(v) * log(2)^2 ===")
    print(f"  V_empirical (unconditional)        = {V_emp:.4f}")
    print(f"  V_Geom(1/2) (i.i.d. baseline)      = {V_geom:.4f}")
    print(f"  V_predicted from empirical C={C_emp:.3f}  = {V_pred_from_C:.4f}")
    print(f"  V_predicted under C=5/2 exactly    = {V_pred_C25:.4f}")

    print()
    print(f"=== Drift check ===")
    print(f"  Empirical drift = log(q) - E[v]*log(2) = {drift_emp:.4f}")
    print(f"  Predicted drift = log(q/4)             = {drift_pred:.4f}")

    print()
    print(f"=== Verdict ===")
    err_C = abs(V_emp - V_pred_from_C) / V_pred_from_C
    err_geom = abs(V_emp - V_geom) / V_geom
    if err_C < 0.10:
        print(f"  V_emp = {V_emp:.4f} matches V_pred_from_C = {V_pred_from_C:.4f} within 10%.")
        print(f"  --> The empirical C = 2.518 is consistent with the gambler's-ruin formula")
        print(f"      using the UNCONDITIONAL v_2 step variance. Trajectory measure is the mechanism.")
    elif err_geom < 0.10:
        print(f"  V_emp = {V_emp:.4f} matches V_Geom(1/2) = {V_geom:.4f} within 10%.")
        print(f"  --> The trajectory measure IS approximately Geom(1/2). The empirical C deviation")
        print(f"      from random-walk prediction must come from something OTHER than step variance")
        print(f"      (correlations, finite-time effects, boundary behavior).")
    else:
        print(f"  V_emp = {V_emp:.4f} matches neither cleanly.")
        print(f"  Geom(1/2) baseline: {V_geom:.4f} (off by {100*err_geom:.1f}%)")
        print(f"  C-based prediction:  {V_pred_from_C:.4f} (off by {100*err_C:.1f}%)")


if __name__ == "__main__":
    main()
