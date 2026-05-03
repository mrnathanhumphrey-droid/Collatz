"""
wh_numerical_check.py — Path A numerical verification of the iid-Geom(1/2)
Wald-Lorden baseline for the Syracuse log-walk's descending ladder.

Computes:
  1. Step moments E[X], E[X²], Var[X] from iid Geom(1/2) approximation.
  2. Direct simulation of strict descending ladder height L⁻ over 10⁶ orbits.
  3. Comparison against the analytic candidate E[X²]/(2|E[X]|) (Lorden).
  4. Gap to empirical W_2 = 7.156 from compute_threads_findings.md.

The walk: X_t = log(3) − v_t · log(2), v_t ~ iid Geom(1/2) on {1,2,3,...}.
Mean E[X] = log(3/4) < 0 (descending). The strict descending ladder time
σ⁻ = inf{n ≥ 1 : S_n < 0} is a.s. finite; L⁻ = −S(σ⁻) > 0 a.s.

Keep < 200 lines per brief.
"""
import numpy as np
import sys

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3   # log(4/3)


def step_moments():
    """Compute E[X], E[X²], Var[X], E[X⁻] for X = log(3) − v log(2), v ~ Geom(1/2)."""
    E_v = 2.0
    E_v2 = 6.0     # E[v²] for Geom(1/2), v ∈ {1,2,...}
    Var_v = E_v2 - E_v ** 2
    E_X = LOG3 - E_v * LOG2
    E_X2 = LOG3 ** 2 - 2 * LOG3 * LOG2 * E_v + LOG2 ** 2 * E_v2
    Var_X = E_X2 - E_X ** 2
    # E[X · 1_{X<0}] = sum_{k≥2} (log 3 - k log 2) · 2^{-k}
    # = log 3 · (1/2) - log 2 · (3/2) = log(3/8)/2
    E_X_times_neg = (LOG3 - 3 * LOG2) / 2.0
    E_X_neg = -E_X_times_neg            # E[max(-X, 0)] = E[-X · 1_{X<0}]
    return {
        "E_v": E_v, "Var_v": Var_v,
        "E_X": E_X, "E_X2": E_X2, "Var_X": Var_X,
        "E_X_neg": E_X_neg, "mu": -E_X,
    }


def simulate_Lneg(n_orbits=10**6, seed=42, max_steps=20_000):
    """Direct simulation of strict descending ladder height L⁻ for the
    Syracuse log-walk under iid Geom(1/2) approximation."""
    rng = np.random.default_rng(seed)
    L_vals = np.zeros(n_orbits)
    failed = 0
    for i in range(n_orbits):
        S = 0.0
        crossed = False
        for _ in range(max_steps):
            v = rng.geometric(0.5)        # support {1, 2, 3, ...}
            S += LOG3 - v * LOG2
            if S < 0:
                L_vals[i] = -S
                crossed = True
                break
        if not crossed:
            failed += 1
    return L_vals[:n_orbits], failed


def main():
    m = step_moments()
    print("=" * 64)
    print("Step moments (iid Geom(1/2) approximation for v_t)")
    print("=" * 64)
    print(f"  E[v]       = {m['E_v']:.6f}")
    print(f"  Var[v]     = {m['Var_v']:.6f}")
    print(f"  E[X]       = {m['E_X']:+.6f}    (= -log(4/3) = {-LOG43:+.6f})")
    print(f"  E[X²]      = {m['E_X2']:.6f}")
    print(f"  Var[X]     = {m['Var_X']:.6f}")
    print(f"  E[X⁻]      = {m['E_X_neg']:.6f}    (mean magnitude of negative jumps)")
    print(f"  log(4/3)   = {LOG43:.6f}")

    # Analytic candidates for E[L⁻]
    cand_X2 = m["E_X2"] / (2 * m["mu"])           # Lorden upper bound / asymptotic residual
    cand_var = m["Var_X"] / (2 * m["mu"])         # alternative form
    cand_Xneg = m["E_X_neg"]                      # mean of single-step negative jump
    print()
    print("=" * 64)
    print("Analytic candidate values (closed forms to compare)")
    print("=" * 64)
    print(f"  E[X²]/(2|E[X]|)   = {cand_X2:.6f} nats = {cand_X2/LOG43:.6f} steps  (Lorden formula)")
    print(f"  Var[X]/(2|E[X]|)  = {cand_var:.6f} nats = {cand_var/LOG43:.6f} steps")
    print(f"  E[X⁻]              = {cand_Xneg:.6f} nats = {cand_Xneg/LOG43:.6f} steps")

    # Simulation
    print()
    print("=" * 64)
    print("Direct simulation of L⁻")
    print("=" * 64)
    n_orbits = 1_000_000
    print(f"  Simulating {n_orbits:,} iid orbits...", flush=True)
    L_vals, failed = simulate_Lneg(n_orbits=n_orbits, seed=42)
    print(f"  Failed (no crossing in max_steps): {failed}")
    E_Lneg = L_vals.mean()
    SE = L_vals.std() / np.sqrt(len(L_vals))
    Var_Lneg = L_vals.var()
    print(f"  E[L⁻]_sim   = {E_Lneg:.6f} ± {1.96*SE:.4f} (95% CI) nats")
    print(f"             = {E_Lneg/LOG43:.6f} steps")
    print(f"  Var[L⁻]_sim = {Var_Lneg:.6f}")
    # Higher moments and median for sanity
    print(f"  Median[L⁻]  = {np.median(L_vals):.6f} nats")

    # Compare to candidates
    print()
    print("=" * 64)
    print("Match to analytic candidates")
    print("=" * 64)
    for name, c_nats in [
        ("E[X²]/(2|E[X]|)", cand_X2),
        ("Var[X]/(2|E[X]|)", cand_var),
        ("E[X⁻]", cand_Xneg),
    ]:
        diff = E_Lneg - c_nats
        z = diff / SE
        print(f"  {name:>20}: candidate {c_nats:.4f} nats vs simulated {E_Lneg:.4f}  "
              f"(diff = {diff:+.4f}, z = {z:+.2f})")

    # Comparison to empirical W_j
    W2_emp_steps = 7.156
    W2_emp_nats = W2_emp_steps * LOG43
    print()
    print("=" * 64)
    print("Comparison to empirical W_j (compute_threads_findings.md)")
    print("=" * 64)
    print(f"  W_iid (simulated)    = {E_Lneg/LOG43:.4f} steps  ({E_Lneg:.4f} nats)")
    print(f"  W_2 empirical         = {W2_emp_steps:.4f} steps  ({W2_emp_nats:.4f} nats)")
    print(f"  W_4 empirical         = -4.755 steps")
    print(f"  W_5 empirical         = +4.590 steps")
    gap_steps = W2_emp_steps - E_Lneg / LOG43
    gap_nats = W2_emp_nats - E_Lneg
    pct = 100 * gap_steps / (E_Lneg / LOG43)
    print(f"\n  Gap (W_2 - W_iid)    = {gap_steps:+.4f} steps  ({gap_nats:+.4f} nats)")
    print(f"  Relative              = {pct:+.2f}%  of W_iid")
    print(f"\n  iid baseline is universal across j (single value W_iid).")
    print(f"  Empirical W_j varies by j (W_2 ≠ W_4 ≠ W_5);")
    print(f"  per-j cross-class structure CANNOT be derived from iid Wiener-Hopf.")
    print(f"  → Path B (Markov-modulated, matrix-WH) needed for j-stratification.")


if __name__ == "__main__":
    main()
