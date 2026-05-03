"""
esscher_duality.py — test the duality identity E_P[L⁻] = E_P*[L⁺*]
                    for the iid Syracuse log-walk under Esscher tilt at w*=1.

Walks:
  P:  v ~ Geom(1/2), step X = log 3 - v·log 2, mean = -log(4/3) < 0
  P*: v ~ Geom(3/4), step X = log 3 - v·log 2, mean = log(3·2^(-4/3)) > 0

Tests:
  (1) E_P[L⁻] (descending ladder mean, original)
  (2) E_P*[L⁺*] (ascending ladder mean, tilted)
  (3) Whether (1) == (2) — user's claimed duality
  (4) E_P[e^(-L⁻)] = q* prediction (Esscher relation)
  (5) Implied E*[L⁺*] from μ*/(1-q*) closed form

Cap < 200 lines.
"""
import numpy as np
import sys
from numba import njit

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2*LOG2 - LOG3      # = log(4/3) ≈ 0.2877
LOG32 = LOG3 - LOG2        # = log(3/2) ≈ 0.4055
MU_STAR = LOG3 - (4/3)*LOG2  # tilted drift, ≈ +0.1744


@njit(cache=True)
def geom_p(rng_state_seed, p, n):
    """Generate n iid Geom(p) on {1,2,...} via inverse CDF."""
    out = np.empty(n, dtype=np.int64)
    rng = np.random
    for i in range(n):
        u = np.random.random()
        out[i] = int(np.ceil(np.log(1 - u) / np.log(1 - p)))
    return out


@njit(cache=True)
def sim_descending(n_orbits, p, max_steps, seed):
    """Simulate strict descending ladder height L⁻ = -S(σ⁻).
    v ~ Geom(p) on {1,2,...}, X = log 3 - v·log 2."""
    np.random.seed(seed)
    L = np.zeros(n_orbits)
    failed = 0
    for i in range(n_orbits):
        S = 0.0
        for _ in range(max_steps):
            u = np.random.random()
            v = int(np.ceil(np.log(1 - u) / np.log(1 - p)))
            S += LOG3 - v * LOG2
            if S < 0:
                L[i] = -S
                break
        else:
            failed += 1
    return L, failed


@njit(cache=True)
def sim_ascending(n_orbits, p, max_steps, seed):
    """Simulate strict ascending ladder height L⁺ = S(σ⁺).
    v ~ Geom(p) on {1,2,...}."""
    np.random.seed(seed)
    L = np.zeros(n_orbits)
    failed = 0
    for i in range(n_orbits):
        S = 0.0
        for _ in range(max_steps):
            u = np.random.random()
            v = int(np.ceil(np.log(1 - u) / np.log(1 - p)))
            S += LOG3 - v * LOG2
            if S > 0:
                L[i] = S
                break
        else:
            failed += 1
    return L, failed


@njit(cache=True)
def sim_q_descending_under_P_star(n_orbits, p_star, max_steps, give_up, seed):
    """Under P* (positive drift), simulate to test if σ⁻* < ∞.
    Returns count of orbits that descend below 0 within give_up steps,
    OR until S > +threshold (then σ⁻* = ∞)."""
    np.random.seed(seed)
    descend_count = 0
    inf_count = 0
    for i in range(n_orbits):
        S = 0.0
        decided = False
        for _ in range(max_steps):
            u = np.random.random()
            v = int(np.ceil(np.log(1 - u) / np.log(1 - p_star)))
            S += LOG3 - v * LOG2
            if S < 0:
                descend_count += 1
                decided = True
                break
            if S > give_up:
                inf_count += 1
                decided = True
                break
        if not decided:
            inf_count += 1
    return descend_count, inf_count


def main():
    print("=" * 70)
    print("Esscher duality test: E_P[L⁻] vs E_P*[L⁺*]")
    print("=" * 70)
    print(f"P:  v ~ Geom(1/2), E[X] = -log(4/3) = {-LOG43:.5f}, NEGATIVE drift")
    print(f"P*: v ~ Geom(3/4), E*[X] = log(3·2^(-4/3)) = {MU_STAR:.5f}, POSITIVE drift")
    print(f"Cramer root w* = 1 (since 2² − 1 = 3)")

    n_orbits = 5_000_000
    print(f"\nSimulating {n_orbits:,} orbits each, numba-accelerated...", flush=True)

    # Original walk: descending ladder
    L_neg, fails_neg = sim_descending(n_orbits, p=0.5, max_steps=100_000, seed=42)
    if fails_neg:
        print(f"  WARN: {fails_neg} descending failures")
    E_Lneg = L_neg.mean()
    SE_neg = L_neg.std() / np.sqrt(n_orbits)
    print(f"\n  E_P[L⁻]   = {E_Lneg:.6f} ± {1.96*SE_neg:.6f}  (95% CI), nats")
    print(f"            = {E_Lneg/LOG43:.6f} step units")

    # Tilted walk: ascending ladder
    L_pos_star, fails_pos = sim_ascending(n_orbits, p=0.75, max_steps=100_000, seed=42)
    if fails_pos:
        print(f"  WARN: {fails_pos} ascending failures")
    E_Lpos_star = L_pos_star.mean()
    SE_pos = L_pos_star.std() / np.sqrt(n_orbits)
    print(f"  E_P*[L⁺*] = {E_Lpos_star:.6f} ± {1.96*SE_pos:.6f}  (95% CI), nats")
    print(f"            = {E_Lpos_star/LOG32:.6f} units of h = log(3/2)")

    # Duality test
    print()
    print("=" * 70)
    print("Duality identity test: E_P[L⁻] =? E_P*[L⁺*]")
    print("=" * 70)
    diff = E_Lneg - E_Lpos_star
    z = diff / np.sqrt(SE_neg**2 + SE_pos**2)
    print(f"  Difference  = {diff:+.6f} nats  ({z:+.2f}σ)")
    if abs(z) < 2:
        print(f"  ✓ Identity HOLDS within 2σ — duality confirmed empirically")
    else:
        print(f"  ✗ Identity FAILS at {abs(z):.2f}σ — duality as stated is INCORRECT")

    # Esscher-relation test: q* = E_P[e^(-L⁻)]
    print()
    print("=" * 70)
    print("Esscher relation: q* = E_P[e^(-L⁻)] (under tilt at w*=1)")
    print("=" * 70)
    E_exp_neg_L = np.exp(-L_neg).mean()
    print(f"  E_P[e^(-L⁻)]   = {E_exp_neg_L:.6f}  (= predicted q*)")
    print(f"  predicted 1-q* = {1 - E_exp_neg_L:.6f}")
    pred_E_Lpos_star = MU_STAR / (1 - E_exp_neg_L)
    print(f"  predicted E*[L⁺*] = μ*/(1-q*) = {pred_E_Lpos_star:.6f}")
    print(f"  empirical E*[L⁺*] = {E_Lpos_star:.6f}")
    print(f"  diff               = {pred_E_Lpos_star - E_Lpos_star:+.6f}")

    # Direct measurement of q* under tilted measure
    print()
    print("Direct measurement of q* under P* (gives ground truth)...")
    descend_count, inf_count = sim_q_descending_under_P_star(
        n_orbits, p_star=0.75, max_steps=100_000, give_up=40.0, seed=137)
    q_star_direct = descend_count / n_orbits
    print(f"  q* (direct)   = {q_star_direct:.6f}")
    print(f"  q* (Esscher)  = E_P[e^(-L⁻)] = {E_exp_neg_L:.6f}")
    diff_q = q_star_direct - E_exp_neg_L
    SE_q_dir = np.sqrt(q_star_direct * (1 - q_star_direct) / n_orbits)
    print(f"  diff = {diff_q:+.6f} ({diff_q/SE_q_dir:+.2f}σ)")

    # Implied E[L⁻] = μ/(1-q) self-consistency
    # We need q from original measure
    print()
    print("=" * 70)
    print("Self-consistency: E[L⁻] = μ/(1-q) for original walk")
    print("=" * 70)
    # 1 - q = μ/E[L⁻] (from sim)
    one_minus_q = LOG43 / E_Lneg
    q_implied = 1 - one_minus_q
    print(f"  1 - q (implied) = μ/E[L⁻] = {LOG43:.4f}/{E_Lneg:.4f} = {one_minus_q:.6f}")
    print(f"  q (implied)     = {q_implied:.6f}")

    print()
    print("=" * 70)
    print("Summary of structural identities")
    print("=" * 70)
    print(f"  E_P[L⁻]              = {E_Lneg:.6f}  (sim, original walk)")
    print(f"  E_P*[L⁺*]            = {E_Lpos_star:.6f}  (sim, tilted walk)")
    print(f"  μ/(1-q) self-check   = {LOG43:.6f}/{one_minus_q:.6f} = {LOG43/one_minus_q:.6f}")
    print(f"  μ*/(1-q*) self-check = {MU_STAR:.6f}/{1-q_star_direct:.6f} = {MU_STAR/(1-q_star_direct):.6f}")
    print(f"  q* = E[e^(-L⁻)] check: {q_star_direct:.6f} vs {E_exp_neg_L:.6f}")


if __name__ == "__main__":
    main()
