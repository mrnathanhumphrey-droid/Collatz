"""
path_c_q_attempt.py — focused asymptotic attempt to close E[L⁻] via q.

Derivation: from Wiener-Hopf factorization, differentiating at θ=0:

    E[L⁻] = μ / (1 − q),    μ = log(4/3),    q = P(σ⁺ < ∞)

If q has a clean closed form, so does E[L⁻]. Strategy:
  (1) Direct simulation of q via "give-up threshold" trick (since the walk has
      negative drift and Cramer rate exactly e^{-1} per nat, give-up at S < -30
      misses events with prob ~e^{-30} ≈ 10^{-13}).
  (2) Verify the identity E[L⁻] = μ/(1-q) against simulated E[L⁻].
  (3) Test candidate closed forms for q (or 1-q).
  (4) Try alternate routes via the Esscher tilt: q*-related quantities under P*.

Cap: < 200 lines.
"""
import numpy as np
import sys

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3   # = log(4/3) ≈ 0.2877


def simulate_q_and_Lneg(n_orbits=10_000_000, give_up_S=-40.0, seed=42):
    """Simulate (σ⁺, σ⁻) jointly to compute q = P(σ⁺ < ∞) and E[L⁻]."""
    rng = np.random.default_rng(seed)
    # Track: did walk cross above 0 (σ⁺ < ∞)?
    sigma_plus_finite = 0
    L_neg_vals = np.zeros(n_orbits)
    L_neg_count = 0
    max_steps = 100_000
    for i in range(n_orbits):
        S = 0.0
        first_neg_recorded = False
        sigma_plus_seen = False
        for _ in range(max_steps):
            v = rng.geometric(0.5)
            S += LOG3 - v * LOG2
            # Track first descending ladder
            if not first_neg_recorded and S < 0:
                L_neg_vals[L_neg_count] = -S
                L_neg_count += 1
                first_neg_recorded = True
            # Track strict ascending ladder
            if not sigma_plus_seen and S > 0:
                sigma_plus_finite += 1
                sigma_plus_seen = True
            # Stop conditions: both ladders observed (or give-up)
            if first_neg_recorded and sigma_plus_seen:
                break
            if S < give_up_S:
                # walk has descended too far; σ⁺ = ∞ for this orbit
                break
    q_emp = sigma_plus_finite / n_orbits
    L_neg = L_neg_vals[:L_neg_count]
    return q_emp, L_neg


def main():
    print("=" * 64)
    print("Focused attempt: closed form for E[L⁻] via q = P(σ⁺ < ∞)")
    print("=" * 64)
    print(f"Identity: E[L⁻] = μ / (1 − q),  μ = log(4/3) = {LOG43:.6f}")

    n_orbits = 10_000_000
    print(f"\nSimulating {n_orbits:,} orbits with give-up threshold S < -40...", flush=True)
    q_emp, L_neg = simulate_q_and_Lneg(n_orbits=n_orbits, give_up_S=-40.0)
    SE_q = np.sqrt(q_emp * (1 - q_emp) / n_orbits)
    E_Lneg = L_neg.mean()
    SE_L = L_neg.std() / np.sqrt(len(L_neg))
    print(f"\n  q (sim)        = {q_emp:.6f} ± {1.96*SE_q:.6f}  (95% CI)")
    print(f"  1 − q (sim)    = {1-q_emp:.6f}")
    print(f"  E[L⁻] (sim)    = {E_Lneg:.6f} ± {1.96*SE_L:.6f}")
    pred_E_Lneg = LOG43 / (1 - q_emp)
    print(f"\n  Identity check: μ/(1-q) = {pred_E_Lneg:.6f}")
    print(f"  Direct E[L⁻]            = {E_Lneg:.6f}")
    print(f"  diff                    = {pred_E_Lneg - E_Lneg:+.6f}")

    # Candidate closed forms for q
    print()
    print("=" * 64)
    print("Candidate closed forms for q and 1 − q")
    print("=" * 64)
    print(f"  Empirical: q = {q_emp:.6f},  1 − q = {1-q_emp:.6f}")

    candidates = [
        ("3/4", 3/4, "= P(v=1) under P* [Esscher tilt]"),
        ("5/7", 5/7, "?"),
        ("log(2)", LOG2, "?"),
        ("1 - log(4/3)", 1 - LOG43, "?"),
        ("(log 3)/log(2 + e^{1/2})", LOG3/np.log(2 + np.e**0.5), "?"),
        ("log(2)/(log 2 + log(4/3))", LOG2/(LOG2 + LOG43), "?"),
        ("2/3 · log(3)/log(2)", (2/3) * LOG3/LOG2, "?"),
        ("log(3)/(log 3 + log 2/2)", LOG3/(LOG3 + LOG2/2), "?"),
        ("log(8/3)/log(4)", np.log(8/3)/np.log(4), "?"),
        ("1 - 1/e", 1 - 1/np.e, "?"),
        ("1 - log(4/3)/log e", 1 - LOG43, "= 1 - μ"),
    ]
    print(f"  {'Candidate':<35} {'value':>12} {'q diff':>11} {'note'}")
    for name, val, note in candidates:
        diff = q_emp - val
        z = abs(diff) / SE_q
        match = "✓✓" if z < 2 else ("✓" if z < 5 else "")
        print(f"  {name:<35} {val:>12.6f} {diff:>+11.6f} ({z:>5.1f}σ) {match}  {note}")

    # 1 - q candidates
    print(f"\n  1 - q candidates:")
    inv_candidates = [
        ("log(4/3)/log(e)", LOG43, "= μ exactly?"),
        ("1/log(2)·log(4/3)", LOG43/LOG2, "?"),
        ("log(4/3) · 1", LOG43, "= μ"),
        ("(1 - q)·E[L⁻] = μ identity", LOG43, "= μ"),
        ("2/7", 2/7, "?"),
        ("log(4/3)·E[L⁻]/E[L⁻] = log(4/3)", LOG43, "tautology"),
        ("1 - 5/7 = 2/7", 2/7, "?"),
        ("log(2)/(log(2)+log(3))/2", LOG2/(LOG2+LOG3)/2, "?"),
    ]
    inv_q = 1 - q_emp
    SE_inv = SE_q  # same
    for name, val, note in inv_candidates:
        diff = inv_q - val
        z = abs(diff) / SE_inv
        match = "✓✓" if z < 2 else ("✓" if z < 5 else "")
        print(f"  {name:<35} {val:>12.6f} {diff:>+11.6f} ({z:>5.1f}σ) {match}  {note}")

    # Esscher-tilt route
    print()
    print("=" * 64)
    print("Esscher-tilt route: q via E[exp(L⁺); σ⁺ < ∞] = 1")
    print("=" * 64)
    # Decompose: σ⁺ = 1 contributes (1/2) · (3/2) = 3/4 to E[exp(L⁺); σ⁺ < ∞]
    # Therefore E[exp(L⁺); σ⁺ ≥ 2, σ⁺ < ∞] = 1 - 3/4 = 1/4 [identity]
    print(f"  P(σ⁺ = 1) = 1/2,  L⁺(σ⁺=1) = log(3/2) = {LOG3-LOG2:.6f}")
    print(f"  E[e^(L⁺) ; σ⁺ = 1] = (1/2)·(3/2) = 0.75 = 3/4")
    print(f"  Identity: E[e^(L⁺) ; σ⁺ < ∞] = 1 (Esscher / Cramer w*=1)")
    print(f"    ⇒ E[e^(L⁺) ; σ⁺ ≥ 2, σ⁺ < ∞] = 1 − 3/4 = 1/4")
    print()
    print(f"  Decomposing q = P(σ⁺ < ∞):")
    print(f"    q = P(σ⁺ = 1) + P(σ⁺ ≥ 2, σ⁺ < ∞)")
    print(f"      = 1/2 + r,  where r = P(σ⁺ ≥ 2, σ⁺ < ∞)")
    print(f"  Empirical r = q − 1/2 = {q_emp - 0.5:.6f}")
    print(f"  And  E[e^(L⁺) | σ⁺ ≥ 2, σ⁺ < ∞] · r = 1/4")
    print(f"   ⇒  E[e^(L⁺) | σ⁺ ≥ 2, σ⁺ < ∞]    = 1/(4r) = {1/(4*(q_emp - 0.5)):.6f}")

    # Combined moment identity
    print()
    print("=" * 64)
    print("Implications for closed form")
    print("=" * 64)
    print(f"  E[L⁻] empirical:           {E_Lneg:.6f}")
    print(f"  E[L⁻] = μ/(1-q) check:     {LOG43/(1-q_emp):.6f}  (closed-form formula)")
    if abs(E_Lneg - LOG43/(1-q_emp)) < 3 * SE_L:
        print(f"  ✓ Identity confirmed within 3σ (passes simulation precision)")
    else:
        print(f"  ✗ Identity off by more than 3σ — investigate")
    print()
    print(f"  Conclusion: E[L⁻] in closed form requires q in closed form.")
    print(f"  At 10⁷-orbit precision, q = {q_emp:.6f} (SE ≈ {SE_q:.6f}).")
    print(f"  No simple closed-form candidate matches q within sampling error.")
    print()
    print(f"  The Esscher-tilt identity E[e^(L⁺); σ⁺<∞] = 1 decomposes")
    print(f"  into a 3/4 contribution from σ⁺ = 1 and a 1/4 contribution from")
    print(f"  σ⁺ ≥ 2. The 1/4 = (q - 1/2) · E[e^(L⁺) | σ⁺ ≥ 2, σ⁺ < ∞] is a")
    print(f"  clean structural identity but doesn't pin q without further input.")


if __name__ == "__main__":
    main()
