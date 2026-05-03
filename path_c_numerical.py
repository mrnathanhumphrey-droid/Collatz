"""
path_c_numerical.py — Path C numerical verification of κ⁻(θ) closed form.

Workflow:
  1. Higher-precision simulation of E[L⁻] and E[σ⁻] (10M orbits).
  2. Direct numerical evaluation of Spitzer's formula
        log(1 - κ⁻(θ)) = -Σ_{n=1}^∞ (1/n) · E[e^(iθ·S_n); S_n < 0]
     for our walk under iid Geom(1/2), comparing to simulation.
  3. Test closed-form conjectures (E[σ⁻] = 7/2, E[L⁻] = (7/2)log(4/3), etc.)
  4. Compute Cramer-Esscher tilt parameters at w* = 1.

Cap: < 200 lines.
"""
import numpy as np
import sys
from scipy.special import gammaln

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3  # log(4/3)
LOG32 = LOG3 - LOG2       # log(3/2)


def simulate_ladder(n_orbits=10_000_000, seed=42, max_steps=200_000):
    """Simulate (L⁻, σ⁻) for the iid Geom(1/2) walk."""
    rng = np.random.default_rng(seed)
    L_vals = np.zeros(n_orbits)
    sigma_vals = np.zeros(n_orbits, dtype=np.int32)
    fails = 0
    for i in range(n_orbits):
        S = 0.0
        for n in range(1, max_steps + 1):
            v = rng.geometric(0.5)
            S += LOG3 - v * LOG2
            if S < 0:
                L_vals[i] = -S
                sigma_vals[i] = n
                break
        else:
            fails += 1
    return L_vals, sigma_vals, fails


def log_negbin(k, n):
    """log P(T_n = k) for T_n ~ NegBin(n, 1/2), k >= n."""
    if k < n:
        return -np.inf
    return gammaln(k) - gammaln(n) - gammaln(k - n + 1) - k * LOG2


def E_neg_S_indicator(n, k_max_factor=10):
    """E[-S_n · 1_{S_n < 0}] for our walk, numerically.
       S_n = n·log3 − T_n·log2, S_n < 0 ↔ T_n > n·log_2(3) ↔ T_n ≥ ⌈n·log_2(3)⌉.
       -S_n · 1_{S_n<0} = (T_n·log2 − n·log3) · 1_{T_n > n·log_2 3}.
    """
    log_2_3 = LOG3 / LOG2
    K_n = int(np.ceil(n * log_2_3))
    if (K_n - n * log_2_3) < 1e-12:
        K_n += 1  # strict inequality
    k_max = max(int(K_n + n * k_max_factor + 200), 500)
    total = 0.0
    for k in range(K_n, k_max + 1):
        log_p = log_negbin(k, n)
        if log_p < -50:
            continue
        p = np.exp(log_p)
        contrib = (k * LOG2 - n * LOG3) * p
        total += contrib
    return total


def spitzer_E_Lneg(n_max=2000):
    """E[L⁻] = Σ_{n=1}^∞ (1/n) · E[-S_n · 1_{S_n<0}]  via Pollaczek-Spitzer."""
    total = 0.0
    for n in range(1, n_max + 1):
        contrib = E_neg_S_indicator(n) / n
        total += contrib
    return total


def main():
    print("=" * 64)
    print("Path C numerical verification")
    print("=" * 64)

    # Step moments
    E_X = LOG3 - 2 * LOG2
    E_X2 = LOG3**2 - 4 * LOG3 * LOG2 + 6 * LOG2**2
    print(f"E[X]    = {E_X:.6f}  (= -log(4/3) = {-LOG43:.6f})")
    print(f"E[X²]   = {E_X2:.6f}")

    # Cramer root
    print(f"\n--- Cramer root analysis ---")
    print(f"E[e^(w·X)] = 3^w / (2·2^w - 1)")
    print(f"E[e^(0)] = 1 (trivial)")
    print(f"E[e^X] (w=1) = 3 / (2·2 - 1) = 3/3 = 1  ← Cramer root w* = 1")
    print(f"  number-theoretic identity: 2² − 1 = 3  (Collatz-structural)")
    # Esscher tilt parameters
    mu_star = LOG3 - (4 / 3) * LOG2
    print(f"\n--- Esscher tilt at w* = 1 ---")
    print(f"μ* = E_P*[X] = log(3) − (4/3)·log(2) = log(3 · 2^(-4/3))")
    print(f"      = {mu_star:.6f}  (positive, as expected)")
    print(f"P*(v = k) = 3 · 4^(-k)  (tilted v geometric)")
    print(f"E_P*[v] = 4/3 = {4.0/3:.6f}  (vs E_P[v] = 2)")

    # Higher-precision simulation
    print(f"\n--- Higher-precision simulation (10M orbits) ---")
    L_vals, sigma_vals, fails = simulate_ladder(n_orbits=10_000_000, seed=42)
    if fails:
        print(f"  WARNING: {fails} orbits failed to descend within max_steps")
    E_L = L_vals.mean()
    SE_L = L_vals.std() / np.sqrt(len(L_vals))
    E_sigma = sigma_vals.mean()
    SE_sigma = sigma_vals.std() / np.sqrt(len(sigma_vals))
    print(f"  E[L⁻]  = {E_L:.6f} ± {1.96*SE_L:.5f}  (95% CI) nats")
    print(f"          = {E_L/LOG43:.6f} step units")
    print(f"  E[σ⁻]  = {E_sigma:.6f} ± {1.96*SE_sigma:.5f}  (95% CI) steps")
    print(f"  Wald check: E[σ⁻]·log(4/3) = {E_sigma * LOG43:.6f}  (should = E[L⁻])")
    print(f"             diff = {E_sigma * LOG43 - E_L:+.6f}")

    # Spitzer formula numerical
    print(f"\n--- Spitzer formula numerical evaluation ---")
    print(f"  E[L⁻] = Σ_n (1/n) · E[-S_n · 1_{{S_n<0}}]   (Pollaczek-Spitzer)")
    print(f"  Truncating at n_max = 2000 ...")
    E_L_spitzer = spitzer_E_Lneg(n_max=2000)
    print(f"  E[L⁻] (Spitzer)  = {E_L_spitzer:.6f} nats = {E_L_spitzer/LOG43:.6f} steps")
    print(f"  E[L⁻] (sim)       = {E_L:.6f} nats = {E_L/LOG43:.6f} steps")
    print(f"  diff (Spitzer - sim) = {E_L_spitzer - E_L:+.6f} nats")

    # Closed-form conjectures
    print(f"\n--- Closed-form conjectures ---")
    candidates = [
        ("(7/2)·log(4/3)", 3.5 * LOG43, "E[σ⁻] = 7/2 ?"),
        ("log(e)", 1.0, "E[L⁻] = 1 nat exactly?"),
        ("log(8/3)/log(4/3) · log(4/3) = log(8/3)", LOG2 + LOG2 - LOG3, "log(8/3) ?"),
        ("log 2 + log(3/4)/2", LOG2 + (LOG3 - 2*LOG2)/2, "?"),
        ("log(2·4/3) = log(8/3)", np.log(8/3), "= log(8/3) double-check"),
        ("log 6 − log 2 = log 3", LOG3, "log(3) ?"),
    ]
    print(f"  {'Candidate':<35} {'Value (nats)':>14} {'vs E[L⁻]':>12}")
    for name, val, note in candidates:
        diff = E_L - val
        z = diff / SE_L
        match = "✓" if abs(z) < 2 else " "
        print(f"  {name:<35} {val:>14.6f} {diff:>+12.6f} ({z:+.1f}σ) {match}  {note}")

    # Specific test of (7/2)·log(4/3)
    cand_72 = 3.5 * LOG43
    print(f"\n  Test E[L⁻] = (7/2)·log(4/3) = {cand_72:.6f} nats")
    print(f"      simulated E[L⁻] = {E_L:.6f} ± {1.96*SE_L:.5f}")
    print(f"      diff = {E_L - cand_72:+.6f} nats ({(E_L-cand_72)/SE_L:+.2f}σ)")

    # Lorden value cross-check (asymptotic residual life)
    print(f"\n--- Lorden asymptotic residual life ---")
    W_lorden = E_X2 / (2 * abs(E_X))
    print(f"  E[X²]/(2|E[X]|) = {W_lorden:.6f} nats = {W_lorden/LOG43:.6f} steps")
    print(f"  This is the asymptotic limit of overshoot at first crossing of -y as y → ∞.")
    print(f"  Differs from strict descending ladder mean E[L⁻] = {E_L:.4f} nats.")

    # Higher moments
    print(f"\n--- Higher moments of L⁻ ---")
    E_L2 = (L_vals**2).mean()
    Var_L = E_L2 - E_L**2
    print(f"  E[L⁻²]  = {E_L2:.6f}")
    print(f"  Var[L⁻] = {Var_L:.6f}")
    print(f"  E[L⁻²]/(2·E[L⁻]) = {E_L2/(2*E_L):.6f} (= asymptotic Lorden formula?)")
    print(f"        compared to E[X²]/(2|E[X]|) = {W_lorden:.6f}")


if __name__ == "__main__":
    main()
