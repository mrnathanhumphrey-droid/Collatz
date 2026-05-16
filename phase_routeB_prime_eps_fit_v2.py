"""
phase_routeB_prime_eps_fit_v2.py — Refined fit: transient + slow-mode CC pair.

Model: ε_k = A·(1/r_trans)^k + B·(1/r_slow)^k · cos(θ k + φ)
  - First term: transient z=2 branch-cut decay (r_trans ≈ 2, |ε|~(1/2)^k).
  - Second term: slow-mode CC pair (r_slow ≈ 1.016 = 1/0.984, θ ≈ 2π/9.2).

Fit (A, r_trans, B, r_slow, θ, φ) to ε_k for k=2..13. Use log-loss if scale varies.
"""
import sys, numpy as np
import json, os
from scipy.optimize import differential_evolution, minimize
sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)

EPS = {
    1:  +2.0000000000e-01,
    2:  +9.5238095238e-03,
    3:  -5.0919863259e-03,
    4:  -2.4522582483e-03,
    5:  -1.1517469151e-03,
    6:  -4.9790566522e-04,
    7:  -1.1752368304e-03,
    8:  -7.4554636729e-04,
    9:  -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
    12: +2.2747137206e-03,
    13: +2.9482473172e-03,
}


def model(params, k_arr):
    A, r_trans, B, r_slow, theta, phi = params
    return A * (1.0 / r_trans) ** k_arr + B * (1.0 / r_slow) ** k_arr * np.cos(theta * k_arr + phi)


def loss(params, k_arr, eps_arr):
    pred = model(params, k_arr)
    return np.sum((pred - eps_arr) ** 2)


def main():
    out = {}
    for (k_lo, k_hi) in [(2, 13), (3, 13), (4, 13)]:
        k_arr = np.arange(k_lo, k_hi + 1, dtype=float)
        eps_arr = np.array([EPS[int(k)] for k in k_arr])
        # bounds:
        #   A free; r_trans near 2; B free; r_slow near 1.0..2.0; theta ∈ (0, π); phi ∈ (-π, π).
        bounds = [
            (-1.0, 1.0),          # A
            (1.5, 4.0),           # r_trans (transient radius)
            (-1.0, 1.0),          # B
            (1.0, 2.5),           # r_slow (slow-mode radius)
            (0.05, np.pi),        # theta (oscillation freq)
            (-np.pi, np.pi),      # phi
        ]
        best = None
        for seed in range(20):
            de = differential_evolution(loss, bounds, args=(k_arr, eps_arr),
                                         seed=seed, tol=1e-15, maxiter=5000, polish=True,
                                         popsize=40, mutation=(0.3, 1.7), recombination=0.9)
            if best is None or de.fun < best.fun:
                best = de
        params = best.x
        A, r_trans, B, r_slow, theta, phi = params

        print(f"\n=== k={k_lo}..{k_hi}, best of 20 seeds ===")
        print(f"  A = {A:+.4e},  r_trans = {r_trans:.6f}  (1/r = {1/r_trans:.6f})")
        print(f"  B = {B:+.4e},  r_slow  = {r_slow:.6f}  (1/r_slow = {1/r_slow:.6f})    [empirical predicted: 1/0.984 = 1.0163]")
        print(f"  theta = {theta:.6f} rad,  period = {2*np.pi/theta:.4f}  [empirical: 9.2]")
        print(f"  phi   = {phi:+.6f}")
        print(f"  loss  = {best.fun:.6e}")
        pred = model(params, k_arr)
        print(f"  Fit table:")
        for k, p, a in zip(k_arr, pred, eps_arr):
            err = p - a
            print(f"    k={int(k):2d}: predicted={p:+.4e}, actual={a:+.4e}, err={err:+.4e}")
        out[f"k={k_lo}..{k_hi}"] = {
            "A": float(A), "r_trans": float(r_trans),
            "B": float(B), "r_slow": float(r_slow),
            "theta": float(theta), "phi": float(phi),
            "period": float(2 * np.pi / theta),
            "decay_rate_slow": float(1.0 / r_slow),
            "loss": float(best.fun),
        }

    with open(os.path.join(OUTDIR, "phase_routeB_prime_eps_fit_v2.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote phase_routeB_prime_eps_fit_v2.json")


if __name__ == "__main__":
    main()
