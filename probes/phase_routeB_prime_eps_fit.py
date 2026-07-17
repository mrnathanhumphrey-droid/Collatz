"""
phase_routeB_prime_eps_fit.py — Route B' direct fit of CC-pair model to ε_k.

Hypothesis: in the post-transient regime k=7..13, ε_k ≈ A·ρ^k cos(θ k + φ) where (ρ, θ)
characterize the dominant slow-mode CC pair.

Fit (A, ρ, θ, φ) via nonlinear least squares. Compare ρ to empirical PADE prediction
0.984 and θ to period 9.2 (i.e., θ ≈ 2π/9.2 = 0.683 rad).

Sub-fits: (a) k=7..13 only (post-transient); (b) k=4..13 (full mid-late regime).
"""
import sys, numpy as np
import json, os
from scipy.optimize import minimize, differential_evolution
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


def model_cc(params, k_arr):
    A, rho, theta, phi = params
    return A * (rho ** k_arr) * np.cos(theta * k_arr + phi)


def loss_cc(params, k_arr, eps_arr):
    pred = model_cc(params, k_arr)
    return np.sum((pred - eps_arr) ** 2)


def model_2mode(params, k_arr):
    """Two modes: dominant CC + secondary geometric. ε_k = A·ρ^k cos(θk + φ) + B·s^k."""
    A, rho, theta, phi, B, s = params
    return A * (rho ** k_arr) * np.cos(theta * k_arr + phi) + B * (s ** k_arr)


def loss_2mode(params, k_arr, eps_arr):
    return np.sum((model_2mode(params, k_arr) - eps_arr) ** 2)


def fit_window(k_lo, k_hi, model_kind="cc"):
    k_arr = np.arange(k_lo, k_hi + 1, dtype=float)
    eps_arr = np.array([EPS[int(k)] for k in k_arr])

    if model_kind == "cc":
        # Global search via differential_evolution then local refinement
        bounds = [(-0.1, 0.1), (0.5, 1.2), (0.1, 1.5), (-np.pi, np.pi)]
        de = differential_evolution(loss_cc, bounds, args=(k_arr, eps_arr),
                                     seed=42, tol=1e-12, maxiter=2000, polish=True)
        params = de.x
        loss_val = de.fun
        A, rho, theta, phi = params
        return {
            "model": "cc",
            "k_range": (k_lo, k_hi),
            "A": float(A), "rho": float(rho), "theta": float(theta), "phi": float(phi),
            "period": float(2 * np.pi / theta) if theta > 0 else None,
            "loss": float(loss_val),
            "predicted": [float(v) for v in model_cc(params, k_arr)],
            "actual": [float(v) for v in eps_arr],
        }
    elif model_kind == "2mode":
        bounds = [(-0.1, 0.1), (0.5, 1.2), (0.1, 1.5), (-np.pi, np.pi), (-0.1, 0.1), (-1.0, 1.0)]
        de = differential_evolution(loss_2mode, bounds, args=(k_arr, eps_arr),
                                     seed=42, tol=1e-12, maxiter=3000, polish=True)
        params = de.x
        loss_val = de.fun
        A, rho, theta, phi, B, s = params
        return {
            "model": "2mode",
            "k_range": (k_lo, k_hi),
            "A": float(A), "rho": float(rho), "theta": float(theta), "phi": float(phi),
            "B": float(B), "s": float(s),
            "period": float(2 * np.pi / theta) if theta > 0 else None,
            "loss": float(loss_val),
            "predicted": [float(v) for v in model_2mode(params, k_arr)],
            "actual": [float(v) for v in eps_arr],
        }


def main():
    out = {}
    for (k_lo, k_hi) in [(7, 13), (5, 13), (4, 13), (3, 13)]:
        for model_kind in ["cc", "2mode"]:
            label = f"k={k_lo}..{k_hi}_{model_kind}"
            print(f"\n=== {label} ===")
            try:
                res = fit_window(k_lo, k_hi, model_kind)
            except Exception as e:
                print(f"  ERROR: {e}")
                continue
            print(f"  model: {res['model']}")
            if res["model"] == "cc":
                print(f"  ρ = {res['rho']:.6f}    (empirical 0.984)")
                print(f"  θ = {res['theta']:.6f} rad   period = {res['period']:.4f}    (empirical 9.2)")
                print(f"  A = {res['A']:.4e}   φ = {res['phi']:+.4f}")
            else:
                print(f"  ρ = {res['rho']:.6f}    (empirical 0.984)")
                print(f"  θ = {res['theta']:.6f} rad   period = {res['period']:.4f}    (empirical 9.2)")
                print(f"  A = {res['A']:.4e}   φ = {res['phi']:+.4f}")
                print(f"  B = {res['B']:.4e}   s = {res['s']:+.6f}")
            print(f"  loss = {res['loss']:.6e}")
            print(f"  fit vs actual:")
            for k, p, a in zip(np.arange(k_lo, k_hi + 1), res["predicted"], res["actual"]):
                print(f"    k={int(k):2d}: predicted={p:+.4e}, actual={a:+.4e}, err={p-a:+.4e}")
            out[label] = res

    with open(os.path.join(OUTDIR, "phase_routeB_prime_eps_fit.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote phase_routeB_prime_eps_fit.json")


if __name__ == "__main__":
    main()
