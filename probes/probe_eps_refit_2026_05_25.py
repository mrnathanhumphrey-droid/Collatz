"""
probe_eps_refit_2026_05_25.py
Two-mode ratio refit on ε_n trajectory (Nate's plan, 2026-05-25).

User's model: r_n = 1/2 + A·μ^n + B·μ_2^n
  - tests whether λ_3 (next subdominant eigenvalue after λ_2 = 1/2) falls out
    cleanly as μ_3 = μ * λ_2 = μ/2 (since correction rate to r_n is λ_3/λ_2).

Variants we ALSO fit (for sensitivity):
  - 1-mode: r_n = 1/2 + A·μ^n (3 params reduced to 2)
  - free leading: r_n = ρ + A·μ^n (3 params)
  - free leading 2-mode: r_n = ρ + A·μ^n + B·μ_2^n (5 params)

We exclude:
  - r_1 (transient, transient pole at n=1 from H1 spectral expansion)
  - r_8, r_9 (ε_9 ≈ 7.5e-6 near zero crossing — ratios blow up: r_8=0.010, r_9=-95.84)
  - We fit on r_2..r_7 (the "settling" regime), then SEPARATELY on r_10..r_15
    (post-crossing slow-mode regime). Report both.
"""
from __future__ import annotations
import json, os, sys
from fractions import Fraction
import numpy as np
from scipy.optimize import least_squares

sys.stdout.reconfigure(encoding="utf-8")

# ----------------------------------------------------------------------------
# Load ε_k data
# ----------------------------------------------------------------------------
OUTDIR = r"C:\Collatz\experiments_output"

# Exact ε_1..ε_8 (Fractions) from R77.7 v2 vec-pool run
with open(os.path.join(OUTDIR, "result_77_7_eps_exact_through_k8_v2_vec_pool.json")) as f:
    eps_exact_raw = json.load(f)
eps_exact = {}
for k_str, entry in eps_exact_raw.items():
    k = int(k_str)
    eps_exact[k] = Fraction(int(entry["num"]), int(entry["den"]))

# Numerical ε_9..ε_16 (from probe_epsilon_*/* JSONs and result_epsilon_*.csv)
# Use the highest-precision numerical values available.
eps_num = {
    9:  -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
    12: +2.2747137206e-03,
    13: +2.9482473172e-03,
    14: +3.5876674275e-03,
    15: +4.1611139460e-03,
    16: +4.6849767261e-03,
}

# Combine — float for fitting
eps = {}
for k in range(1, 9):
    eps[k] = float(eps_exact[k])
for k in range(9, 17):
    eps[k] = eps_num[k]

# Ratios r_n = ε_{n+1}/ε_n
ratios = {n: eps[n+1] / eps[n] for n in range(1, 16)}

print("=" * 70)
print("ε_k values (k=1..8 exact, k=9..16 numerical):")
print("=" * 70)
for k in range(1, 17):
    src = "exact" if k <= 8 else "num"
    print(f"  ε_{k:2d} = {eps[k]:+.6e}   ({src})")

print()
print("=" * 70)
print("Ratios r_n = ε_{n+1}/ε_n:")
print("=" * 70)
for n in range(1, 16):
    flag = ""
    if n == 1: flag = "  (transient)"
    elif n in (8, 9): flag = "  (near zero-crossing, ε_9 ≈ 7.5e-6 → ratio pathological)"
    elif n == 7: flag = "  (envelope jump at k=7)"
    print(f"  r_{n:2d} = {ratios[n]:+.6f}{flag}")

# ----------------------------------------------------------------------------
# Model definitions
# ----------------------------------------------------------------------------

def model_1mode_fixed(params, n_arr):
    """r_n = 1/2 + A·μ^n"""
    A, mu = params
    return 0.5 + A * mu**n_arr

def model_2mode_fixed(params, n_arr):
    """r_n = 1/2 + A·μ^n + B·μ_2^n"""
    A, mu, B, mu2 = params
    return 0.5 + A * mu**n_arr + B * mu2**n_arr

def model_1mode_free(params, n_arr):
    """r_n = ρ + A·μ^n (free leading rate)"""
    rho, A, mu = params
    return rho + A * mu**n_arr

def model_2mode_free(params, n_arr):
    """r_n = ρ + A·μ^n + B·μ_2^n (free leading + two corrections)"""
    rho, A, mu, B, mu2 = params
    return rho + A * mu**n_arr + B * mu2**n_arr

def fit_model(model, params0, n_arr, r_arr, bounds=None):
    def residuals(p):
        return model(p, n_arr) - r_arr
    kwargs = {}
    if bounds is not None:
        kwargs["bounds"] = bounds
    result = least_squares(residuals, params0, **kwargs)
    sse = float(np.sum(result.fun**2))
    return {
        "params": result.x.tolist(),
        "sse": sse,
        "n_data": len(r_arr),
        "n_params": len(params0),
        "residuals": result.fun.tolist(),
        "predictions": model(result.x, n_arr).tolist(),
        "success": bool(result.success),
    }

def summarize_fit(label, fit, names, n_arr, r_arr):
    print(f"\n--- {label} (data n={list(n_arr)}, K={len(r_arr)} points) ---")
    for name, val in zip(names, fit["params"]):
        print(f"  {name:>6} = {val:+.6f}")
    print(f"  SSE = {fit['sse']:.6e}    AIC ~ K·ln(SSE/K) + 2·p = "
          f"{len(r_arr)*np.log(fit['sse']/len(r_arr)) + 2*fit['n_params']:.2f}")
    print(f"  residuals: {['%+.4f' % r for r in fit['residuals']]}")

# ----------------------------------------------------------------------------
# WINDOW A: settling regime r_2..r_7 (Nate's "5-6 ratios")
# ----------------------------------------------------------------------------
print()
print("=" * 70)
print("WINDOW A — settling regime r_n for n=2..7 (6 ratios)")
print("=" * 70)
print("Note: r_2 = -0.535 is the sign-flip ratio (ε_2 +, ε_3 -). Including it.")

n_A = np.array([2, 3, 4, 5, 6, 7], dtype=float)
r_A = np.array([ratios[int(n)] for n in n_A])

# 1-mode fixed leading
fit_1f = fit_model(
    model_1mode_fixed,
    [-0.1, 0.5],
    n_A, r_A,
    bounds=([-10, 0.01], [10, 0.99]),
)
summarize_fit("1-mode, ρ=1/2 fixed:   r_n = 1/2 + A·μ^n",
              fit_1f, ["A", "μ"], n_A, r_A)

# 2-mode fixed leading (Nate's primary model)
fit_2f = fit_model(
    model_2mode_fixed,
    [-0.5, 0.5, 0.1, 0.7],
    n_A, r_A,
    bounds=([-100, 0.01, -100, 0.01], [100, 0.99, 100, 0.99]),
)
summarize_fit("2-mode, ρ=1/2 fixed:   r_n = 1/2 + A·μ^n + B·μ_2^n  ★",
              fit_2f, ["A", "μ", "B", "μ_2"], n_A, r_A)

# 1-mode free leading
fit_1F = fit_model(
    model_1mode_free,
    [0.4, -1.0, 0.5],
    n_A, r_A,
    bounds=([-2, -100, 0.01], [2, 100, 1.5]),
)
summarize_fit("1-mode, ρ free:        r_n = ρ + A·μ^n",
              fit_1F, ["ρ", "A", "μ"], n_A, r_A)

# ----------------------------------------------------------------------------
# WINDOW B — extended early regime r_2..r_6 (drops r_7 envelope jump)
# ----------------------------------------------------------------------------
print()
print("=" * 70)
print("WINDOW B — pre-jump regime r_n for n=2..6 (5 ratios, excludes r_7=2.36)")
print("=" * 70)
n_B = np.array([2, 3, 4, 5, 6], dtype=float)
r_B = np.array([ratios[int(n)] for n in n_B])

fit_1f_B = fit_model(
    model_1mode_fixed,
    [-0.1, 0.5],
    n_B, r_B,
    bounds=([-10, 0.01], [10, 0.99]),
)
summarize_fit("1-mode, ρ=1/2 fixed", fit_1f_B, ["A", "μ"], n_B, r_B)

fit_2f_B = fit_model(
    model_2mode_fixed,
    [-0.5, 0.5, 0.1, 0.7],
    n_B, r_B,
    bounds=([-100, 0.01, -100, 0.01], [100, 0.99, 100, 0.99]),
)
summarize_fit("2-mode, ρ=1/2 fixed  ★", fit_2f_B, ["A", "μ", "B", "μ_2"], n_B, r_B)

fit_1F_B = fit_model(
    model_1mode_free,
    [0.4, -1.0, 0.5],
    n_B, r_B,
    bounds=([-2, -100, 0.01], [2, 100, 1.5]),
)
summarize_fit("1-mode, ρ free", fit_1F_B, ["ρ", "A", "μ"], n_B, r_B)

# ----------------------------------------------------------------------------
# WINDOW C — drop r_2 (sign flip dominates) — r_3..r_6
# ----------------------------------------------------------------------------
print()
print("=" * 70)
print("WINDOW C — sign-stable regime r_n for n=3..6 (4 ratios, drops sign flip)")
print("=" * 70)
n_C = np.array([3, 4, 5, 6], dtype=float)
r_C = np.array([ratios[int(n)] for n in n_C])

fit_1f_C = fit_model(
    model_1mode_fixed,
    [-0.1, 0.5],
    n_C, r_C,
    bounds=([-10, 0.01], [10, 0.99]),
)
summarize_fit("1-mode, ρ=1/2 fixed", fit_1f_C, ["A", "μ"], n_C, r_C)

fit_1F_C = fit_model(
    model_1mode_free,
    [0.4, -1.0, 0.5],
    n_C, r_C,
    bounds=([-2, -100, 0.01], [2, 100, 1.5]),
)
summarize_fit("1-mode, ρ free", fit_1F_C, ["ρ", "A", "μ"], n_C, r_C)

# ----------------------------------------------------------------------------
# WINDOW D — post-crossing slow-mode regime r_10..r_15
# ----------------------------------------------------------------------------
print()
print("=" * 70)
print("WINDOW D — post-crossing slow-mode r_n for n=10..15 (6 ratios)")
print("=" * 70)
print("ε flipped sign at k=9→10. Ratios r_10..r_15: 2.084, 1.514, 1.296, 1.217, 1.160, 1.126")
print("Decelerating downward. NOT compatible with 1/2 leading.")

n_D = np.array([10, 11, 12, 13, 14, 15], dtype=float)
r_D = np.array([ratios[int(n)] for n in n_D])

fit_1F_D = fit_model(
    model_1mode_free,
    [1.0, 1.0, 0.9],
    n_D, r_D,
    bounds=([-2, -100, 0.01], [2, 100, 1.5]),
)
summarize_fit("1-mode, ρ free", fit_1F_D, ["ρ", "A", "μ"], n_D, r_D)

# ----------------------------------------------------------------------------
# EXTRAPOLATION CHECK — use Window B 2-mode fit to predict r_7, r_8, ...
# ----------------------------------------------------------------------------
print()
print("=" * 70)
print("EXTRAPOLATION — does Window B 2-mode fit predict r_7..r_15?")
print("=" * 70)
A, mu, B, mu2 = fit_2f_B["params"]
print(f"  Window B 2-mode: A={A:+.4f}, μ={mu:.4f}, B={B:+.4f}, μ_2={mu2:.4f}")
print(f"  {'n':>3} {'r_n actual':>13} {'r_n predicted':>15} {'residual':>12}")
for n in range(2, 16):
    pred = 0.5 + A * mu**n + B * mu2**n
    actual = ratios[n]
    print(f"  {n:>3} {actual:>+13.6f} {pred:>+15.6f} {actual-pred:>+12.6f}")

# ----------------------------------------------------------------------------
# λ_3 BACKOUT — under the model assumption λ_2 = 1/2, the correction rate
# (λ_3/λ_2) is what we measure as μ. So λ_3_candidate = μ * 1/2 = μ/2.
# ----------------------------------------------------------------------------
print()
print("=" * 70)
print("λ_3 CANDIDATES (under assumption λ_2 = 1/2)")
print("=" * 70)
print("Correction rate to r_n is λ_3/λ_2, so λ_3_candidate = μ · 1/2")
print()
print(f"  Window A 2-mode μ  = {fit_2f['params'][1]:.4f}  → λ_3_cand = {fit_2f['params'][1]/2:.4f}")
print(f"  Window A 2-mode μ_2= {fit_2f['params'][3]:.4f}  → λ_4_cand = {fit_2f['params'][3]/2:.4f}")
print(f"  Window B 2-mode μ  = {fit_2f_B['params'][1]:.4f}  → λ_3_cand = {fit_2f_B['params'][1]/2:.4f}")
print(f"  Window B 2-mode μ_2= {fit_2f_B['params'][3]:.4f}  → λ_4_cand = {fit_2f_B['params'][3]/2:.4f}")

# ----------------------------------------------------------------------------
# DUMP results to JSON
# ----------------------------------------------------------------------------
out_path = os.path.join(OUTDIR, "probe_eps_refit_2026_05_25.json")
with open(out_path, "w") as f:
    json.dump({
        "eps": {str(k): eps[k] for k in eps},
        "ratios": {str(n): ratios[n] for n in ratios},
        "window_A_n2_to_7": {
            "n": n_A.tolist(),
            "ratios": r_A.tolist(),
            "fits": {
                "1mode_fixed": fit_1f,
                "2mode_fixed": fit_2f,
                "1mode_free": fit_1F,
            },
        },
        "window_B_n2_to_6": {
            "n": n_B.tolist(),
            "ratios": r_B.tolist(),
            "fits": {
                "1mode_fixed": fit_1f_B,
                "2mode_fixed": fit_2f_B,
                "1mode_free": fit_1F_B,
            },
        },
        "window_C_n3_to_6": {
            "n": n_C.tolist(),
            "ratios": r_C.tolist(),
            "fits": {
                "1mode_fixed": fit_1f_C,
                "1mode_free": fit_1F_C,
            },
        },
        "window_D_n10_to_15": {
            "n": n_D.tolist(),
            "ratios": r_D.tolist(),
            "fits": {
                "1mode_free": fit_1F_D,
            },
        },
    }, f, indent=2)
print(f"\nDumped: {out_path}")
