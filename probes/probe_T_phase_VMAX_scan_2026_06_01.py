"""
V_MAX sensitivity scan of phase-twisted T_3 + PSLQ identification of eigenvalues.

V_MAX values: 20, 30, 40, 53 (float64 underflow boundary).
At each V_MAX, build T and extract top 18 eigenvalues.
Compare across V_MAX to identify which are V_MAX-sensitive vs converged.
For converged eigenvalues, PSLQ against:
  - Simple rationals {1/k for k=1..30}
  - Sqrt rationals
  - log(2), log(3), pi/k
  - 3/(4 chi(2) - 1) family at q=3
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np
from mpmath import mp, mpf, pslq, log as mplog, pi as mppi, sqrt as mpsqrt
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30

N = 3
N_dom = 3 ** N
N_cod = 3 ** (N + 1)
coprime_n = [r for r in range(N_dom) if r % 3 != 0]
idx_n = {r: i for i, r in enumerate(coprime_n)}
N_eta = len(coprime_n)
DIM = N_eta * N_dom

def build_T(V_MAX):
    Z = 1.0 - 2.0 ** (-V_MAX)
    T = np.zeros((DIM, DIM), dtype=np.complex128)
    inv2_full = [pow(2, -v, N_cod) for v in range(V_MAX + 1)]
    inv2_mod3 = [pow(2, -v, 3) if v > 0 else 1 for v in range(V_MAX + 1)]
    pow2_dom = [pow(2, v, N_dom) for v in range(V_MAX + 1)]
    pow2_inv_dom = [pow(2, -v, N_dom) for v in range(V_MAX + 1)]
    coprime_np = [r for r in range(N_cod) if r % 3 != 0]
    for eta_prime in coprime_np:
        eta_target = eta_prime % N_dom
        if eta_target == 0 or eta_target % 3 == 0:
            continue
        eta_target_idx = idx_n[eta_target]
        eta_prime_mod3 = eta_prime % 3
        for v in range(1, V_MAX + 1):
            inv2_v_full = inv2_full[v]
            inv2_v_mod3 = inv2_mod3[v]
            pow2_v_dom = pow2_dom[v]
            for vp in range(1, V_MAX + 1):
                inv2_vp_full = inv2_full[vp]
                inv2_vp_mod3 = inv2_mod3[vp]
                required_mod3 = (eta_prime_mod3 * inv2_vp_mod3 - inv2_v_mod3) % 3
                if v >= vp:
                    pow2_vmvp_dom = pow2_dom[v - vp]
                else:
                    pow2_vmvp_dom = pow2_inv_dom[vp - v]
                eta_input = (eta_target * pow2_vmvp_dom) % N_dom
                if eta_input == 0 or eta_input % 3 == 0:
                    continue
                eta_input_idx = idx_n[eta_input]
                weight = 3.0 * (2.0 ** (-v - vp)) / (Z * Z)
                if weight < 1e-300:
                    continue  # underflow guard
                for delta_prime in range(required_mod3, N_cod, 3):
                    num = (delta_prime + inv2_v_full - eta_prime * inv2_vp_full) % N_cod
                    if num % 3 != 0:
                        continue
                    num_div3 = num // 3
                    delta_input = (num_div3 * pow2_v_dom) % N_dom
                    delta_target = delta_prime % N_dom
                    row = eta_target_idx * N_dom + delta_target
                    col = eta_input_idx * N_dom + delta_input
                    T[row, col] += weight
    T_avg = T / 9.0
    eigs = np.linalg.eigvals(T_avg)
    return eigs

print(f"V_MAX sensitivity scan at n=3, DIM={DIM}")
print()

V_MAX_list = [20, 30, 40, 53]
results = {}
for V_MAX in V_MAX_list:
    t0 = time.time()
    eigs = build_T(V_MAX)
    # Get real parts of top 18 (after dominant 1)
    abs_eigs = np.abs(eigs)
    order = np.argsort(-abs_eigs)
    top = eigs[order][:25]
    # Filter near-zero
    top_real = np.array([e.real for e in top if abs(e) > 1e-8])
    results[V_MAX] = top_real
    print(f"V_MAX={V_MAX}  build+eig time={time.time()-t0:.1f}s")
    for i, e in enumerate(top_real[:18]):
        print(f"  {i:>3}: {e:.15f}")
    print()

# Now compare across V_MAX: for each eigenvalue index, see how it changes
print(f"\nConvergence check — top 18 eigenvalues across V_MAX:")
print(f"  {'idx':>3} {'V_MAX=20':>17} {'V_MAX=30':>17} {'V_MAX=40':>17} {'V_MAX=53':>17}  {'Δ(53-40)':>12}")
for i in range(18):
    vals = [results[v][i] if i < len(results[v]) else float('nan') for v in V_MAX_list]
    delta = vals[3] - vals[2] if not (np.isnan(vals[2]) or np.isnan(vals[3])) else float('nan')
    print(f"  {i:>3} {vals[0]:>17.13f} {vals[1]:>17.13f} {vals[2]:>17.13f} {vals[3]:>17.13f}  {delta:>+12.3e}")

# PSLQ each V_MAX=53 eigenvalue against a basis
print(f"\n=== PSLQ identification of V_MAX=53 eigenvalues ===")
eigs_53 = results[53][:18]
basis = [mpf(1), mpf(1)/2, mpf(1)/3, mpf(1)/4, mpf(1)/5, mpf(1)/6, mpf(1)/7, mpf(1)/8, mpf(1)/9, mpf(1)/10, mpf(1)/11, mpf(1)/13, mpf(1)/15, mpf(1)/17, mpf(1)/21, mpf(1)/27]
bnames = ["1", "1/2", "1/3", "1/4", "1/5", "1/6", "1/7", "1/8", "1/9", "1/10", "1/11", "1/13", "1/15", "1/17", "1/21", "1/27"]

for i, eig in enumerate(eigs_53):
    eig_mp = mpf(eig)
    # Direct check: does it match a single basis element?
    best_match = None
    best_diff = mpf('inf')
    for j, b in enumerate(basis):
        d = abs(eig_mp - b)
        if d < best_diff:
            best_diff = d
            best_match = (j, b, bnames[j])
    j, b, name = best_match
    # PSLQ: eig vs full basis (look for small integer relation)
    try:
        rel = pslq([eig_mp] + basis, tol=mpf(10)**(-15), maxcoeff=100)
        if rel is not None and rel[0] != 0:
            terms = [f"({rel[0]:+d})*λ"]
            for c, n_ in zip(rel[1:], bnames):
                if c != 0:
                    terms.append(f"({c:+d})*{n_}")
            pslq_str = " ".join(terms)
        else:
            pslq_str = "no relation found"
    except Exception as e:
        pslq_str = f"PSLQ error: {e}"
    diff_to_best = abs(eig - float(b))
    marker = " <-- EXACT MATCH" if diff_to_best < 1e-10 else (" (close)" if diff_to_best < 1e-3 else "")
    print(f"  λ[{i}]={eig:.13f}  closest={name}={float(b):.10f}  diff={diff_to_best:.3e}{marker}")
    print(f"    PSLQ: {pslq_str[:120]}")

# Specifically look for closed forms of the non-obvious eigenvalues
# 0.806, 0.5166, 0.2323, 0.1756, 0.1240, 0.1142
print(f"\n=== Try algebraic forms for non-trivial eigenvalues ===")
suspicious_eigs = [results[53][i] for i in range(min(18, len(results[53])))]
extended_basis = [
    mpf(1), mpf(1)/2, mpf(1)/3, mpf(1)/5, mpf(1)/7, mpf(1)/9,
    mpsqrt(mpf(2))/2, mpsqrt(mpf(3))/2, mpsqrt(mpf(5))/2,
    mplog(mpf(2)), mplog(mpf(3)),
    mpf(1)/mpsqrt(mpf(2)), mpf(1)/mpsqrt(mpf(3)),
    mpf(2)/3, mpf(3)/5, mpf(4)/5, mpf(5)/6, mpf(7)/9, mpf(8)/9,
]
ext_names = ["1", "1/2", "1/3", "1/5", "1/7", "1/9",
             "√2/2", "√3/2", "√5/2",
             "log2", "log3",
             "1/√2", "1/√3",
             "2/3", "3/5", "4/5", "5/6", "7/9", "8/9"]
print(f"Extended basis identification (eigs that don't match simple 1/k):")
for i, eig in enumerate(suspicious_eigs):
    eig_mp = mpf(eig)
    best = (None, mpf('inf'), "")
    for b, n_ in zip(extended_basis, ext_names):
        d = abs(eig_mp - b)
        if d < best[1]:
            best = (n_, d, float(b))
    if best[1] < 1e-3:
        marker = " <-- match" if best[1] < 1e-6 else " (close)"
        print(f"  λ[{i}]={eig:.13f}  closest in extended: {best[0]} ≈ {best[2]:.10f}  diff={float(best[1]):.3e}{marker}")
    else:
        print(f"  λ[{i}]={eig:.13f}  no extended-basis match within 1e-3")

# Save
output = {
    "N": N, "DIM": DIM,
    "V_MAX_scan": {str(v): results[v].tolist() for v in V_MAX_list},
}
with open("C:/Collatz/T_phase_VMAX_scan_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: C:/Collatz/T_phase_VMAX_scan_2026_06_01.json")
