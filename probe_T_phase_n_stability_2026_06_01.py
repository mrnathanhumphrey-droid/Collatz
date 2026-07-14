"""
n-stability test of V_n^M_phase spectrum.

Build T_n at n=2 (dim 54), n=3 (dim 486), n=4 (dim 4374) all at V_MAX=40.
Compare top eigenvalues across n.

Outcome:
  - If top eigenvalues stable across n=2,3,4 -> V_n^M_phase IS a carrier of
    the asymptotic rate; Paper 3 Section 5's negative claim needs updating
    for this operator.
  - If they drift -> finite-truncation artifact; Section 5 stands.

Critical comparison values from n=3 V_MAX=53:
  0.8056527712398  <-- candidate asymptotic rate
  0.5165763600509
  0.2322660628363
  0.1756051838662
  0.1240052027887
  0.1141712757429
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

V_MAX = 40

def build_T_n(n_target):
    N_dom = 3 ** n_target
    N_cod = 3 ** (n_target + 1)
    coprime_n = [r for r in range(N_dom) if r % 3 != 0]
    coprime_np = [r for r in range(N_cod) if r % 3 != 0]
    idx_n = {r: i for i, r in enumerate(coprime_n)}
    DIM = len(coprime_n) * N_dom

    Z = 1.0 - 2.0 ** (-V_MAX)
    T = np.zeros((DIM, DIM), dtype=np.complex128)
    inv2_full = [pow(2, -v, N_cod) for v in range(V_MAX + 1)]
    inv2_mod3 = [pow(2, -v, 3) if v > 0 else 1 for v in range(V_MAX + 1)]
    pow2_dom = [pow(2, v, N_dom) for v in range(V_MAX + 1)]
    pow2_inv_dom = [pow(2, -v, N_dom) for v in range(V_MAX + 1)]

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
                    continue
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
    return T_avg

print(f"V_n^M_phase n-stability test, V_MAX={V_MAX}")
print()

results = {}
for n in [2, 3, 4]:
    DIM = (2 * 3 ** (n - 1)) * 3 ** n
    print(f"--- n={n}  DIM={DIM} ---", flush=True)
    t0 = time.time()
    T = build_T_n(n)
    print(f"  Built in {time.time()-t0:.1f}s", flush=True)
    t1 = time.time()
    eigs = np.linalg.eigvals(T)
    print(f"  Diagonalized in {time.time()-t1:.1f}s", flush=True)
    order = np.argsort(-np.abs(eigs))
    eigs_sorted = eigs[order]
    # Top 25 by |λ|
    top = [(complex(e.real, e.imag), float(abs(e)), float(np.angle(e)*180/np.pi)) for e in eigs_sorted[:25]]
    results[n] = top
    print(f"  Top 25 |λ|:")
    for i, (lam, abs_lam, arg_deg) in enumerate(top):
        complex_marker = " [COMPLEX]" if abs(lam.imag) > 1e-9 else ""
        print(f"    {i:>3}: λ={lam.real:+.13f}{lam.imag:+.13f}j  |λ|={abs_lam:.13f}  arg={arg_deg:+8.3f}°{complex_marker}")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print()

# Cross-n comparison: top 15 magnitudes at each n
print(f"=" * 80)
print(f"CROSS-N EIGENVALUE STABILITY (|λ| only)")
print(f"=" * 80)
print(f"  {'idx':>3} {'n=2':>17} {'n=3':>17} {'n=4':>17}  {'Δ(4-3)':>14}  {'Δ(3-2)':>14}")
n_show = min(len(results[2]), len(results[3]), len(results[4]), 18)
for i in range(n_show):
    vals = [results[n][i][1] for n in [2, 3, 4]]
    delta_43 = vals[2] - vals[1]
    delta_32 = vals[1] - vals[0]
    marker = ""
    if abs(delta_43) < 1e-6:
        marker = " <-- stable"
    elif abs(delta_43) < 1e-3:
        marker = " <-- mostly stable"
    elif abs(delta_43) > 0.05:
        marker = " <-- DRIFT"
    print(f"  {i:>3} {vals[0]:>17.13f} {vals[1]:>17.13f} {vals[2]:>17.13f}  {delta_43:>+14.3e}  {delta_32:>+14.3e}{marker}")

# Specifically check the headline values from n=3 V_MAX=53
print(f"\n=== Cross-n check of n=3 V_MAX=53 reference values ===")
n3_targets = [0.8056527712398, 0.5165763600509, 0.2322660628363,
              0.1756051838662, 0.1240052027887, 0.1141712757429]
for target in n3_targets:
    print(f"\n  Target = {target:.13f}:")
    for n in [2, 3, 4]:
        best_diff = float('inf')
        best_lam = None
        for lam, abs_lam, arg_deg in results[n]:
            d = abs(abs_lam - target)
            if d < best_diff:
                best_diff = d
                best_lam = (lam, abs_lam, arg_deg)
        if best_lam:
            lam, abs_lam, arg = best_lam
            marker = " <-- present" if best_diff < 1e-3 else (" (drift)" if best_diff < 0.05 else " <-- ABSENT")
            print(f"    n={n}: closest |λ|={abs_lam:.13f}  diff={best_diff:+.3e}  arg={arg:+.2f}°{marker}")

# Save
output = {
    "V_MAX": V_MAX,
    "results": {
        str(n): [{"re": l[0].real, "im": l[0].imag, "abs": l[1], "arg_deg": l[2]} for l in results[n]]
        for n in [2, 3, 4]
    },
    "n3_targets": n3_targets,
}
with open("C:/Collatz/T_phase_n_stability_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: C:/Collatz/T_phase_n_stability_2026_06_01.json")
