"""
Build the FULL phase-twisted T operator on V_n^M_phase at q=3, diagonalize.

V_n^M_phase: complex-valued functions on (Z/3^n)* x (Z/3^n).
Dim = 2 * 3^(n-1) * 3^n = 2 * 3^(2n-1).
  n=2: 54
  n=3: 486
  n=4: 4374

T_n is a self-map via projection (eta_prime, delta_prime) -> (eta_prime mod 3^n, delta_prime mod 3^n).

Recursion (derived from extending T_M_truncated_spectrum.py to allow output phase != 0):

  M_{n+1}_phase(eta', delta') = (3/Z^2) sum_{v,v'} 2^(-v-v')
      * 1[delta' ≡ eta'*2^(-v') - 2^(-v) (mod 3)]
      * M_n_phase(eta_n * 2^(v-v') mod 3^n, delta_n)

where:
  eta_n = eta' mod 3^n
  delta_n = ((delta' + 2^(-v) - eta'*2^(-v')) / 3) * 2^v mod 3^n

To get T_n as self-map, project (eta', delta') -> (eta_target, delta_target) at level n
by mod 3^n. Each (eta_target, delta_target) gets contributions from 9 lifts at level n+1.

DIAGONALIZE and look for: complex pair near rho = 0.19, theta = +-41°
(matching Model 3 fit of eps_k data).
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

N = 3  # target level
V_MAX = 20

N_dom = 3 ** N         # 3^n  (level-n domain size)
N_cod = 3 ** (N + 1)   # 3^(n+1)  (level-n+1 codomain size)

coprime_n = [r for r in range(N_dom) if r % 3 != 0]
coprime_np = [r for r in range(N_cod) if r % 3 != 0]
idx_n = {r: i for i, r in enumerate(coprime_n)}
N_eta = len(coprime_n)
DIM = N_eta * N_dom

print(f"Building T_phase at q=3, n={N}, V_MAX={V_MAX}")
print(f"  |coprime_n|={N_eta}, N_dom={N_dom}, N_cod={N_cod}")
print(f"  V_n^M_phase dim = {DIM}")
print()

Z = 1.0 - 2.0 ** (-V_MAX)
T = np.zeros((DIM, DIM), dtype=np.complex128)

# Precompute inverses mod N_cod
inv2_full = [pow(2, -v, N_cod) for v in range(V_MAX + 1)]  # inv2_full[v] = 2^(-v) mod N_cod, inv2_full[0]=1
inv2_mod3 = [pow(2, -v, 3) if v > 0 else 1 for v in range(V_MAX + 1)]
pow2_dom = [pow(2, v, N_dom) for v in range(V_MAX + 1)]
pow2_inv_dom = [pow(2, -v, N_dom) for v in range(V_MAX + 1)]

t0 = time.time()
n_terms = 0

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

            # Required mod-3 residue of delta_prime
            required_mod3 = (eta_prime_mod3 * inv2_vp_mod3 - inv2_v_mod3) % 3

            # Compute eta_input = (eta_target) * 2^(v-vp) mod N_dom
            if v >= vp:
                pow2_vmvp_dom = pow2_dom[v - vp]
            else:
                pow2_vmvp_dom = pow2_inv_dom[vp - v]
            eta_input = (eta_target * pow2_vmvp_dom) % N_dom
            if eta_input == 0 or eta_input % 3 == 0:
                continue
            eta_input_idx = idx_n[eta_input]

            weight = 3.0 * (2.0 ** (-v - vp)) / (Z * Z)

            # Iterate over delta_prime values satisfying mod-3 constraint
            for delta_prime in range(required_mod3, N_cod, 3):
                # Compute the "numerator" — must be divisible by 3 by mod-3 constraint
                num = (delta_prime + inv2_v_full - eta_prime * inv2_vp_full) % N_cod
                if num % 3 != 0:
                    continue  # paranoid; shouldn't happen
                num_div3 = num // 3
                delta_input = (num_div3 * pow2_v_dom) % N_dom
                delta_target = delta_prime % N_dom

                row = eta_target_idx * N_dom + delta_target
                col = eta_input_idx * N_dom + delta_input
                T[row, col] += weight
                n_terms += 1

print(f"Built T in {time.time()-t0:.1f}s, n_terms accumulated: {n_terms}", flush=True)

# Project: each (eta_target, delta_target) receives contributions from 9 lifts.
# Divide T row-wise by 9 (the lift count) to average — matching T_M_truncated's averaging convention.
T_avg = T / 9.0

# Diagonalize T_avg
print(f"Diagonalizing T_avg ({DIM}x{DIM} complex)...", flush=True)
t1 = time.time()
eigs = np.linalg.eigvals(T_avg)
print(f"Eigenvalue compute: {time.time()-t1:.1f}s", flush=True)

# Sort by |λ|
order = np.argsort(-np.abs(eigs))
eigs_sorted = eigs[order]

# Report top eigenvalues
print(f"\nTop 30 eigenvalues of T_avg (sorted by |λ|):")
print(f"  {'idx':>3} {'Re(λ)':>14} {'Im(λ)':>14} {'|λ|':>14} {'arg deg':>12}  match")
for i, lam in enumerate(eigs_sorted[:30]):
    abs_lam = abs(lam)
    arg_deg = np.angle(lam) * 180 / np.pi
    marker = ""
    if 0.18 <= abs_lam <= 0.20 and 38 <= abs(arg_deg) <= 44:
        marker = "  <-- COTAESCU MATCH (Model 3 rho=0.19, theta=41°)"
    elif 0.23 <= abs_lam <= 0.25 and 14 <= abs(arg_deg) <= 18:
        marker = "  <-- COTAESCU MATCH (Model 1 rho=0.24, theta=16°)"
    elif 0.45 <= abs_lam <= 0.55 and abs(arg_deg) < 5:
        marker = "  <-- R77.2 rate-1/2 (real)"
    print(f"  {i:>3} {lam.real:>+14.6e} {lam.imag:>+14.6e} {abs_lam:>14.6e} {arg_deg:>+12.4f}°{marker}")
print()

# Compare to truncated T_M (should match dominant: 1/3)
print(f"Truncated T_M_trunc had dominant eigenvalue: 1/3 = {1/3:.6f}")
print(f"Our T_avg dominant eigenvalue: {abs(eigs_sorted[0]):.6f}")
print()

# Specific check: search for eigenvalues matching Cotaescu fit zones
print(f"Cotaescu Model 3 target: |λ| ≈ 0.190, arg ≈ ±41°")
print(f"Cotaescu Model 1 target: |λ| ≈ 0.236, arg ≈ ±16°")
hits_m3 = []
hits_m1 = []
hits_05_real = []
for i, lam in enumerate(eigs_sorted):
    abs_lam = abs(lam)
    arg_deg = np.angle(lam) * 180 / np.pi
    if 0.15 <= abs_lam <= 0.25 and 35 <= abs(arg_deg) <= 47:
        hits_m3.append((i, lam))
    if 0.20 <= abs_lam <= 0.28 and 10 <= abs(arg_deg) <= 22:
        hits_m1.append((i, lam))
    if 0.45 <= abs_lam <= 0.55 and abs(arg_deg) < 5:
        hits_05_real.append((i, lam))

print(f"\nModel 3 zone hits ({len(hits_m3)}):")
for i, lam in hits_m3[:10]:
    print(f"  idx={i}: |λ|={abs(lam):.6f}, arg={np.angle(lam)*180/np.pi:+.3f}°")
print(f"\nModel 1 zone hits ({len(hits_m1)}):")
for i, lam in hits_m1[:10]:
    print(f"  idx={i}: |λ|={abs(lam):.6f}, arg={np.angle(lam)*180/np.pi:+.3f}°")
print(f"\nReal 1/2 zone hits ({len(hits_05_real)}):")
for i, lam in hits_05_real[:10]:
    print(f"  idx={i}: |λ|={abs(lam):.6f}, arg={np.angle(lam)*180/np.pi:+.3f}°")

# Save
output = {
    "N": N, "V_MAX": V_MAX, "DIM": DIM,
    "dominant_eig": {"re": float(eigs_sorted[0].real), "im": float(eigs_sorted[0].imag), "abs": float(abs(eigs_sorted[0])), "arg_deg": float(np.angle(eigs_sorted[0])*180/np.pi)},
    "top_30_eigs": [{"re": float(l.real), "im": float(l.imag), "abs": float(abs(l)), "arg_deg": float(np.angle(l)*180/np.pi)} for l in eigs_sorted[:30]],
    "cotaescu_m3_hits": [{"idx": i, "re": float(l.real), "im": float(l.imag), "abs": float(abs(l)), "arg_deg": float(np.angle(l)*180/np.pi)} for i, l in hits_m3],
    "cotaescu_m1_hits": [{"idx": i, "re": float(l.real), "im": float(l.imag), "abs": float(abs(l)), "arg_deg": float(np.angle(l)*180/np.pi)} for i, l in hits_m1],
    "real_half_hits": [{"idx": i, "re": float(l.real), "im": float(l.imag), "abs": float(abs(l)), "arg_deg": float(np.angle(l)*180/np.pi)} for i, l in hits_05_real],
}
with open("C:/Collatz/T_phase_twisted_full_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: C:/Collatz/T_phase_twisted_full_2026_06_01.json")
print(f"Total wall time: {time.time()-t0:.1f}s")
