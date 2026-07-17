"""
w4_diagnostic_M4.py — definitive test of W4's "1/√3 spectral" vs "counting tautology" reading.

Computes M_k(η) and ‖d_k‖² = Σ_{η ≠ 1} |M_k(η)|² for k = 1, 2, 3, 4.

Then computes level-to-level ratios in three forms:
1. Raw L²: ‖d_{k+1}‖ / ‖d_k‖
2. Normalized by √(mode count): (‖d_{k+1}‖ / √n_{k+1}) / (‖d_k‖ / √n_k)
3. Per-mode: ‖d_{k+1}‖² / ‖d_k‖² (geometric per-level decay)

Reads:
- If raw ‖d_k‖ stays bounded or decays → spectral reading (W4's claim survives)
- If raw ‖d_k‖ grows ~ √n_k (i.e., proportional to √mode-count) → counting tautology (W4 dies)

Per auditor's "deadlock-breaker #1 + #2" recommendation in W4_ADVERSARIAL_AUDIT.md.
"""
import sys
import os
sys.path.insert(0, r'C:\Collatz')
sys.stdout.reconfigure(encoding='utf-8')

from bilinear_pair_operator import build_markov_rational, stationary_rational, compute_M_n

print("# W4 deadlock-breaker: M_k(η) for k=1..4")
print()

records = []

for k in [1, 2, 3, 4]:
    print(f"## Computing k={k} ...")
    K, coprime = build_markov_rational(k)
    pi = stationary_rational(K)
    M_dict = compute_M_n(pi, coprime, k)

    S_k = M_dict[1].real

    # ‖d_k‖² = Σ_{η ≠ 1, η coprime to 3} |M_k(η)|²
    d_norm_sq = sum(abs(M_dict[eta])**2 for eta in coprime if eta != 1)
    d_norm = d_norm_sq**0.5

    # Number of off-diagonal modes
    n_modes = len(coprime) - 1

    # Normalized
    if n_modes > 0:
        normalized = d_norm / (n_modes**0.5)
    else:
        normalized = 0.0

    records.append({
        'k': k,
        'S_k': S_k,
        'd_norm_sq': d_norm_sq,
        'd_norm': d_norm,
        'n_modes': n_modes,
        'normalized': normalized,
    })

    print(f"  S_k = {S_k:.8f}")
    print(f"  ‖d_k‖² = {d_norm_sq:.8f}")
    print(f"  ‖d_k‖ = {d_norm:.8f}")
    print(f"  n_modes (= |(Z/3^k)*| - 1) = {n_modes}")
    print(f"  normalized = ‖d_k‖ / √n_modes = {normalized:.8f}")
    print()

print("# Level-to-level ratios")
print()
print(f"{'k→k+1':>8s} {'raw ratio ‖d_{k+1}‖/‖d_k‖':>30s} {'normalized ratio':>22s} {'n_{k+1}/n_k':>14s}")
for i in range(len(records) - 1):
    r1 = records[i]
    r2 = records[i+1]
    raw_ratio = r2['d_norm'] / r1['d_norm']
    norm_ratio = r2['normalized'] / r1['normalized']
    n_ratio = r2['n_modes'] / r1['n_modes']
    print(f"{r1['k']}→{r2['k']:1d}     {raw_ratio:>30.6f} {norm_ratio:>22.6f} {n_ratio:>14.6f}")

print()
print("# Target comparisons")
print(f"  1/√3 = {(1/3)**0.5:.6f}")
print(f"  √3   = {3**0.5:.6f}")
print(f"  Mode-count growth √(n_{{k+1}}/n_k):  ~√3 = {3**0.5:.6f} (since |(Z/3^k)*| = 2·3^(k-1))")
print()

print("# Interpretation guide:")
print("  - If raw ratio ≈ √3 ≈ 1.732 and normalized ratio ≈ 1 → counting tautology (W4 dies)")
print("  - If raw ratio ≈ 1 (bounded) and normalized ratio ≈ 1/√3 → consistent spectral decay (W4 lives)")
print("  - If raw ratio < 1 → genuine spectral decay (W4 strongly lives)")
print("  - If raw ratio > √3 → mode count + amplitude both grow (W4 strongly dies)")
print()

# Save JSON
import json
out = {
    'description': 'W4 deadlock-breaker M_k(η) computation for k=1..4',
    'records': records,
    'ratios': [
        {
            'from_k': records[i]['k'],
            'to_k': records[i+1]['k'],
            'raw_ratio': records[i+1]['d_norm'] / records[i]['d_norm'],
            'normalized_ratio': records[i+1]['normalized'] / records[i]['normalized'],
            'mode_count_ratio': records[i+1]['n_modes'] / records[i]['n_modes'],
        }
        for i in range(len(records) - 1)
    ],
    'targets': {
        '1/sqrt(3)': (1/3)**0.5,
        'sqrt(3)': 3**0.5,
    }
}
out_path = r'C:\Collatz\experiments_output\w4_diagnostic_M4.json'
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2, default=lambda o: float(o) if hasattr(o, 'real') else str(o))
print(f"[save] {out_path}")
