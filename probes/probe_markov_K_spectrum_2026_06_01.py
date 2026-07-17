"""
Direct diagonalization of the Markov kernel K at q=3 levels n=1..5.

K is real => eigenvalues come in {real} or {complex conjugate pairs}.

The bilinear pair operator (which generates ε_k) has spectrum in K ⊗ K^*.
So eigenvalues of K ⊗ K^* are products lambda_i(K) * lambda_j(K)^*.

Search for:
  - rho ≈ 0.190, arg ≈ ±41°  -- ε_k Model 3 fit
    Could come from K eigenvalue directly OR from lambda_i lambda_j* = rho e^(itheta)
  - rho ≈ 0.236, arg ≈ ±15.7°  -- ε_k Model 1 fit (single complex, no asymptote)

Targets for K eigenvalue:
  Direct: |λ_K| = 0.190, arg = ±41°  →  λ_K = 0.143 ± 0.125i
  Via square: |λ_K| = sqrt(0.190) = 0.436, arg = ±20.5°  →  λ_K = 0.408 ± 0.153i
  Via 1 × λ_K: same as direct
  Via λ_K × λ_K (i.e., K^2 eigenvalue): same as via square
"""
from __future__ import annotations
import sys, os, time, json
from fractions import Fraction
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:/Collatz")

from T_lead_operator import build_markov_rational

N_MAX = 5

print(f"Diagonalizing Markov kernel K at q=3, n=1..{N_MAX}", flush=True)
print()

all_specs = {}
for n in range(1, N_MAX + 1):
    t0 = time.time()
    K_frac, coprime = build_markov_rational(n)
    K = np.array([[float(K_frac[i][j]) for j in range(len(K_frac))] for i in range(len(K_frac))], dtype=np.float64)
    dim = K.shape[0]
    eigs = np.linalg.eigvals(K)
    # Sort by |λ| descending
    idx = np.argsort(-np.abs(eigs))
    eigs_sorted = eigs[idx]
    print(f"n={n}  dim={dim}  build time={time.time()-t0:.2f}s")
    print(f"  Full spectrum (sorted by |λ|):")
    for i, lam in enumerate(eigs_sorted):
        if abs(lam) < 1e-10:
            break
        marker = ""
        if 0.18 <= abs(lam) <= 0.20 and 38 <= abs(np.angle(lam)*180/np.pi) <= 44:
            marker = " <-- MATCH Model 3 direct"
        if 0.42 <= abs(lam) <= 0.45 and 18 <= abs(np.angle(lam)*180/np.pi) <= 23:
            marker = " <-- MATCH Model 3 via square"
        if 0.23 <= abs(lam) <= 0.25 and 14 <= abs(np.angle(lam)*180/np.pi) <= 18:
            marker = " <-- MATCH Model 1"
        print(f"    {i:>3}  Re={lam.real:+.8f}  Im={lam.imag:+.8f}  |λ|={abs(lam):.6f}  arg={np.angle(lam)*180/np.pi:+8.3f}°{marker}")
    # also show top zeros count
    n_nonzero = sum(1 for lam in eigs_sorted if abs(lam) > 1e-10)
    print(f"  Nonzero eigenvalues (>1e-10): {n_nonzero} of {dim}")
    print()
    all_specs[n] = {
        "dim": dim,
        "eigvals": [{"re": float(l.real), "im": float(l.imag), "abs": float(abs(l)), "arg_deg": float(np.angle(l)*180/np.pi)} for l in eigs_sorted],
        "n_nonzero": n_nonzero,
    }

# Now: K ⊗ K* eigenvalue search
print(f"=" * 72)
print(f"Search K ⊗ K* spectrum (products λ_i · λ_j^*) at n=3")
print(f"=" * 72)
K_frac, coprime = build_markov_rational(3)
K = np.array([[float(K_frac[i][j]) for j in range(len(K_frac))] for i in range(len(K_frac))], dtype=np.float64)
eigs_K = np.linalg.eigvals(K)
eigs_K = eigs_K[np.abs(eigs_K) > 1e-10]
products = []
for i, l1 in enumerate(eigs_K):
    for j, l2 in enumerate(eigs_K):
        p = l1 * np.conj(l2)
        products.append((p, i, j, l1, l2))
products.sort(key=lambda x: -abs(x[0]))
print(f"  Top 20 K⊗K* eigenvalues by |λ|:")
print(f"  {'i':>3} {'j':>3} {'lambda':>30} {'|λ|':>12} {'arg°':>10}")
for k, (p, i, j, l1, l2) in enumerate(products[:20]):
    marker = ""
    if 0.18 <= abs(p) <= 0.20 and 38 <= abs(np.angle(p)*180/np.pi) <= 44:
        marker = " <-- MATCH Model 3"
    if 0.23 <= abs(p) <= 0.25 and 14 <= abs(np.angle(p)*180/np.pi) <= 18:
        marker = " <-- MATCH Model 1"
    print(f"  {i:>3} {j:>3} {str(p)[:30]:>30} {abs(p):>12.6f} {np.angle(p)*180/np.pi:>+10.4f}°{marker}")
print()

# Save
with open("C:/Collatz/markov_K_spectrum_2026_06_01.json", "w") as f:
    json.dump(all_specs, f, indent=2)
print("Saved: C:/Collatz/markov_K_spectrum_2026_06_01.json")
