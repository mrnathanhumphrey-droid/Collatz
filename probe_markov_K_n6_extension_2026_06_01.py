"""
Extend Markov K diagonalization to n=6, look at (a) dominant subdominant |λ_2(n)|
growth pattern — extrapolate to limit, (b) presence of arg ≈ ±41° eigenvalue.

At n=3: |λ_2| = 3.19e-6
At n=4: |λ_2| = 2.95e-5  (10× growth per n)
At n=5: |λ_2| = 3.38e-4  (10× growth per n)
Extrapolation: n=6 should give ~3e-3 if linear, or saturate near 0.19 if it's the real subdominant.
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:/Collatz")

from T_lead_operator import build_markov_rational

for n in [3, 4, 5, 6]:
    t0 = time.time()
    K_frac, coprime = build_markov_rational(n)
    K = np.array([[float(K_frac[i][j]) for j in range(len(K_frac))] for i in range(len(K_frac))], dtype=np.float64)
    dim = K.shape[0]
    print(f"n={n}  dim={dim}  build time={time.time()-t0:.1f}s", flush=True)
    t1 = time.time()
    eigs = np.linalg.eigvals(K)
    idx = np.argsort(-np.abs(eigs))
    eigs_sorted = eigs[idx]
    print(f"  eig time={time.time()-t1:.1f}s", flush=True)
    # Top 15 nonzero
    print(f"  Top 15 |λ|:")
    cnt = 0
    for lam in eigs_sorted:
        if cnt >= 15:
            break
        if abs(lam) > 1e-12:
            marker = ""
            if 39 <= abs(np.angle(lam)*180/np.pi) <= 43:
                marker = " <-- arg ≈ 41°"
            if 0.18 <= abs(lam) <= 0.20:
                marker += " <-- |λ| matches Model 3"
            print(f"    Re={lam.real:+.6e}  Im={lam.imag:+.6e}  |λ|={abs(lam):.6e}  arg={np.angle(lam)*180/np.pi:+8.3f}°{marker}")
            cnt += 1
    print()
