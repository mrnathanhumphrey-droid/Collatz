"""
probe_sigma_chain_matrix_2026_05_30.py

Build the σ-chain transition matrix EMPIRICALLY from FFT data:
  T_m(s, s') = P(σ_{m+1} = s' | σ_m = s, v_q(D) >= m+1)
extracted by reading consecutive q-adic digits of P_D at level n ≥ m+2.

For each starting state s ∈ {1,...,16}, sum P_D over d with leading non-zero
digit s at q-position m and next digit anything at q-position m+1.

Then:
  - Check if T_1, T_2, T_3 are nearly identical (asymptotic Markov).
  - Compute stationary π of T_∞ (left null vector of T - I).
  - c_∞ = sum_s π(s) χ_2(s).

If the chain IS Markov, this gives c_∞ at ~10-digit precision from float64 FFT.
"""
from __future__ import annotations
import sys, gc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 17

def offset_distribution(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a)
        p = (p * inv2) % N
    P_U /= P_U.sum()
    u_support = np.nonzero(P_U)[0]; u_weight = P_U[u_support]
    P_S = P_U.copy()
    for j in range(n - 1, 0, -1):
        v_idx = (1 + q * arange) % N
        P_V = np.zeros(N, dtype=np.float64)
        np.add.at(P_V, v_idx, P_S)
        P_new = np.zeros(N, dtype=np.float64)
        for u, w in zip(u_support, u_weight):
            P_new[(u * arange) % N] += w * P_V
        P_S = P_new
    return P_S

def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

chi_q = np.array([legendre(x, q) for x in range(q)], dtype=np.float64)

print(f"Computing P_X at n=5 (length {q**5}, A_MAX=200)...")
import time
t0 = time.time()
n = 5
N = q ** n
P_X = offset_distribution(q, n, 200)
print(f"  done in {time.time()-t0:.1f}s")

# Compute P_D via FFT
mu = np.fft.fft(P_X)
PD = np.real(np.fft.ifft(np.abs(mu)**2))
del mu
gc.collect()

# Extract joint P(σ_m, σ_{m+1}) for m=1,2,3 from PD.
# σ_m is the leading non-zero q-adic digit of D at q-position m.
# For v_q(D) = m, D = s_m q^m + s_{m+1} q^{m+1} + (higher) q^{m+2} + ...
# For σ_m = s_m to be non-zero, s_m ∈ {1..16}.
# σ_{m+1} can be anything in {0..16}.

# At FFT level n=5, d ∈ {0, ..., 17^5 - 1}. The (s_m, s_{m+1}) joint:
# P(σ_m = s_m, σ_{m+1} = s_{m+1}) (conditional on v_q(D) ≥ m):
#   = sum over d ≡ s_m·q^m + s_{m+1}·q^{m+1} mod q^{m+2} with d != 0 mod q^m·... oh wait
# Actually conditional on v_q(D) = m exactly (σ_m ≠ 0): just need the (s_m, s_{m+1}) digits.
# Sum P_D(d) over all d with: d mod q^m == 0, (d / q^m) mod q == s_m, (d / q^{m+1}) mod q == s_{m+1}

def joint_sigma_m_m1(m, q, n, PD):
    """Compute P(σ_m, σ_{m+1}) as a (q-1) x q matrix (s_m ∈ {1..q-1}, s_{m+1} ∈ {0..q-1})."""
    N = q ** n
    qm = q ** m
    qmp1 = q ** (m + 1)
    qmp2 = q ** (m + 2) if m + 2 <= n else N  # cap at N
    d_vals = np.arange(N)
    # Condition: d % qm == 0 (so v_q >= m)
    mask_vqgem = (d_vals % qm == 0)
    # Extract s_m and s_{m+1}
    s_m = (d_vals // qm) % q       # leading digit at q-pos m
    s_mp1 = (d_vals // qmp1) % q   # next digit

    joint = np.zeros((q-1, q))  # joint[s_m-1, s_{m+1}]
    for sm in range(1, q):  # s_m must be non-zero for σ_m = s_m
        for smp1 in range(0, q):
            sel = mask_vqgem & (s_m == sm) & (s_mp1 == smp1)
            joint[sm-1, smp1] = PD[sel].sum()
    return joint

# Compute joint for m=1, 2, 3
results = {}
for m in (1, 2, 3):
    print(f"\nComputing joint P(σ_{m}, σ_{m+1})...")
    joint = joint_sigma_m_m1(m, q, n, PD)
    # Marginal P(σ_m)
    P_sm = joint.sum(axis=1)
    P_smp1 = joint.sum(axis=0)
    # T_m[s, s'] = joint[s-1, s'] / P_sm[s-1]
    T = joint / P_sm[:, None]
    results[m] = {'joint': joint, 'P_sm': P_sm, 'P_smp1': P_smp1, 'T': T}
    # Validate: row sums of T should be 1
    print(f"  row sums of T_{m}: max dev from 1 = {np.abs(T.sum(axis=1) - 1).max():.2e}")
    # c(m) from marginal: sum chi(s) * P(s)
    c_m = sum(chi_q[s] * P_sm[s-1] for s in range(1, q))
    # but normalize by total mass at v_q = m exactly
    total_mass = P_sm.sum()
    c_m /= total_mass
    print(f"  c({m}) from marginal: {c_m:.12f}")

# === Compare T_1 vs T_2 vs T_3 ===
print(f"\n=== T_m matrix comparison (T_m vs T_3 as 'asymptotic') ===")
T1 = results[1]['T']
T2 = results[2]['T']
T3 = results[3]['T']
print(f"  ||T_1 - T_2|| (Frobenius) = {np.linalg.norm(T1 - T2):.6e}")
print(f"  ||T_2 - T_3|| (Frobenius) = {np.linalg.norm(T2 - T3):.6e}")
print(f"  ||T_1 - T_3|| (Frobenius) = {np.linalg.norm(T1 - T3):.6e}")

# === Use T_3 as the asymptotic T ===
# But T has shape (q-1, q) — extends to (q-1, q) including σ' = 0 (terminating).
# For Markov chain on {1..q-1}, restrict to (q-1, q-1) — drop σ' = 0 column.
# Or σ' = 0 means deepening terminates (no further depth). Probably rare.
T_full = T3
print(f"\nT_3[s, 0] (terminating probability per row): {T_full[:, 0]}")
print(f"  Max P(σ_{{m+1}} = 0 | σ_m = s) = {T_full[:, 0].max():.6e}")
# If small, can ignore. Otherwise need to handle.

# Build Markov chain on {1..q-1} ignoring s'=0 termination
T_markov = T_full[:, 1:q]  # shape (q-1, q-1)
# Renormalize rows to sum to 1
row_sums = T_markov.sum(axis=1)
T_markov_renorm = T_markov / row_sums[:, None]
print(f"\nT_markov row sums (pre-renorm, deviation from 1): {1 - row_sums}")

# Find left stationary: π·T = π
# Equivalently, eigenvector of T^T for eigenvalue 1
evals, evecs = np.linalg.eig(T_markov_renorm.T)
print(f"\nT_markov eigenvalues:")
for i, ev in enumerate(evals):
    print(f"  λ_{i}: {ev}  |λ|={abs(ev):.6f}")

# The eigenvalue closest to 1 is the stationary one
idx_one = np.argmin(np.abs(evals - 1))
pi_raw = evecs[:, idx_one]
# Normalize: π is a probability vector, sum = 1
pi = np.real(pi_raw / pi_raw.sum())
print(f"\nStationary π (sum=1):")
for s in range(1, q):
    print(f"  π({s}) = {pi[s-1]:+.10f}")

# c_inf from stationary
c_inf_stationary = sum(pi[s-1] * chi_q[s] for s in range(1, q))
print(f"\nc_∞ from T_3 stationary = {c_inf_stationary:.12f}")
print(f"c_∞ from damped osc model = 0.152989120606")
print(f"c_∞ Shanks 3,4,5         = 0.152988994428")
print(f"  diff (stationary - Shanks) = {c_inf_stationary - 0.152988994428:+.6e}")

# Also: coset probabilities
in_QR = np.array([1 if legendre(s, q) == 1 else 0 for s in range(1, q)])
P_QR = (pi * in_QR).sum()
P_NQR = (pi * (1 - in_QR)).sum()
print(f"\n  P(σ_∞ in QR)  = {P_QR:.10f}  (=(1 + c_∞)/2 = {(1 + c_inf_stationary)/2:.10f})")
print(f"  P(σ_∞ in NQR) = {P_NQR:.10f}")

# Within-coset distribution
print(f"\nWithin-coset distribution (relative):")
qr_states = [s for s in range(1, q) if legendre(s, q) == 1]
nqr_states = [s for s in range(1, q) if legendre(s, q) == -1]
print(f"  QR states: {qr_states}")
print(f"  P(σ_∞ = s | σ_∞ ∈ QR) (uniform = {1/8:.6f}):")
for s in qr_states:
    print(f"    s={s}: {pi[s-1] / P_QR:.10f}")
print(f"  NQR states: {nqr_states}")
print(f"  P(σ_∞ = s | σ_∞ ∈ NQR):")
for s in nqr_states:
    print(f"    s={s}: {pi[s-1] / P_NQR:.10f}")
