"""
probe_aerial_dye_cameras_2026_05_30.py

Build the K=2 deepening kernel for the deep-collision shift chain at q=17. This single
object yields:
  - Aerial view: leading eigenvector pi_K2 of the 289x289 sub-stochastic matrix.
  - Animals: confirm eigenvalue lambda ~= 1/3 (sub-stochastic structure of conditioned chain).
  - Cameras: c_inf computed from pi_K2, error O(q^-2) ~ 0.003.
  - Marginalization to K=1: 17x17 mod-q transition graph for aerial structural analysis.
  - Dye: T^m delta_s for several starting states.

Setup. State s in Z/q^2 = Z/289. Each step:
  draw (a,b) ~ Geom(2)^2, draw X*, Y* ~ Syracuse offset mod q^3 (= P_3, 4913 entries).
  D = 2^-a X* - 2^-b Y* mod q^3 (need q^3 precision to extract new s' mod q^2).
  if (s + D) mod q != 0: TERMINATE with leading digit (s + D) mod q.
  else: new state s' = ((s + D) mod q^3) / q mod q^2.
"""
from __future__ import annotations
import numpy as np
import sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17; N2 = q**2; N3 = q**3
A_MAX = 60

def legendre(x, q):
    x = int(x) % q
    return 0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)

# Offset distribution mod q^3
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

print(f"Computing Syracuse offset distribution P_3 (mod {N3})...")
P3 = offset_distribution(q, 3, A_MAX)
print(f"  P3 sum = {P3.sum():.10f}, support size = {np.count_nonzero(P3)}")

print(f"\nBuilding U_dist = sum_a 2^-a * dist(2^-a X*) where X* ~ P3...")
U_dist = np.zeros(N3)
for a in range(1, A_MAX + 1):
    pow_2a = pow(2, a, N3)
    perm = (pow_2a * np.arange(N3)) % N3
    U_dist += 2.0 ** (-a) * P3[perm]
U_dist /= U_dist.sum()
print(f"  U_dist sum = {U_dist.sum():.10f}")

print(f"Computing D = U - V distribution via FFT correlation...")
F_U = np.fft.fft(U_dist)
D_dist = np.real(np.fft.ifft(F_U * np.conj(F_U)))
print(f"  D_dist sum = {D_dist.sum():.10f}")

# --- Build K=2 deepening kernel ---
print(f"\nBuilding K=2 deepening kernel ({N2}x{N2})...")
K = np.zeros((N2, N2))
p_term = np.zeros(N2)
N_chi = np.zeros(N2)  # E[chi(leading digit) | s, terminate] * p_term(s)
for s in range(N2):
    s_arr = np.arange(N3)
    sum_sd = (s + s_arr) % N3
    # deepen: sum_sd mod q == 0; else terminate with leading digit sum_sd mod q
    mod_q = sum_sd % q
    deep_mask = (mod_q == 0)
    term_mask = ~deep_mask
    # Deepening: new state = sum_sd // q mod q^2
    if deep_mask.any():
        new_s = (sum_sd[deep_mask] // q) % N2
        deltas = s_arr[deep_mask]
        np.add.at(K[s], new_s, D_dist[deltas])
    # Termination: accumulate p_term and N_chi
    if term_mask.any():
        td = s_arr[term_mask]
        digits = mod_q[term_mask]
        chi_vals = np.array([legendre(int(d), q) for d in digits])
        weights = D_dist[td]
        p_term[s] = float(np.sum(weights))
        N_chi[s] = float(np.sum(weights * chi_vals))

row_sums = K.sum(axis=1)
print(f"  Row sums (deepen probability) mean = {row_sums.mean():.6f}, std = {row_sums.std():.6f}")
print(f"  Row sums range: [{row_sums.min():.6f}, {row_sums.max():.6f}]")
print(f"  Mean termination prob: {p_term.mean():.6f}  (expected ~ 2/3)")

# --- ANIMALS: leading eigenvalue of K ---
print("\nDiagonalizing K (animals diagnostic)...")
eigvals, eigvecs = np.linalg.eig(K)
order = np.argsort(-np.abs(eigvals))
print(f"  Top 5 |eigvals|: {[f'{abs(eigvals[order[i]]):.8f}' for i in range(5)]}")
print(f"  Top eigval (complex): {eigvals[order[0]]}")
lam = float(np.real(eigvals[order[0]]))
print(f"  Leading lambda = {lam:.10f}  (compared to expected 1/3 = {1/3:.10f}, diff = {lam - 1/3:.2e})")

# Leading right eigenvector (pi_K2 quasi-stationary)
v = np.real(eigvecs[:, order[0]])
# need v sign chosen so positive
if v.sum() < 0: v = -v
pi_K2 = v / v.sum()
print(f"  pi_K2 sum = {pi_K2.sum():.6f}, min = {pi_K2.min():.4e}, max = {pi_K2.max():.4e}")
print(f"  positive entries: {(pi_K2 > 0).sum()} / {N2}")

# --- CAMERAS: c_inf via pi_K2 ---
print("\nCAMERAS: c_inf at K=2 precision")
# c_inf = sum_s pi(s) N_chi(s) / sum_s pi(s) p_term(s)
num = float(np.sum(pi_K2 * N_chi))
den = float(np.sum(pi_K2 * p_term))
c_inf_K2 = num / den
print(f"  c_inf (K=2) = {c_inf_K2:.12f}")
print(f"  c_inf (FFT n=6 Shanks) = 0.152988994")
print(f"  diff = {c_inf_K2 - 0.152988994:.2e}")

# --- AERIAL VIEW: marginalize to K=1 (mod-q) ---
print("\nAERIAL VIEW: marginalize K=2 to K=1 (mod-q transitions)")
# pi_K2 grouped by s mod q
pi_K1 = np.array([pi_K2[i::q].sum() for i in range(q)])  # sum over fiber (q lifts per mod-q class)
# Wait: states s = 0..N2-1; s mod q = s % q.
pi_K1 = np.array([pi_K2[(np.arange(N2) % q) == r].sum() for r in range(q)])
print(f"  pi_K1 (stationary on Z/{q}):")
for r in range(q):
    print(f"    s mod {q} = {r}: pi = {pi_K1[r]:.6f}  chi(r) = {legendre(r,q):+d}")
print(f"  pi_K1 symmetry check: pi(s) vs pi(-s mod {q}):")
for s in range(1, q//2 + 1):
    print(f"    pi({s}) = {pi_K1[s]:.6f},  pi({q-s}) = {pi_K1[q-s]:.6f},  match = {abs(pi_K1[s] - pi_K1[q-s]) < 1e-8}")

# Aerial view: K=1 transition matrix on Z/q
K1 = np.zeros((q, q))
for s in range(N2):
    for sp in range(N2):
        K1[s % q, sp % q] += pi_K2[s] * K[s, sp]
# Normalize rows? K1[r] is mass flowing OUT of class r to class r', conditional on starting in fiber distributed by pi_K2.
# Row sums should equal mass deepening from class r
print(f"\n  K1 row sums (deepen prob by mod-q class):")
for r in range(q):
    print(f"    row {r}: sum = {K1[r].sum():.6f}")

# --- DYE experiments: T^m delta_s for selected s, in K=2 state space ---
print("\nDYE: T^m starting from various mod-q^2 states")
for s_start in [0, 1, 17, 50, 144]:  # 0, mod 17 = 0; 1, mod 17 = 1; 17, mod 17 = 0; 50, mod 17 = 16; 144 mod 17 = 8
    e = np.zeros(N2); e[s_start] = 1
    print(f"  start s = {s_start} (mod {q} = {s_start % q}, chi = {legendre(s_start % q, q):+d}):")
    for m in (1, 2, 3, 4, 5):
        e = K @ e
        mass = e.sum()
        # marginal mod q
        marg = np.array([e[(np.arange(N2) % q) == r].sum() for r in range(q)])
        print(f"    m={m}: mass={mass:.6f}, top-5 mod-q targets: {sorted(enumerate(marg), key=lambda kv: -kv[1])[:5]}")

# --- ANIMALS confirmation ---
print(f"\nANIMALS: K is sub-stochastic, leading eigenvalue {lam:.8f} (theoretical 1/3 = {1/3:.8f})")
print(f"  diff = {lam - 1/3:.4e}")
print(f"  -> conditioned chain has survival rate 1/3 per deepening step (consistent with empirical lambda).")
