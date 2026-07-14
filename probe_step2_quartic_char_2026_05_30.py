"""
probe_step2_quartic_char_2026_05_30.py

Step 2: Compute c_4(m) = E[χ_4((X-Y)/q^m mod q) | v_q(X-Y) = m]
where χ_4 is the quartic character mod 17 (order 4).
Test damped oscillation:
  Δ_4_m = c_4(m) - c_4(0)
  Re(Δ_4_m) and Im(Δ_4_m) each fit damped oscillation?
  Same ρ as the χ_2 (Legendre) fit?
If YES → quartic Dirichlet origin confirmed.

Also report the full c_4(m) complex sequence for m=0..5 at FFT precision.

χ_4 table (using g=3 as primitive root of 17, χ_4(3^k) = i^k):
  +1 on {1, 4, 13, 16}    (kernel: 4th-power residues)
  +i on {3, 5, 12, 14}
  -1 on {2, 8, 9, 15}     (quadratic-residue-not-quartic)
  -i on {6, 7, 10, 11}
Note χ_4² = χ_2 (Legendre).
"""
from __future__ import annotations
import sys, gc, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 17

# Build χ_4 table using g=3 as primitive root
g = 3
chi4 = np.zeros(q, dtype=np.complex128)
chi4[0] = 0
x = 1
for k in range(q - 1):
    chi4[x] = 1j ** k   # i^k
    x = (x * g) % q

print(f"χ_4 table (g={g} primitive root):")
for x in range(q):
    print(f"  χ_4({x:2}) = {chi4[x]}")

# Verify χ_4^2 = Legendre
def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

chi2_check = np.array([(chi4[x]**2).real for x in range(q)])
chi2_true = np.array([legendre(x, q) for x in range(q)], dtype=float)
print(f"\nχ_4² = Legendre? {np.allclose(chi2_check, chi2_true)}")
# Re part is the real of (chi4)^2 which should match Legendre
# (chi4)^2 takes values in {1, -1} since i^2=-1, (-1)^2=1, (-i)^2=-1
# Verify
print(f"  (Im part of chi4^2): {[(chi4[x]**2).imag for x in range(q)]}")  # should be all 0

# === FFT machinery ===
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

def c_chi_m(n, m, chi_table, A_MAX=200):
    """Compute E[chi((D)/q^m mod q) | v_q(D)=m] at FFT level n with complex chi."""
    N = q ** n
    P = offset_distribution(q, n, A_MAX)
    mu = np.fft.fft(P)
    a2 = np.abs(mu) ** 2
    PD = np.real(np.fft.ifft(a2))
    qm = q ** m
    qmp1 = q ** (m + 1)
    d_vals = np.arange(N)
    mask = (d_vals % qm == 0) & (d_vals % qmp1 != 0)
    leading_digit = (d_vals[mask] // qm) % q
    chi_vals = chi_table[leading_digit]
    num = (PD[mask] * chi_vals).sum()
    den = PD[mask].sum()
    return complex(num / den)

# Compute c_4(m) for m=0..5
print(f"\n=== Computing c_4(m) at A_MAX=200 ===")
c4_vals = {}
for m in (0, 1, 2, 3, 4, 5):
    n = max(m + 1, 2)
    t0 = time.time()
    c4_vals[m] = c_chi_m(n=n, m=m, chi_table=chi4, A_MAX=200)
    print(f"  c_4({m}) (n={n}) = {c4_vals[m]:.12f}   ({time.time()-t0:.1f}s)", flush=True)
    gc.collect()

# Note: c_4(0) = E[χ_4((X-Y) mod q) | X != Y mod q].
# Probably nonzero complex. Use it as the anchor.
c4_0 = c4_vals[0]
print(f"\nc_4(0) = {c4_0}")
print(f"|c_4(0)| = {abs(c4_0):.12f}")
print(f"arg(c_4(0)) = {np.angle(c4_0):.6f} rad = {np.degrees(np.angle(c4_0)):.4f}°")

# Δ_4_m = c_4(m) - c_4(0), complex
print(f"\n=== Δ_4_m sequence (complex) ===")
delta4 = {m: c4_vals[m] - c4_0 for m in c4_vals}
for m in sorted(delta4):
    d = delta4[m]
    print(f"  Δ_4_{m} = {d.real:+.12e}{d.imag:+.12e}j  |Δ|={abs(d):.6e}  arg={np.degrees(np.angle(d)):+.4f}°")

# === Fit damped oscillation in COMPLEX form ===
# Δ_4_m = A_4 + C·z^m + D·z̄^m  (5 complex parameters, but z̄ = conj of z so really 4)
# OR assume z FROM χ_2 fit: z = 0.0334 + 0.0683i and just fit A_4 + C·z^m + D·z̄^m
# OR fit independently: Δ_4_m satisfies same recurrence
#   Δ_4_{m+2} = u·Δ_4_{m+1} - v·Δ_4_m + K_4
# with same (u, v) as χ_2 fit if same dynamics, K_4 different.

# First: extract (u, v, K_4) from χ_4 data independently
# 3 equations, 3 unknowns (complex):
# Δ_4_3 = u·Δ_4_2 - v·Δ_4_1 + K_4
# Δ_4_4 = u·Δ_4_3 - v·Δ_4_2 + K_4
# Δ_4_5 = u·Δ_4_4 - v·Δ_4_3 + K_4
M4 = np.array([
    [delta4[2], -delta4[1], 1.0],
    [delta4[3], -delta4[2], 1.0],
    [delta4[4], -delta4[3], 1.0],
], dtype=np.complex128)
rhs4 = np.array([delta4[3], delta4[4], delta4[5]], dtype=np.complex128)
u_4, v_4, K_4 = np.linalg.solve(M4, rhs4)
print(f"\n=== χ_4 recurrence fit ===")
print(f"  u_4 = {u_4}  (should match χ_2 u = 0.0668618 if same dynamics)")
print(f"  v_4 = {v_4}")
print(f"  K_4 = {K_4}")
print(f"\n  χ_2 fit (reference): u = 0.0668618, v = 0.00577844")

# If u_4 ≈ 0.0668618 and v_4 ≈ 0.00577844 (real), same dynamics.
diff_u = u_4 - 0.0668617916014197
diff_v = v_4 - 5.7784400779138535e-03
print(f"\n  u_4 - u_χ2 = {diff_u}")
print(f"  v_4 - v_χ2 = {diff_v}")

# Asymptote
denom = 1 - u_4 + v_4
A_4 = K_4 / denom
print(f"\n  A_4 (asymptote) = {A_4}")
print(f"  c_4_∞ from model = c_4(0) + A_4 = {c4_0 + A_4}")

# Verify model on all 6 points
print(f"\n=== Model verification ===")
for m in (1, 2, 3, 4, 5):
    if m+2 <= 5:
        predicted_m2 = u_4 * delta4[m+1] - v_4 * delta4[m] + K_4
        actual = delta4[m+2]
        print(f"  Δ_4_{m+2}: actual={actual.real:+.10e}{actual.imag:+.10e}j  "
              f"model={predicted_m2.real:+.10e}{predicted_m2.imag:+.10e}j  "
              f"diff={abs(actual-predicted_m2):.3e}")

# Test: if u_4, v_4 are REAL and match χ_2 values, then same dynamics with different observable.
# In that case, the (u, v) eigenvalues z = ρ·e^{iθ} are the SAME, and the quartic-vs-quadratic
# differs only in amplitude C (complex).
print(f"\n=== Same-z hypothesis test ===")
u_real = u_4.real; u_imag = u_4.imag
v_real = v_4.real; v_imag = v_4.imag
print(f"  u_4 imag part: {u_imag:.6e}  (should be ~0 if same dynamics)")
print(f"  v_4 imag part: {v_imag:.6e}  (should be ~0 if same dynamics)")
ratio_u = u_real / 0.0668617916014197
ratio_v = v_real / 5.7784400779138535e-03
print(f"  u_4 real / u_χ2 = {ratio_u:.8f}  (should be ~1 if same dynamics)")
print(f"  v_4 real / v_χ2 = {ratio_v:.8f}  (should be ~1 if same dynamics)")
