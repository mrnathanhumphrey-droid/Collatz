"""
T-A2b' (NEW): q=41, ≡ 1 mod 8, splits in Z[i] as 41 = (5+4i)(5-4i).
Second confirmation of the q ≡ 1 mod 8 slow-decay prediction (after q=17).

Predict:
- c(0) ≠ 0 (clean rational)
- c_∞ ≠ 0 (non-trivial limit, regulator-class)
- Damped oscillation with ρ small (analog of q=17's 0.076)
- Eigenvalue direction along Gaussian prime (5+4i) or (5-4i)

Limit: 41^5 = 115M (10 min), 41^6 = 4.75B (too big). So c(0)..c(4) feasible.
"""
from __future__ import annotations
import sys, gc, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 41
def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

print(f"q = {q}, q mod 4 = {q % 4}, q mod 8 = {q % 8}")
print(f"Legendre(-1/{q}) = {legendre(-1, q)} (predicted +1 for q ≡ 1 mod 4)")
print(f"Legendre(2/{q})  = {legendre(2, q)} (predicted +1 for q ≡ 1 mod 8)")
ord_2 = 1; x = 2 % q
while x != 1: ord_2 += 1; x = (x * 2) % q
print(f"ord_{q}(2) = {ord_2}, index in (Z/{q})* = {(q-1)//ord_2}")
print(f"41 = 25 + 16 = (5+4i)(5-4i) in Z[i]")

chi_q = np.array([legendre(x, q) for x in range(q)], dtype=np.float64)

def offset_distribution(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a); p = (p * inv2) % N
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

def c_at_n_m(n, m, A_MAX=200):
    N = q ** n
    P = offset_distribution(q, n, A_MAX)
    mu = np.fft.fft(P)
    PD = np.real(np.fft.ifft(np.abs(mu)**2))
    qm = q ** m; qmp1 = q ** (m + 1)
    d_vals = np.arange(N)
    mask = (d_vals % qm == 0) & (d_vals % qmp1 != 0)
    leading_digit = (d_vals[mask] // qm) % q
    num = (PD[mask] * chi_q[leading_digit]).sum()
    den = PD[mask].sum()
    return float(num / den)

print(f"\nComputing c(m) at q={q} for m=0..4...")
results = {}
for m in range(0, 5):
    n = max(m + 1, 2)
    t0 = time.time()
    print(f"  m={m}, n={n}, N={q**n}...", flush=True)
    try:
        results[m] = c_at_n_m(n=n, m=m, A_MAX=200)
        print(f"    c({m}) = {results[m]:+.16f}  ({time.time()-t0:.1f}s)")
    except (MemoryError, ValueError) as e:
        print(f"    fail: {e}"); break
    gc.collect()

print(f"\n=== q={q} results ===")
for m in sorted(results):
    print(f"  c({m}) = {results[m]:+.16f}")

# Damping behavior — should be SLOW (not -3/5), small ρ analogous to q=17's 0.076
print(f"\n=== Successive differences c(m+1) - c(m) ===")
for m in range(4):
    if m+1 in results:
        d = results[m+1] - results[m]
        print(f"  c({m+1}) - c({m}) = {d:+.6e}")

print(f"\n=== Damping ratios c(m+1)/c(m) ===")
for m in range(4):
    if m+1 in results and abs(results[m]) > 1e-10:
        r = results[m+1] / results[m]
        print(f"  c({m+1})/c({m}) = {r:+.6f}  (predicted near +1 with small oscillation; NOT -3/5)")

# Try the recurrence fit
if 4 in results:
    deltas = [results[m] - results[0] for m in range(5)]
    try:
        M = np.array([
            [deltas[2], -deltas[1], 1.0],
            [deltas[3], -deltas[2], 1.0],
            [deltas[4], -deltas[3], 1.0],
        ])
        rhs = np.array([deltas[3], deltas[4], deltas[4]])  # only 5 points; underdetermined
        # try recurrence using only 4 points
        M = np.array([
            [deltas[1], -1.0, 1.0],
            [deltas[2], -deltas[1], 1.0],
            [deltas[3], -deltas[2], 1.0],
        ])
        # alternative: solve damped osc system from deltas
        M2 = np.array([
            [deltas[2], -deltas[1], 1.0],
            [deltas[3], -deltas[2], 1.0],
        ])
        # underdetermined; just report eigenvalue heuristically
        # Compute eigenvalue z from c(m) - c_∞ as c(m+2) = u·c(m+1) - v·c(m) + K
        if 4 in results:
            # 3 eq, 3 unknowns
            M = np.array([
                [results[1], -results[0], 1.0],
                [results[2], -results[1], 1.0],
                [results[3], -results[2], 1.0],
            ])
            rhs = np.array([results[2], results[3], results[4]])
            u, v, K = np.linalg.solve(M, rhs)
            print(f"\n=== Recurrence fit ===")
            print(f"  u = {u:.8f}, v = {v:.10f}, K = {K:.10f}")
            if v > 0:
                rho = np.sqrt(v)
                print(f"  ρ = {rho:.6f}  (predicted small, NOT 3/5; analog of q=17's 0.076)")
                if abs(u / (2*rho)) <= 1:
                    theta = np.arccos(u / (2*rho))
                    print(f"  θ = {theta:.6f} rad = {np.degrees(theta):.4f}°")
                    z = rho * np.exp(1j * theta)
                    print(f"  eigenvalue z = {z}")
                    # Test alignment with (5+4i)/(some scale)
                    # (5+4i) has |.|=sqrt(41), arg=arctan(4/5)=38.66°
                    print(f"  arctan(4/5) = 38.66° (predicted direction along (5+4i) Gaussian prime)")
                    print(f"  z direction: arg = {np.degrees(theta):.4f}° vs Gaussian prime (5+4i) arg = 38.66°")
            else:
                print(f"  v = {v:.6e} < 0 (real eigenvalues)")
                disc = u**2 - 4*v
                if disc > 0:
                    z1 = (u + np.sqrt(disc))/2; z2 = (u - np.sqrt(disc))/2
                    print(f"  z1 = {z1}, z2 = {z2}")
    except Exception as e:
        print(f"  recurrence fit failed: {e}")

print(f"\n=== Verdict ===")
if 0 in results and abs(results[0]) > 1e-10:
    print(f"  c(0) ≠ 0 ✓ (={results[0]:.6f})")
    if 4 in results and abs(results[4]) > 0.5 * abs(results[0]):
        print(f"  c(m) not collapsing to 0 ✓ (slow decay consistent with q ≡ 1 mod 8)")
        print(f"  *** T-A2b' PASS: q=41 shows slow-decay regime ***")
    elif 4 in results and abs(results[4]) < 0.1 * abs(results[0]):
        print(f"  c(m) collapsing (FAST decay) — INCONSISTENT with q ≡ 1 mod 8 prediction")
        print(f"  *** T-A2b' FAIL ***")
    else:
        print(f"  *** T-A2b' AMBIGUOUS — need more points ***")
else:
    print(f"  *** T-A2b' FAIL: c(0) = 0 unexpectedly ***")
