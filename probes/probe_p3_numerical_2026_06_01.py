"""
Compute p_3(σ) distribution NUMERICALLY via offset_distribution at N = q^4.

Same algorithm as probe_c5_highprec but with n=4 (so N = 17^4 = 83,521).
Trivial memory, runs in seconds.

Extract full distribution P_X mod q^4, then:
  - c(3) = Σ_j χ_L(j) · P(D mod q^4 = j·q^3) / Σ_j P(D mod q^4 = j·q^3)
  - p_3(σ) = P(D mod q^4 = σ·q^3) / Σ_σ P(D mod q^4 = σ·q^3)
  - |<ω, p_3>|^2 and the diff from c_inf

This is the decisive depth-3 test for |<ω, μ_∞>|^2 = c_∞ conjecture.
"""
from __future__ import annotations
import sys, gc, time
from mpmath import mp, mpf, mpc
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30
q = 17
n = 4   # for p_3 we need D mod q^(3+1) = q^4
N = q ** n   # 83521
A_MAX = 100  # 2^-200 tail

print(f"Computing p_3 distribution via offset_distribution: q={q}, n={n}, N={N}, A_MAX={A_MAX}, dps={mp.dps}")

def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

chi_L_table = [legendre(x, q) for x in range(q)]

# omega = chi_4 (order 4 character)
omega_re = {1: 1, 2: -1, 4: 1, 8: -1, 9: -1, 13: 1, 15: -1, 16: 1}
omega_im = {3: 1, 5: 1, 6: -1, 7: -1, 10: -1, 11: -1, 12: 1, 14: 1}

def offset_distribution_mp(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    P_U = [mpf(0)] * N
    p_val = inv2
    half = mpf(1) / mpf(2)
    cur_w = half
    for a in range(1, A_MAX + 1):
        P_U[p_val] += cur_w
        p_val = (p_val * inv2) % N
        cur_w *= half
    total = sum(P_U[i] for i in range(N))
    P_U = [x / total for x in P_U]

    u_sup = [i for i in range(N) if P_U[i] != 0]
    u_wt  = [P_U[i] for i in u_sup]
    print(f"    |u_sup| = {len(u_sup)}", flush=True)

    P_S = P_U
    for j in range(n - 1, 0, -1):
        t_j = time.time()
        P_V = [mpf(0)] * N
        for i in range(N):
            P_V[(1 + q * i) % N] += P_S[i]
        P_new = [mpf(0)] * N
        for u, w in zip(u_sup, u_wt):
            for i in range(N):
                P_new[(u * i) % N] += w * P_V[i]
        P_S = P_new
        del P_V
        gc.collect()
        print(f"    j={j} done in {time.time()-t_j:.1f}s", flush=True)
    return P_S

t0 = time.time()
P_X = offset_distribution_mp(q, n, A_MAX)
print(f"P_X done in {time.time()-t0:.1f}s", flush=True)

# Extract p_3(σ): for σ ∈ 1..q-1, sum P(D mod N = σ·q^3 + something with all lower digits zero)
# Actually we need P(d_0=d_1=d_2=0 AND d_3=σ). D mod q^4 must equal σ·q^3 exactly.
qm = q ** 3  # 4913

t1 = time.time()
p3_sigma = [mpf(0)] * q  # σ = 0..16
for j in range(1, q):  # σ = 1..16
    d = j * qm
    Pd = mpf(0)
    for y in range(N):
        Pd += P_X[y] * P_X[(y - d) % N]
    p3_sigma[j] = Pd
print(f"P_D extraction in {time.time()-t1:.1f}s", flush=True)

# Normalize
total_mass = sum(p3_sigma)
p3 = [v / total_mass for v in p3_sigma]
print(f"\n=== p_3(σ) distribution ===")
for sig in range(q):
    chi_val = chi_L_table[sig]
    print(f"  σ={sig:2d}  p_3 = {float(p3[sig]):.15f}  χ_L={chi_val:+d}")

# c(3) = Σ χ_L(σ) p_3(σ)
c3 = sum(mpf(chi_L_table[sig]) * p3[sig] for sig in range(q))
print(f"\nc(3) numerical = {c3}")
print(f"c(3) reference ≈ 0.15300533169151")
print(f"c(3) diff      = {float(c3 - mpf('0.1530053316915140267018082388619123587748')):+.2e}")

# |<omega, p_3>|^2
Re = sum(omega_re[s] * p3[s] for s in range(1, q) if s in omega_re)
Im = sum(omega_im[s] * p3[s] for s in range(1, q) if s in omega_im)
mag_sq = Re * Re + Im * Im

print(f"\n=== |<omega, p_3>|^2 ===")
print(f"  Re<omega, p_3> = {Re}")
print(f"  Im<omega, p_3> = {Im}")
print(f"  |<omega, p_3>|^2 = {mag_sq}")

c_inf_ref = mpf("0.15298912060588517527891674877413229926086222622334")
diff = mag_sq - c_inf_ref
print(f"\n  c_inf reference = {c_inf_ref}")
print(f"  |<omega, p_3>|^2 - c_inf = {float(diff):+.6e}")
print(f"\n  Predicted if 30× decay continues: ≈ {2.27e-6/30:+.1e}")
print(f"  Predicted plateau if coincidence:  ≈ +2.27e-6")
print(f"  Ratio to previous diff: {float(diff/mpf('2.27e-6')):.4f}")

print(f"\nTotal time: {time.time()-t0:.1f}s")

# Save
import json
out = {
    "q": q, "n": n, "N": N, "dps": mp.dps,
    "p_3_numerical": {str(sig): float(p3[sig]) for sig in range(q)},
    "c_3": float(c3),
    "Re_omega_p3": float(Re),
    "Im_omega_p3": float(Im),
    "omega_inner_sq": float(mag_sq),
    "diff_from_c_inf": float(diff),
}
with open("C:/Collatz/p3_numerical_2026_06_01.json", "w") as f:
    json.dump(out, f, indent=2)
print("Saved to C:/Collatz/p3_numerical_2026_06_01.json")
