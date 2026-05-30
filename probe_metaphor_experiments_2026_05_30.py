"""
probe_metaphor_experiments_2026_05_30.py

The metaphor experiments translated to math operations on the σ chain at q=17.

Key prior claim to verify: in the dominant a=b regime, σ_{m+1} = 2^-a σ_m preserves coset.
=> P(σ_m in <2>) = P(σ_1 in <2>) for all m >= 1 (modulo O(2^-136) corrections).
If true, c_inf = c(1) EXACTLY (giant rational) and not non-elementary.

Experiments:
  (A) FLICK: spectral decomposition of within-coset transfer op (via characters of <2>=Z/8).
  (B) TEMPERATURE: P_inf as function of Geom(p) parameter.
  (C) REVERSIBILITY: detailed balance test for the σ chain.
  (D) TUBE MATERIAL: exact c(m) at varying FFT level to check convergence pattern.
  (E) BOUNDARY: depth-0-to-1 kernel decomposition (dominant vs subdominant contributions).
"""
from __future__ import annotations
import numpy as np
from fractions import Fraction
import sys, math
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8

def legendre(x, q):
    x = int(x) % q
    return 0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)

# <2> = QR mod 17
H = []
x = 1
while True:
    H.append(x)
    x = (x * 2) % q
    if x == 1: break
print(f"<2> = QR mod {q} = {sorted(H)}")
H_set = set(H)
discrete_log = {h: H.index(h) for h in H}  # log_2 on <2>

# =======================================================================
# (A) FLICK: spectral decomposition of within-coset transfer operator
# =======================================================================
print("\n" + "="*70)
print("(A) FLICK: spectral decomposition of sigma chain on <2>")
print("="*70)
print("T f(sigma) = E_a[f(2^-a sigma)], a ~ Geom(1/2)")
print("Characters chi_k(sigma) = omega^(k log_2 sigma), omega = e^(2 pi i / 8)")
print("Eigenvalue: lambda_k = E[chi_k(2^-a)] = sum_a 2^-a * omega^(-a k)")
print("           = (omega^-k / 2) / (1 - omega^-k/2) = omega^-k / (2 - omega^-k)\n")
print(f"  {'k':>3} {'omega^-k':>20} {'eigenvalue lambda_k':>30} {'|lambda_k|':>12}")
for k in range(ord_2):
    omega_minus_k = np.exp(-2j * np.pi * k / ord_2)
    lam = omega_minus_k / (2 - omega_minus_k)
    print(f"  {k:>3} {f'{omega_minus_k.real:+.4f}{omega_minus_k.imag:+.4f}i':>20} "
          f"{f'{lam.real:+.6f}{lam.imag:+.6f}i':>30} {abs(lam):>12.6f}")
print()
print("  k=0 (trivial char): lambda=1 => Haar invariant (uniform on coset).")
print("  k=4 (sign char): omega^-4=-1, lambda = -1/3 EXACTLY.")
print("  Other k: |lambda| < 1, mix within coset.")
print(f"  Leading non-trivial within-coset eigenvalue: |lambda| = {1/np.sqrt(5 - 2*np.sqrt(2)):.6f}")
print("  Coincidence note: |lambda_4| = 1/3 matches the depth-collision rate 1/3 numerically,")
print("  but these are different objects (one is within-coset mixing, other is cross-depth survival).")

# =======================================================================
# (B) TEMPERATURE SWEEP: P_1 as function of p in Geom(p)
# =======================================================================
print("\n" + "="*70)
print("(B) TEMPERATURE: P_1 = P(sigma_1 in <2>) vs p in Geom(p)")
print("="*70)
print("sigma_1 in dominant case = 2^-a_1 (x-y). coset(sigma_1) = coset(x-y mod q).")
print("(x-y mod q) distribution = distribution of (X-Y mod q) for X, Y iid Tao Syracuse.")
print()
print("Tao Syracuse mod q has first digit 2^-a_1 mod q for a_1 ~ Geom(p).")
print("Vary p, compute P(X-Y in <2> mod q | X != Y mod q) for the dominant case.\n")
print(f"  {'p':>8} {'mean Geom':>10} {'P(X-Y in <2>)':>15} {'c(0)=2P-1':>10}")

inv2_mod_q = pow(2, q - 2, q)
pow_inv2_q = [pow(inv2_mod_q, k, q) for k in range(ord_2)]

for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
    # Distribution of X mod q under Geom(p): X mod q = 2^{-a} for a ~ Geom(p)
    # P(a = k) = (1-p)^(k-1) * p, k >= 1
    # P(X mod q = 2^-r mod q) for r in {1..ord_2}: sum over a with a ≡ r mod ord_2
    # = sum_{j >= 0} (1-p)^(r-1+j*ord_2) * p  for r in {1..ord_2}
    P_X = {}  # P(X mod q = value)
    for r in range(1, ord_2 + 1):
        v = pow_inv2_q[r % ord_2]
        prob = sum((1-p)**(r-1 + j*ord_2) * p for j in range(200))
        P_X[v] = P_X.get(v, 0) + prob
    # Verify probs sum to 1
    # print(f"    p={p}: P_X sums to {sum(P_X.values()):.6f}")
    # Distribution of (X - Y) mod q
    P_D = {d: 0.0 for d in range(q)}
    for vx, px in P_X.items():
        for vy, py in P_X.items():
            d = (vx - vy) % q
            P_D[d] += px * py
    p_nontrivial = sum(P_D[d] for d in range(1, q))
    p_in_H = sum(P_D[d] for d in range(1, q) if d in H_set) / p_nontrivial
    c0 = 2 * p_in_H - 1
    print(f"  {p:>8.2f} {1/p:>10.3f} {p_in_H:>15.8f} {c0:>10.6f}")

print()
print("  At p=1/2 (our case): P_0 = 73/127 = 0.5748, c(0) = 19/127 = 0.1496.")
print("  Boiling limit p->1 (a=1 always): all draws are 2^-1=9 in <2>, P->1 (full <2>).")
print("  Freezing limit p->0 (a->inf): 2^-a uniformly on <2>, P->1/2 (symmetric).")

# =======================================================================
# (C) REVERSIBILITY: detailed balance test
# =======================================================================
print("\n" + "="*70)
print("(C) REVERSIBILITY: detailed balance for the sigma chain on <2>")
print("="*70)
print("Detailed balance: pi(s) K(s, s') = pi(s') K(s', s)")
print("Our chain: s -> 2^-a s, a ~ Geom(2). K(s, 2^-a s) = 2^-a.")
print("Under Haar pi (uniform on <2>, prob 1/8 each):")
print("  pi(s) K(s, s') = (1/8) * 2^-a where 2^-a s = s'  i.e.  s' = s/2^a")
print("  pi(s') K(s', s) = (1/8) * 2^-a' where 2^-a' s' = s  i.e.  s = s'/2^a' => a' = -a mod 8")
print("  But a' must be >= 1 (Geom). For a = 1: a' = -1 mod 8 = 7. So a' = 7 corresponds to a = 1.")
print("  Detailed balance needs 2^-1 = 2^-7, i.e., 1/2 = 1/128, FALSE.")
print("  => CHAIN IS NOT REVERSIBLE. P_inf requires full transfer operator, not detailed balance.")
# Numerical check
sigma1, sigma2 = 1, 9  # 1 -> 9 via 2^-? Actually 2^-1 * 1 = 9. So a=1.
# K(1, 9) = 2^-1 = 1/2
# K(9, 1) = ?. 2^-a * 9 = 1 means 9 * 2^-a = 1 mod 17, so 2^-a = 1/9 mod 17. 1/9 = 9^-1 mod 17. 9 * 2 = 18 = 1, so 9^-1 = 2. So 2^-a = 2 means a = -1, but a >= 1, so a = -1 mod 8 = 7. K(9, 1) = 2^-7 = 1/128.
print(f"\n  Numerical: K(1, 9) = 1/2 = 0.5; K(9, 1) = 2^-7 = {1/128:.6f}")
print(f"  ratio K(s,s')/K(s',s) = {0.5 / (1/128):.1f} != 1 => not reversible.")

# =======================================================================
# (D) TUBE MATERIAL: exact c(m) convergence + decomposition
# =======================================================================
print("\n" + "="*70)
print("(D) TUBE MATERIAL: c(m) values, do they actually converge or stay constant?")
print("="*70)
print("c(m) values measured at FFT n=6:")
data = {0: 0.149606299213, 1: 0.153178230055, 2: 0.153247920779,
        3: 0.153005331692, 4: 0.152988709047, 5: 0.152988999414}
print(f"  {'m':>3} {'c(m)':>15} {'c(m)-c(1)':>15}")
for m in sorted(data):
    diff = data[m] - data[1]
    print(f"  {m:>3} {data[m]:>15.12f} {diff:>+15.6e}")
print()
print("OBSERVATION: c(m) varies by ~1e-4 for m >= 1. If my coset-preservation analysis")
print("were correct (c(m) constant for m>=1 modulo 1e-41 corrections), this variation must")
print("be a finite-n FFT artifact OR my analysis is missing something.")
print()
print("Hypothesis to check: at FFT level n=6, the c(m) at boundary m=5 is MOST biased.")
print("If true infinite-precision c(m) is constant = c(1) = 0.153178, then c(5) ~ 0.152989")
print("would be a -0.0002 bias due to truncating at n=6.")
print()
print("Test: re-extract c(2) at level n=3 (boundary) vs n=6 (well-interior). If identical,")
print("then c(m) is NOT a finite-n artifact and the variation is REAL structure.")
# Note: probe_coset_collision_op n=3 gave c(2)=0.153247920779
# probe_c_inf_depth_extrap n=6 also gives c(2)=0.153247920779
# They are IDENTICAL (verified from agent's output)
print("  c(2) at n=3 boundary: 0.153247920779")
print("  c(2) at n=6 interior: 0.153247920779")
print("  IDENTICAL => c(m) variation is REAL structure, not finite-n bias.")
print()
print("CONCLUSION: my dominant-regime coset-preservation analysis is INCOMPLETE.")
print("Something at depth m >= 2 DOES move sigma across cosets at the 1e-4 level per step.")
print("This isn't subdominant 8q^m|(a-b) (which is ~1e-41); it's something else.")

# =======================================================================
# (E) BOUNDARY: depth-0-to-1 exact P_1 computation
# =======================================================================
print("\n" + "="*70)
print("(E) BOUNDARY: depth-0-to-1 transition gives P_1 exactly")
print("="*70)
# From c(1) exact = 265011804960406635465672455997699/1730087916969634762193659498034425
c1_num = 265011804960406635465672455997699
c1_den = 1730087916969634762193659498034425
c1 = Fraction(c1_num, c1_den)
print(f"c(1) exact = (33-digit num) / (34-digit den) = {float(c1):.12f}")
P1 = (1 + c1) / 2
print(f"P_1 = (1 + c(1))/2 = {float(P1):.12f}")
print()
print(f"Note: c(0) = 19/127, so P_0 = (1 + 19/127)/2 = 146/254 = 73/127 = {73/127:.12f}")
print(f"P_1 - P_0 = {float(P1) - 73/127:+.8f}  (boundary shift)")
print()
print("c_inf (Shanks extrapolation from FFT) ~ 0.152988994")
P_inf_shanks = (1 + 0.152988994) / 2
print(f"P_inf (Shanks) ~ {P_inf_shanks:.10f}")
print(f"P_1 - P_inf = {float(P1) - P_inf_shanks:+.8f}  (asymptotic shift)")
print()
print("BOTTOM LINE: c(m) does vary across m, suggesting genuine cross-coset transitions")
print("happen at every depth (not just depth 0-to-1). This contradicts the agent's simple")
print("coset-preservation claim. The 'tube' is leakier than I claimed; some material crosses")
print("at every depth, not just the boundary. Identifying that mechanism is the open question.")
