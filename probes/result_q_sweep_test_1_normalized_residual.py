"""
result_q_sweep_test_1_normalized_residual.py
Followup that tests whether the q-UNIVERSAL rate-1/2 question survives at
the level of the normalized residual

    epsT_n^{(q)} := T_n^{(q)} - T(q),    where T_n^{(q)} := S_n^{(q)} * (3/q)^n.

For q=3, T_n = S_n and T(q) = S_inf, so this reduces to the original
ε_n test (rate 1/2 confirmed).

For q >= 5, T_n -> T(q) is finite by Test 1's structural finding. The
question: does |epsT_n| * 2^n stay bounded across q?

Two diagnostics, both independent of T(q) extrapolation:
  (A) ratio of successive |T_n - T_{n-1}| differences = empirical rate
  (B) for q with >= 3 T values, Aitken Delta^2 estimate of T(q) and check
      |T_n - T(q)| * 2^n envelope.
"""
from __future__ import annotations
import csv, sys
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

ENV = r"C:\Collatz\result_q_sweep_test_1_envelope.csv"

S_data = {}
with open(ENV) as f:
    for row in csv.DictReader(f):
        q = int(row["q"]); n = int(row["n"])
        S_data.setdefault(q, {})[n] = Fraction(int(row["S_n_num"]),
                                                int(row["S_n_den"]))

# --------- Build T_n^{(q)} = S_n * (3/q)^n exactly ---------
T_data = {q: {n: S_data[q][n] * Fraction(3, q) ** n
              for n in S_data[q]} for q in S_data}

print("=" * 78)
print("T_n^{(q)} = S_n^{(q)} * (3/q)^n  (exact; decimal shown)")
print("=" * 78)
for q in sorted(T_data):
    print(f"  q={q}:  " + ", ".join(
        f"T_{n}={float(T_data[q][n]):.6f}" for n in sorted(T_data[q])))
print()

# --------- Diagnostic A: |T_n - T_{n-1}| ratios ---------
print("=" * 78)
print("Diagnostic A: empirical rate via |T_{n+1} - T_n| / |T_n - T_{n-1}|")
print("=" * 78)
print("(For rate r convergence, this ratio approaches r as n grows.)")
print()
for q in sorted(T_data):
    ks = sorted(T_data[q])
    diffs = [(k, T_data[q][k] - T_data[q][k-1]) for k in ks if k-1 in T_data[q]]
    print(f"  q={q}:  T_n - T_{{n-1}} for n = {[d[0] for d in diffs]}:")
    for k, d in diffs:
        print(f"      n={k}: {float(d):+.6e}")
    if len(diffs) >= 2:
        ratios = []
        for i in range(1, len(diffs)):
            d_prev = diffs[i-1][1]
            d_curr = diffs[i][1]
            r = d_curr / d_prev if d_prev != 0 else None
            ratios.append((diffs[i][0], r))
            print(f"      ratio at n={diffs[i][0]}:  "
                  f"({float(d_curr):+.4e}) / ({float(d_prev):+.4e}) "
                  f"= {float(r):+.6f}" if r is not None else "      ratio: undefined")
        last_r = float(ratios[-1][1])
        print(f"      LAST RATIO: {last_r:+.6f}    (rate-1/2 reference: 0.5)")
    print()

# --------- Diagnostic B: Aitken estimate of T(q) and envelope ---------
print("=" * 78)
print("Diagnostic B: Aitken T(q) estimate + |T_n - T(q)| * 2^n envelope")
print("=" * 78)
print()
def aitken(s0, s1, s2):
    den = s2 - 2*s1 + s0
    if den == 0: return None
    return s2 - (s2 - s1)**2 / den

for q in sorted(T_data):
    ks = sorted(T_data[q])
    if len(ks) < 3:
        # Try direct fit using only 2 T values: cannot Aitken
        print(f"  q={q}:  only {len(ks)} T values; cannot Aitken. "
              f"|T_2 - T_1| = {float(abs(T_data[q][ks[-1]] - T_data[q][ks[0]])):.4e}")
        continue
    a, b, c = ks[-3], ks[-2], ks[-1]
    Ta, Tb, Tc = T_data[q][a], T_data[q][b], T_data[q][c]
    Tinf = aitken(Ta, Tb, Tc)
    if Tinf is None:
        print(f"  q={q}: Aitken denom zero, skipping")
        continue
    print(f"  q={q}: T_inf estimate = {float(Tinf):.10f}  "
          f"[Aitken on T_{a}, T_{b}, T_{c}]")
    print(f"          |T_n - T_inf| * 2^n envelope:")
    for n in ks:
        eps = T_data[q][n] - Tinf
        env = abs(float(eps)) * (2 ** n)
        print(f"             n={n}: |T_n - T_inf|={float(abs(eps)):.4e},  "
              f"|.|*2^n = {env:.4e}")
    print()

# --------- Honest summary ---------
print("=" * 78)
print("Summary")
print("=" * 78)
print("""
For q=3 (5 T values, T(3) = 7/15): the |epsT_n|*2^n envelope is flat near
0.04 and the empirical |T_{n+1}-T_n|/|T_n-T_{n-1}| ratio approaches 0.5
in the last available step (0.4924 between n=4..5). Rate 1/2 confirmed.

For q=5, q=7 (3 T values each): only ONE empirical difference-ratio
available, computed at low n. Empirical ratios:
  q=5: ratio at n=3 = +0.063  (= 1/16, NOT 1/2)
  q=7: ratio at n=3 = -0.011  (sign flip; |ratio| ~ 1/87)
These are far from 1/2. But: with only one ratio per q, we cannot
distinguish "asymptotic rate is q-specific and not 1/2" from "rate-1/2
emerges at higher n with a long transient." Test 2/3 compute bounds
limit n at q=5, 7 to k=3, so the question is open.

For q=11, q=13 (2 T values each): no rate-difference-ratio computable.
|T_2 - T_1| values are small (4e-4 and 7e-5) suggesting fast convergence,
but ONE point cannot establish a rate.

NET: the rate-1/2 envelope at the *normalized* level is confirmed at q=3,
INCONCLUSIVE for q >= 5 with the current compute bounds. A genuine
q-universal rate finding would require pushing q=5 to k>=5 and q=7 to k>=4
to give multi-point envelope tests at the normalized level.
""")
