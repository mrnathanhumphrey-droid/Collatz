"""
result_q_sweep_test_1_followup.py — quick exact-rational followup that
verifies the structural finding hidden in Test 1's data:

    S_n^{(q)} grows by factor q/3 per refinement step,
    with q=3 the critical case (q/3 = 1) where S_n converges.

This means:
   T_n^{(q)} := S_n^{(q)} * (3/q)^n   should be approximately constant in n.

Computes T_n^{(q)} as exact rationals for the data already in
result_q_sweep_test_1_envelope.csv and prints the table.
"""
from __future__ import annotations
import csv
import sys
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

ENV = r"C:\Collatz\result_q_sweep_test_1_envelope.csv"

S_data = {}  # q -> {n: Fraction}
with open(ENV) as f:
    for row in csv.DictReader(f):
        q = int(row["q"])
        n = int(row["n"])
        Sn = Fraction(int(row["S_n_num"]), int(row["S_n_den"]))
        S_data.setdefault(q, {})[n] = Sn

print("=" * 78)
print("Followup: structural normalization S_n^{(q)} * (3/q)^n")
print("=" * 78)
print()
print("Hypothesis (q-universal growth law):")
print("    S_n^{(q)} ~ T(q) * (q/3)^n  with q-dependent constant T(q)")
print()
print("Verification: compute T_n^{(q)} := S_n^{(q)} * (3/q)^n exactly.")
print()

print(f"  {'q':>3}  {'n':>2}  {'S_n':>15}  {'(q/3)^n':>15}  "
      f"{'T_n=S_n*(3/q)^n':>20}  {'ratio S_n/S_{n-1}':>20}")
prev = {}
for q in sorted(S_data):
    for n in sorted(S_data[q]):
        Sn = S_data[q][n]
        scale = Fraction(3, q) ** n
        Tn = Sn * scale
        ratio = (Sn / prev.get(q)) if q in prev else None
        prev[q] = Sn
        ratio_str = (f"{float(ratio):.6f}" if ratio is not None else "(none)")
        print(f"  {q:>3}  {n:>2}  {float(Sn):>15.6f}  {float((Fraction(q,3))**n):>15.6f}  "
              f"{float(Tn):>20.6f}  {ratio_str:>20}")
    print()

print()
print("=" * 78)
print("Q-universal growth ratio test")
print("=" * 78)
print()
print(f"  {'q':>3}  {'q/3':>8}  {'last empirical ratio':>22}  {'rel err':>10}")
for q in sorted(S_data):
    ks = sorted(S_data[q])
    if len(ks) < 2:
        print(f"  {q:>3}  {float(Fraction(q,3)):>8.4f}  {'<2 points':>22}")
        continue
    last_ratio = S_data[q][ks[-1]] / S_data[q][ks[-2]]
    qover3 = Fraction(q, 3)
    rel = abs(float(last_ratio) - float(qover3)) / float(qover3)
    print(f"  {q:>3}  {float(qover3):>8.4f}  {float(last_ratio):>22.6f}  {rel:>10.2e}")

# T(q) limit estimates: average of T_n for n>=2 (skip transient at n=1)
print()
print("=" * 78)
print("T(q) prefactor estimates (average over n>=2)")
print("=" * 78)
print()
print(f"  {'q':>3}  {'T(q) ~ avg T_n for n>=2':>30}")
for q in sorted(S_data):
    Ts = []
    for n in sorted(S_data[q]):
        if n < 2:
            continue
        Ts.append(float(S_data[q][n] * Fraction(3, q) ** n))
    if not Ts:
        print(f"  {q:>3}  (no n>=2 data)")
        continue
    avg = sum(Ts) / len(Ts)
    print(f"  {q:>3}  {avg:>30.6f}    (n's used: {[n for n in sorted(S_data[q]) if n>=2]})")
