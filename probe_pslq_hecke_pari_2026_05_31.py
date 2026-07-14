"""
PSLQ c_inf against PARI Hecke L-values for Q(i) at conductors above {2, 5, 13, 17}.

Strategy:
  1. Load 142 L1, L2 values from PARI dump.
  2. Build basis: {1, logs of small primes, sqrt-reciprocals, pi} + all Hecke L Re/Im
  3. PSLQ at progressively finer tolerances (10^-10 down to 10^-50).
  4. Print any nontrivial relation with bounded coefficients.

If null at 50 digits over this basis: the c_inf is NOT a finite linear combination of
Hecke L_1, L_2 values over Q(i) at conductors dividing primes {2, 5, 13, 17} with
small integer coefficients. That rules out the Q(i)-Hecke-L hypothesis.
"""
from __future__ import annotations
import sys, csv
from mpmath import mp, mpf, mpc, log, sqrt, pi, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 60

c_inf = mpf("0.15298912060588517527891674877413229926086222622334")

# Load PARI output
lvals = []
with open("C:/Collatz/hecke_pari_lvalues_2026_05_31.csv") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("label"):
            continue
        parts = line.split(";")
        if len(parts) != 6:
            continue
        try:
            label = parts[0]
            chi = parts[1]
            L1re = mpf(parts[2])
            L1im = mpf(parts[3])
            L2re = mpf(parts[4])
            L2im = mpf(parts[5])
            lvals.append((label, chi, L1re, L1im, L2re, L2im))
        except Exception as e:
            print(f"skip: {line[:80]} -> {e}")

print(f"Loaded {len(lvals)} L-value rows")

# Deduplicate near-identical values (within 1e-50)
def dedup(vals, tol=mpf(10)**(-50)):
    kept = []
    kept_names = []
    for name, v in vals:
        if abs(v) < tol:
            continue
        is_dup = False
        for kv in kept:
            if abs(v - kv) < tol or abs(v + kv) < tol:
                is_dup = True
                break
        if not is_dup:
            kept.append(v)
            kept_names.append(name)
    return kept_names, kept

# Build candidate basis
basis_raw = []
basis_raw.append(("1", mpf(1)))
basis_raw.append(("log2", log(mpf(2))))
basis_raw.append(("log3", log(mpf(3))))
basis_raw.append(("log5", log(mpf(5))))
basis_raw.append(("log13", log(mpf(13))))
basis_raw.append(("log17", log(mpf(17))))
basis_raw.append(("pi", pi))
for n in [2, 3, 5, 13, 17, 65, 85, 221]:
    basis_raw.append((f"1/sqrt({n})", 1 / sqrt(mpf(n))))

# Add Hecke L values
for label, chi, L1re, L1im, L2re, L2im in lvals:
    if abs(L1re) > mpf(10)**(-50):
        basis_raw.append((f"L1re[{label};{chi}]", L1re))
    if abs(L1im) > mpf(10)**(-50):
        basis_raw.append((f"L1im[{label};{chi}]", L1im))
    if abs(L2re) > mpf(10)**(-50):
        basis_raw.append((f"L2re[{label};{chi}]", L2re))
    if abs(L2im) > mpf(10)**(-50):
        basis_raw.append((f"L2im[{label};{chi}]", L2im))

print(f"Raw basis size: {len(basis_raw)}")

# Dedupe
names_dd, vals_dd = dedup(basis_raw, tol=mpf(10)**(-50))
print(f"Deduped basis size: {len(vals_dd)}")

# ====== PSLQ ======
print("\n=== PSLQ c_inf against full basis ===")
print(f"  c_inf = {c_inf}")
print()

for tol_exp in [10, 15, 20, 25, 30, 35, 40, 45]:
    tol = mpf(10) ** (-tol_exp)
    print(f"--- tol = 10^-{tol_exp}, maxcoeff = 10^4 ---")
    try:
        rel = pslq([c_inf] + vals_dd, tol=tol, maxcoeff=10**4)
        if rel is None:
            print("  no relation")
        else:
            terms = [f"({rel[0]:+d})*c_inf"]
            for c, name in zip(rel[1:], names_dd):
                if c != 0:
                    terms.append(f"({c:+d})*{name}")
            residual = sum(c*v for c, v in zip(rel, [c_inf]+vals_dd))
            print(f"  RELATION: {' '.join(terms)}")
            print(f"  residual = {float(residual):.2e}")
    except Exception as e:
        print(f"  pslq error: {e}")

# Try single-ratio scan
print("\n=== Single-ratio scan: c_inf / b for b in basis ===")
from fractions import Fraction
hits = []
for name, v in zip(names_dd, vals_dd):
    if abs(v) < mpf("1e-30"):
        continue
    r = c_inf / v
    rf = float(r)
    if abs(rf) > 1e10 or abs(rf) < 1e-10:
        continue
    f = Fraction(rf).limit_denominator(1000)
    approx = mpf(f.numerator)/mpf(f.denominator)
    diff = r - approx
    if abs(diff) < mpf("1e-8"):
        hits.append((name, f, float(r), float(diff)))

if hits:
    print(f"  Found {len(hits)} suspicious ratios (diff < 1e-8):")
    for name, frac, ratio, diff in hits[:20]:
        print(f"    c_inf / ({name}) = {ratio:.20f} ~ {frac}, diff = {diff:.2e}")
else:
    print("  No suspicious single-ratio hits")

# Try pair-ratio: c_inf = a*L_i + b*L_j with small a, b
print("\n=== Pair PSLQ: c_inf vs each (L_i, L_j) ===")
pair_hits = []
n = len(vals_dd)
for i in range(min(20, n)):  # only check Hecke pairs against few candidates
    if not names_dd[i].startswith("L"):
        continue
    for j in range(i+1, n):
        if not names_dd[j].startswith("L"):
            continue
        try:
            rel = pslq([c_inf, vals_dd[i], vals_dd[j], mpf(1)], tol=mpf(10)**(-40), maxcoeff=200)
            if rel is None:
                continue
            if rel[0] == 0:
                continue
            residual = rel[0]*c_inf + rel[1]*vals_dd[i] + rel[2]*vals_dd[j] + rel[3]
            if abs(residual) < mpf(10)**(-40):
                pair_hits.append((names_dd[i], names_dd[j], rel))
        except Exception:
            pass

if pair_hits:
    print(f"  Found {len(pair_hits)} pair relations")
    for ni, nj, rel in pair_hits[:10]:
        print(f"    ({rel[0]})*c_inf + ({rel[1]})*({ni}) + ({rel[2]})*({nj}) + ({rel[3]}) = 0")
else:
    print("  No pair relations under maxcoeff=200")

print("\n=== Summary ===")
print(f"Tested c_inf against {len(vals_dd)} basis elements")
print(f"Basis includes: {sum(1 for n in names_dd if n.startswith('L'))} Hecke L-values + {sum(1 for n in names_dd if not n.startswith('L'))} elementary constants")
