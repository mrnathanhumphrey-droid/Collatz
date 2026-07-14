"""
PSLQ c_inf against 1188 Q(i) grossencharacter L-values from PARI gchar,
covering infinity types (1,0), (2,0), (1,1), (0,1), (0,2), (-1,0), (1,-1),
(2,1), (1,2), (3,0), (0,3) at 14 conductors.

This is the FULL Hecke L test the agent recommended (modulo L'(0) which is
folded into the functional eq from L(1) of dual).
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 60

c_inf = mpf("0.15298912060588517527891674877413229926086222622334")

lvals = []
with open("C:/Collatz/hecke_gchar_lvalues_2026_05_31.csv") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("label"):
            continue
        parts = line.split(";")
        if len(parts) != 7:
            continue
        try:
            label, inftype, chi, L1re, L1im, L2re, L2im = parts
            lvals.append((label, inftype, chi, mpf(L1re), mpf(L1im), mpf(L2re), mpf(L2im)))
        except Exception:
            pass

print(f"Loaded {len(lvals)} gchar L-value rows")

# Build basis
basis_raw = []
basis_raw.append(("1", mpf(1)))
for name, v in [("log2", log(mpf(2))), ("log3", log(mpf(3))), ("log5", log(mpf(5))),
                ("log13", log(mpf(13))), ("log17", log(mpf(17))), ("pi", pi)]:
    basis_raw.append((name, v))
for n in [2, 3, 5, 13, 17, 65, 85, 221]:
    basis_raw.append((f"1/sqrt({n})", 1 / sqrt(mpf(n))))

for label, inftype, chi, L1re, L1im, L2re, L2im in lvals:
    tag = f"{label};{inftype};{chi}"
    if abs(L1re) > mpf(10)**(-50):
        basis_raw.append((f"L1re[{tag}]", L1re))
    if abs(L1im) > mpf(10)**(-50):
        basis_raw.append((f"L1im[{tag}]", L1im))
    if abs(L2re) > mpf(10)**(-50):
        basis_raw.append((f"L2re[{tag}]", L2re))
    if abs(L2im) > mpf(10)**(-50):
        basis_raw.append((f"L2im[{tag}]", L2im))

print(f"Raw basis size: {len(basis_raw)}")

# Dedupe with stricter tol
def dedup(items, tol=mpf(10)**(-45)):
    kept_names, kept_vals = [], []
    for name, v in items:
        if abs(v) < tol:
            continue
        dup = False
        for kv in kept_vals:
            if abs(v - kv) < tol or abs(v + kv) < tol:
                dup = True
                break
        if not dup:
            kept_names.append(name)
            kept_vals.append(v)
    return kept_names, kept_vals

names, vals = dedup(basis_raw)
print(f"Deduped basis size: {len(vals)}")
n_lvals_in_basis = sum(1 for n in names if n.startswith("L"))
print(f"  Hecke L-values: {n_lvals_in_basis}")
print(f"  Elementary: {len(names) - n_lvals_in_basis}")

# Full-basis PSLQ at multiple tols
print()
print("=== Full-basis PSLQ ===")
print(f"c_inf = {c_inf}")
print()

found_nontrivial = False
for tol_exp in [10, 15, 20, 25, 30, 35, 40, 45]:
    tol = mpf(10) ** (-tol_exp)
    rel = pslq([c_inf] + vals, tol=tol, maxcoeff=10**4)
    if rel is None:
        print(f"  tol=10^-{tol_exp}: no relation")
        continue
    # Print relation
    if rel[0] != 0:
        found_nontrivial = True
        terms = [f"({rel[0]:+d})*c_inf"]
        for c, n in zip(rel[1:], names):
            if c != 0:
                terms.append(f"({c:+d})*{n}")
        residual = sum(c*v for c, v in zip(rel, [c_inf]+vals))
        print(f"  tol=10^-{tol_exp}: NONTRIVIAL relation involving c_inf!")
        print(f"    {' '.join(terms)}")
        print(f"    residual = {float(residual):.2e}")
    else:
        nonzero = [(c, n) for c, n in zip(rel[1:], names) if c != 0]
        if len(nonzero) <= 6:
            print(f"  tol=10^-{tol_exp}: trivial rel (c_inf coef 0, {len(nonzero)} terms)")
        else:
            print(f"  tol=10^-{tol_exp}: trivial rel (c_inf coef 0, {len(nonzero)} terms — L-value dependency)")

# Subset focused on (1,0) and (1,1) — agent's recommended types
print()
print("=== Focused PSLQ: c_inf vs L-values at infty (1,0) and (1,1) only ===")
focused_names = []
focused_vals = []
focused_names.append("1"); focused_vals.append(mpf(1))
focused_names.append("log2"); focused_vals.append(log(mpf(2)))
focused_names.append("log17"); focused_vals.append(log(mpf(17)))
focused_names.append("pi"); focused_vals.append(pi)
focused_names.append("1/sqrt(17)"); focused_vals.append(1/sqrt(mpf(17)))
focused_names.append("1/sqrt(5)"); focused_vals.append(1/sqrt(mpf(5)))
focused_names.append("1/sqrt(85)"); focused_vals.append(1/sqrt(mpf(85)))

for label, inftype, chi, L1re, L1im, L2re, L2im in lvals:
    if inftype not in ["(1,0)", "(1,1)"]:
        continue
    tag = f"{label};{inftype};{chi[:20]}"
    for n, v in [(f"L1re_{tag}", L1re), (f"L1im_{tag}", L1im),
                 (f"L2re_{tag}", L2re), (f"L2im_{tag}", L2im)]:
        if abs(v) > mpf(10)**(-45):
            focused_names.append(n)
            focused_vals.append(v)

focused_names, focused_vals = dedup(list(zip(focused_names, focused_vals)))
print(f"  Focused basis size: {len(focused_vals)}")
for tol_exp in [15, 25, 35, 45]:
    tol = mpf(10) ** (-tol_exp)
    rel = pslq([c_inf] + focused_vals, tol=tol, maxcoeff=10**4)
    if rel is None:
        print(f"  tol=10^-{tol_exp}: no relation")
    elif rel[0] != 0:
        print(f"  tol=10^-{tol_exp}: NONTRIVIAL c_inf relation, see below")
        terms = [f"({rel[0]:+d})*c_inf"] + [f"({c:+d})*{n}" for c, n in zip(rel[1:], focused_names) if c != 0]
        print(f"    {' '.join(terms)}")
    else:
        nonzero = sum(1 for c in rel[1:] if c != 0)
        print(f"  tol=10^-{tol_exp}: trivial (c_inf coef 0, {nonzero} terms)")

# Single-ratio scan against weight-(1,0) L-values
print()
print("=== Single-ratio scan: c_inf / L for L at infty (1,0) or (1,1) ===")
from fractions import Fraction
hits = []
for label, inftype, chi, L1re, L1im, L2re, L2im in lvals:
    if inftype not in ["(1,0)", "(1,1)", "(0,1)"]:
        continue
    tag = f"{label};{inftype};{chi[:20]}"
    for n, v in [(f"L1re_{tag}", L1re), (f"L1im_{tag}", L1im),
                 (f"L2re_{tag}", L2re), (f"L2im_{tag}", L2im)]:
        if abs(v) < mpf("1e-30"):
            continue
        r = c_inf / v
        rf = float(r)
        if abs(rf) > 1e8 or abs(rf) < 1e-8:
            continue
        f = Fraction(rf).limit_denominator(500)
        approx = mpf(f.numerator)/mpf(f.denominator)
        diff = r - approx
        if abs(diff) < mpf("1e-10"):
            hits.append((n, f, float(r), float(diff)))

if hits:
    print(f"  {len(hits)} suspicious ratios (diff < 1e-10, denom <= 500):")
    for name, frac, ratio, diff in sorted(hits, key=lambda x: abs(x[3]))[:20]:
        print(f"    {ratio:.20f} = {frac}, name={name}, diff={diff:.2e}")
else:
    print("  No suspicious single-ratio hits at 1e-10")

print()
print("=== Summary ===")
print(f"PSLQ c_inf vs {len(vals)} basis elements ({n_lvals_in_basis} Hecke L), 50-digit precision")
print(f"Nontrivial c_inf relation found: {found_nontrivial}")
