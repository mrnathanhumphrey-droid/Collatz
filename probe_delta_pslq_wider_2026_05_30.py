"""
probe_delta_pslq_wider_2026_05_30.py

Widen the PSLQ search basis. Constants of interest in this problem:
  - log(2), log(17), log(255), log(256), log(255/256), log(15/17)
  - polylogs Li_s(1/2), Li_s(1/256), Li_s(1/q) for q=17
  - sqrt(17), 1/127, 19/127, 15/127
  - hypergeometric pFq evaluations

Also extrapolate c_∞ via Shanks on c(2), c(3), c(4) and PSLQ Δ_∞ alone.

Tighter PSLQ tolerance (10^-13) to filter out noise relations.
"""
from __future__ import annotations
import sys, gc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

# Use the previously computed high-precision values
delta_vals = {
    1: 3.5719308420872459e-03,    # EXACT from c(1) - 19/127
    2: 3.6416215664888670e-03,
    3: 3.3990324789156260e-03,
    4: 3.3824098342225830e-03,
}

c_vals = {
    1: 0.15317823005468550,
    2: 0.15324792077908730,
    3: 0.15300533169151410,
    4: 0.15298870904682100,
}

c0 = 19 / 127

# Shanks extrapolation of c(m) → c_∞ from m=2,3,4
def shanks(a, b, c):
    den = c - 2*b + a
    return c - (c - b)**2 / den if abs(den) > 1e-30 else c

c_inf_shanks = shanks(c_vals[2], c_vals[3], c_vals[4])
delta_inf = c_inf_shanks - c0
print(f"Shanks c_∞ (from c(2),c(3),c(4)) = {c_inf_shanks:.15f}")
print(f"           Δ_∞                   = {delta_inf:+.15e}")

# Iterated Shanks if we had c(5)
# (We don't have it at high precision from this run; the metaphor probe gave c(5)≈0.152988999414)
c5 = 0.152988999414
c6_shanks = shanks(c_vals[3], c_vals[4], c5)
delta_6_shanks = c6_shanks - c0
print(f"Shanks c_∞ (from c(3),c(4),c(5)) = {c6_shanks:.15f}")
print(f"           Δ_∞                   = {delta_6_shanks:+.15e}")

# More refined: double-Shanks
c_inf_dbl = shanks(c_inf_shanks, c6_shanks, c_vals.get(5, c5))  # rough
print(f"\nNote: c_∞ ≈ 0.152989 region; Δ_∞ ≈ 3.383e-3")

# Now PSLQ on Δ_1, Δ_2, Δ_3, Δ_4 against a much wider basis
try:
    from mpmath import mp, mpf, mpc, pslq, identify, polylog, sqrt, log, pi, e, exp
    from mpmath import lerchphi, hyper, gamma, zeta
    mp.dps = 30

    d = {m: mpf(repr(delta_vals[m])) for m in (1, 2, 3, 4)}
    d_inf = mpf(repr(delta_inf))

    # rich basis
    basis = {
        "1": mpf(1),
        "1/q": mpf(1)/mpf(17),
        "1/q^2": mpf(1)/mpf(17)**2,
        "1/255": mpf(1)/mpf(255),
        "1/256": mpf(1)/mpf(256),
        "log(2)": log(mpf(2)),
        "log(17)": log(mpf(17)),
        "log(255)": log(mpf(255)),
        "log(256)": log(mpf(256)),
        "log(255/256)": log(mpf(255)/mpf(256)),
        "log(15/16)": log(mpf(15)/mpf(16)),
        "log(15/17)": log(mpf(15)/mpf(17)),
        "log(1-1/q)": log(1 - mpf(1)/mpf(17)),
        "log(1-1/256)": log(1 - mpf(1)/mpf(256)),
        "log(1-2/q)": log(1 - mpf(2)/mpf(17)),
        "sqrt(17)": sqrt(mpf(17)),
        "Li_2(1/256)": polylog(2, mpf(1)/mpf(256)),
        "Li_3(1/256)": polylog(3, mpf(1)/mpf(256)),
        "Li_2(1/q)": polylog(2, mpf(1)/mpf(17)),
        "Li_3(1/q)": polylog(3, mpf(1)/mpf(17)),
        "Li_2(1/2)": polylog(2, mpf(1)/mpf(2)),
        "Li_3(1/2)": polylog(3, mpf(1)/mpf(2)),
        "Li_2(-1/256)": polylog(2, -mpf(1)/mpf(256)),
        "Li_2(-1/q)": polylog(2, -mpf(1)/mpf(17)),
    }

    basis_names = list(basis.keys())
    basis_vals = list(basis.values())

    print(f"\n=== PSLQ with wider basis, tol=10^-15, maxcoeff=10^6 ===")
    print(f"Basis size: {len(basis_vals)}")
    for cand_name, cand_val in [("Δ_1", d[1]), ("Δ_2", d[2]), ("Δ_3", d[3]), ("Δ_4", d[4]), ("Δ_∞ (Shanks)", d_inf)]:
        rel = pslq([cand_val] + basis_vals, maxcoeff=10**6, tol=mpf(10)**(-15))
        if rel is None:
            print(f"  {cand_name}: None")
        else:
            terms = [f"{rel[0]}·{cand_name}"]
            for c, name in zip(rel[1:], basis_names):
                if c != 0:
                    terms.append(f"{c:+d}·{name}")
            print(f"  {cand_name}: {' '.join(terms)}")

    # Also: Δ_m × q^m for various m (maybe Δ_m has q^{-m} scaling)
    print(f"\n=== Test Δ_m * q^m for hidden scale ===")
    for m in (1, 2, 3, 4):
        x = d[m] * mpf(17)**m
        print(f"  Δ_{m} * 17^{m} = {x}")
        ident = identify(x, ["log(2)", "log(17)", "log(15)", "sqrt(17)", "log(255/256)"])
        print(f"    identify: {ident}")

    # Also test Δ_m * 256^m  (since 2^8 = 256 appears repeatedly)
    print(f"\n=== Test Δ_m * 256^m ===")
    for m in (1, 2, 3, 4):
        x = d[m] * mpf(256)**m
        print(f"  Δ_{m} * 256^{m} = {x}")

    # PSLQ specifically: differences Δ_{m+1} - Δ_m
    print(f"\n=== PSLQ on c(m+1) - c(m) sequence ===")
    diffs = {
        "c(2)-c(1)": d[2] - d[1],
        "c(3)-c(2)": d[3] - d[2],
        "c(4)-c(3)": d[4] - d[3],
    }
    for name, val in diffs.items():
        print(f"  {name} = {val}")
        rel = pslq([val] + basis_vals, maxcoeff=10**5, tol=mpf(10)**(-14))
        if rel is not None:
            terms = [f"{rel[0]}·({name})"]
            for c, n in zip(rel[1:], basis_names):
                if c != 0:
                    terms.append(f"{c:+d}·{n}")
            print(f"    PSLQ: {' '.join(terms)}")

except Exception as e:
    print(f"PSLQ section failed: {e}")
    import traceback; traceback.print_exc()
