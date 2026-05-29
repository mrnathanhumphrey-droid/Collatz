"""
probe_qx1_multicoset_2026_05_29.py

Richer multi-coset case: q=31. ord_2 mod 31 = 5 (2^5=32=1) => <2>={1,2,4,8,16}, phi=30,
INDEX 6 => SIX cosets of <2> in (Z/31^n)*. And -1=30 not in <2> (ord_2 odd), so this is the
q=7 *family* (conjugation-paired) but with 6 cosets.

Structure: the quotient (Z/31^n)*/<2> ~= Z/6; conjugation xi->-xi is its order-2 element, so it
pairs the 6 cosets into 3 CONJUGATE pairs. Reality of P (|mu_hat(-xi)|=|mu_hat(xi)|) => the two
cosets in each pair carry EQUAL Plancherel mass. So at most 3 distinct mass levels.

QUESTION: are the 3 pair-levels distinct (a gradient), and does the <2>-pair (the powers-of-2
resonance, containing 1) carry the most? Prediction: yes to both -- 3 distinct levels, <2>-pair top.
"""
from __future__ import annotations
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
from probe_qx1_coset_plancherel_2026_05_29 import offset_distribution, subgroup_pow2

def coset_labels(q, n, Helems):
    """Assign each unit a coset id of <2> via BFS multiplication; non-units get -1."""
    N = q ** n
    lab = np.full(N, -1, dtype=np.int64)
    Helems = np.asarray(Helems, dtype=np.int64)
    nxt = 0
    for u in range(1, N):
        if u % q == 0 or lab[u] != -1:
            continue
        members = (u * Helems) % N
        lab[members] = nxt
        nxt += 1
    return lab, nxt

def main():
    q = 31
    print("=" * 100)
    print(f"q={q}: <2> index 6 -> 6 cosets, -1 not in <2> -> 3 conjugate-equal pairs")
    print("=" * 100)
    for n in (2, 3, 4):
        N = q ** n
        P = offset_distribution(q, n)
        assert abs(P.sum() - 1.0) < 1e-9
        mu = np.fft.fft(P)
        a2 = np.abs(mu) ** 2
        H, ord2 = subgroup_pow2(q, n)
        Helems = np.nonzero(H)[0]
        phi = q ** (n - 1) * (q - 1)
        idx = phi // ord2
        lab, ncos = coset_labels(q, n, Helems)
        units = (np.arange(N) % q != 0)
        mtot = a2[units].sum()
        # per-coset mass + max
        coset_mass = np.array([a2[lab == c].sum() / mtot for c in range(ncos)])
        coset_max = np.array([np.abs(mu[lab == c]).max() for c in range(ncos)])
        H_coset = lab[1 % N]   # coset containing 1 = <2> itself
        # conjugate partner of each coset
        conj_of = np.array([lab[(N - np.nonzero(lab == c)[0][0]) % N] for c in range(ncos)])
        print(f"\n n={n}  N={N}  ord2={ord2}  phi={phi}  index={idx}  #cosets={ncos}  <2>=coset {H_coset}")
        print(f"  {'coset':>5} {'massfrac':>10} {'max|mu|':>10} {'conj_partner':>12} {'is <2>?':>8}")
        order = np.argsort(-coset_mass)
        for c in order:
            print(f"  {c:>5} {coset_mass[c]:>10.6f} {coset_max[c]:>10.5f} {conj_of[c]:>12} "
                  f"{'YES' if c == H_coset else '':>8}")
        # group into conjugate pairs, list distinct levels
        seen = set(); pairs = []
        for c in range(ncos):
            if c in seen: continue
            d = conj_of[c]; seen.add(c); seen.add(d)
            pairs.append((c, d, coset_mass[c], coset_mass[d]))
        print("  conjugate pairs (mass_c, mass_partner) -- should match within pair:")
        for c, d, mc, md in pairs:
            tag = " <-- <2>-pair" if (c == H_coset or d == H_coset) else ""
            print(f"    cosets ({c},{d}): {mc:.6f} , {md:.6f}   |diff|={abs(mc-md):.2e}{tag}")
        levels = sorted({round(mc, 6) for c, d, mc, md in pairs}, reverse=True)
        print(f"  distinct pair-levels (rounded 6dp): {levels}   ({len(levels)} distinct)")

    print("\nReading: 3 conjugate pairs, equal mass within each pair (reality of P).")
    print("  If 3 distinct levels with <2>-pair highest => graded multi-coset resonance, as predicted.")

if __name__ == "__main__":
    main()
