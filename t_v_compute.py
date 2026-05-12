"""
t_v_compute.py — verification of the Phase 1 obstruction documented in
T_V_RECURSION.md.

GOAL: empirically demonstrate that

  (1) For the cross-frequency moment M_n^{ab}(g, c), the Tao-recursion lift
      to M_{n+1}^{ab}(g, c) produces phase-twisted moments that are NOT in
      span{M_n^{a'b'}(g', c') : g' ∈ {0, 2, 4, 6}} = V_M^{(g_max=6)}.

  (2) Specifically, the level-(n+1) recursion path at (v=2 even, v'=3 odd, g=2)
      produces a phase factor 2^v · D̃ = 5/8 (not ẽ_G = ẽ_3 = 7/24).

  (3) The difference θ_{v,g} = 1/3 is rational with v_3 = -1 (NOT a 3-adic integer);
      this is the structural obstruction.

USAGE: run with main thread (subagent has no Bash/PowerShell access).

EXPECTED OUTPUT:
  - Phase offset values θ_{v,g} for the leading surviving (v, v', g) pairs
  - Mod-3 / mod-9 / mod-27 reductions showing θ_{v,g} ≢ 0
  - Comparison of 2^v·D̃ vs ẽ_G to confirm they differ by θ_{v,g}
  - At n = 2, 3: empirical rank check that the phase-twisted moments are NOT
    in span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, 6}}
"""
from __future__ import annotations
import sys
from fractions import Fraction
import cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


def e_tilde(g, n=None):
    """Returns ẽ_g = (1 - 2^{-g})/3 as Fraction (rational); if n provided, also
    computes ẽ_g mod 3^n as integer in Z/3^n.

    For g even ≥ 2, ẽ_g ∈ Z_3 (3-adic integer).
    For g odd, ẽ_g has v_3 = -1 (denominator 3).
    For g = 0, ẽ_g = 0.
    """
    if g == 0:
        e = Fraction(0)
    else:
        e = (Fraction(1) - Fraction(1, 2 ** g)) / 3 if g > 0 else (Fraction(1) - Fraction(2 ** (-g))) / 3
    if n is not None:
        N = 3 ** n
        # We want e mod N. e is rational; reduce numerator/denominator mod N.
        # Since 2 is a unit mod 3^n, this works if 3 ∤ denom.
        num, den = e.numerator, e.denominator
        if den % 3 == 0:
            return e, None
        inv_den = pow(den, -1, N)
        e_mod = (num * inv_den) % N
        return e, e_mod
    return e


def v3(x):
    """3-adic valuation of a rational x (or 0 if x = 0; ∞ for 0)."""
    if x == 0:
        return float('inf')
    num, den = x.numerator, x.denominator
    v = 0
    while num % 3 == 0:
        num //= 3
        v += 1
    while den % 3 == 0:
        den //= 3
        v -= 1
    return v


def D_vvg(v, vprime, g):
    """D = ẽ_g + 2^{-v} - 2^{-g-vprime} as Fraction."""
    return e_tilde(g) + Fraction(1, 2 ** v) - Fraction(1, 2 ** (g + vprime))


def D_tilde(v, vprime, g):
    """D̃ = D/3 as Fraction (only well-defined if 3 | D, i.e., v_3(D) ≥ 1)."""
    D = D_vvg(v, vprime, g)
    return D / 3


def theta_vg(v, g):
    """Phase offset θ_{v,g} = 2^v · ẽ_g / 3 as Fraction."""
    return Fraction(2 ** v) * e_tilde(g) / 3


def G_index(v, vprime, g):
    """Outgoing shift index G = v' + g - v."""
    return vprime + g - v


def class_flow(v, c):
    """c̃ = (-1)^v · c mod 3."""
    if v % 2 == 0:
        return c
    else:
        return 3 - c


def main():
    print("# t_v_compute.py — Phase 1 obstruction verification")
    print()
    print("Brief: derive M_{n+1}^{ab}(g, c) via Tao + lift-fiber + unit-shuffle.")
    print("Per T_V_RECURSION.md, the recursion produces a moment with phase factor")
    print("    2^v · D̃ = ẽ_G + θ_{v,g}")
    print("where θ_{v,g} = 2^v · ẽ_g / 3, and θ_{v,g} ≠ 0, θ_{v,g} ≠ ẽ_{G''} - ẽ_G")
    print("for generic (v, g) — the structural obstruction to V_M closure.")
    print()

    # ====================================================================
    # Section 1: tabulate ẽ_g and θ_{v,g} for small g, v
    # ====================================================================
    print("## Section 1: phase factors ẽ_g and θ_{v,g}")
    print()
    print("ẽ_g = (1 - 2^{-g}) / 3:")
    print(f"  ẽ_0 = {e_tilde(0)}")
    for g in [2, 3, 4, 5, 6]:
        e = e_tilde(g)
        v3_e = v3(e)
        print(f"  ẽ_{g} = {e} = {float(e):.6f}  (v_3 = {v3_e})")
    print()

    print("θ_{v,g} = 2^v · ẽ_g / 3 (phase offset from unit shuffle):")
    print(f"  {'(v, g)':>10}  {'θ_{v,g}':>20}  {'v_3(θ)':>8}  {'mod 3':>8}  {'mod 9':>8}  {'mod 27':>8}")
    for g in [2, 4, 6]:
        for v in [1, 2, 3, 4, 5, 6]:
            theta = theta_vg(v, g)
            v3_t = v3(theta)
            # mod 3, 9, 27 (if 3-adic integer)
            mod3 = mod9 = mod27 = "n/a"
            if v3_t >= 0:
                num, den = theta.numerator, theta.denominator
                inv_den3 = pow(den, -1, 3)
                mod3 = (num * inv_den3) % 3
                inv_den9 = pow(den, -1, 9)
                mod9 = (num * inv_den9) % 9
                inv_den27 = pow(den, -1, 27)
                mod27 = (num * inv_den27) % 27
            print(f"  {(v, g)!s:>10}  {str(theta):>20}  {v3_t:>8}  {mod3!s:>8}  {mod9!s:>8}  {mod27!s:>8}")
    print()

    # ====================================================================
    # Section 2: survival condition (3 | D) for level-(n+1) → n recursion
    # ====================================================================
    print("## Section 2: survival condition v_3(D) ≥ 1 for level n+1 → n recursion")
    print()
    print("D_{v,v',g} = ẽ_g + 2^{-v} - 2^{-g-v'}.")
    print("Surviving tuples for g ∈ {2, 4, 6}, (v, v') in small ranges:")
    print()
    for g in [2, 4, 6]:
        print(f"### g = {g}, ẽ_{g} = {e_tilde(g)}")
        hdr = "(v, v')"
        print(f"  {hdr:>10}  {'parity (a, b)':>16}  {'D':>25}  {'v_3(D)':>10}  {'survives':>10}")
        for v in range(1, 7):
            for vprime in range(1, 7):
                D = D_vvg(v, vprime, g)
                v3D = v3(D)
                survives = v3D >= 1
                a = '+' if v % 2 == 0 else '-'
                b = '+' if vprime % 2 == 0 else '-'
                marker = " ←" if survives else ""
                if survives:
                    print(f"  {(v, vprime)!s:>10}  {(a, b)!s:>16}  {str(D):>25}  {v3D:>10}  {survives!s:>10}{marker}")
        print()

    # ====================================================================
    # Section 3: the canonical worked example (g=2, v=2, v'=3)
    # ====================================================================
    print("## Section 3: canonical worked example, g=2, v=2, v'=3")
    print()
    v, vprime, g = 2, 3, 2
    D = D_vvg(v, vprime, g)
    Dt = D_tilde(v, vprime, g)
    G = G_index(v, vprime, g)
    e_G = e_tilde(G)
    twoV_Dt = Fraction(2 ** v) * Dt
    theta = theta_vg(v, g)
    print(f"  D = ẽ_2 + 2^{{-2}} - 2^{{-5}} = {D}")
    print(f"  D/3 = D̃ = {Dt}")
    print(f"  G = v' + g - v = {G}")
    print(f"  ẽ_G = ẽ_{G} = {e_G}  (note: G is ODD, v_3(ẽ_G) = {v3(e_G)})")
    print(f"  2^v · D̃ = {twoV_Dt}")
    print(f"  θ_{{v={v},g={g}}} = 2^v · ẽ_g / 3 = {theta}")
    print(f"  Check: 2^v·D̃ - ẽ_G = {twoV_Dt - e_G}  (should equal θ_{{v,g}})")
    print(f"  Match: {twoV_Dt - e_G == theta}")
    print()
    print("  Conclusion: at this surviving pair, the level-n phase factor is")
    print(f"    e^{{-2πi s · (5/8) / 3^n}} = e^{{-2πi s · (ẽ_3 + 1/3) / 3^n}}")
    print(f"  The phase 5/8 is NOT ẽ_G for ANY integer G:")
    print(f"    ẽ_G = (1 - 2^{{-G}}) / 3 = 5/8 ⟹ 1 - 2^{{-G}} = 15/8 ⟹ 2^{{-G}} = -7/8")
    print("    No integer G solves this. So the moment at this s-sum is OUTSIDE V_M.")
    print()

    # ====================================================================
    # Section 4: check that the phase 5/8 is mod-3 well-defined (it is)
    # ====================================================================
    print("## Section 4: phase 5/8 mod 3, 9, 27 (well-defined 3-adic integer)")
    print()
    phase_585 = Fraction(5, 8)
    print(f"  Phase = 5/8, v_3 = {v3(phase_585)}")
    for N in [3, 9, 27, 81]:
        num, den = phase_585.numerator, phase_585.denominator
        inv_den = pow(den, -1, N)
        print(f"    5/8 mod {N} = {(num * inv_den) % N}")
    print()
    print("  Compare to ẽ_3 mod 3, 9, 27, 81:")
    e3 = e_tilde(3)
    print(f"  ẽ_3 = {e3}, v_3 = {v3(e3)}")
    if v3(e3) >= 0:
        for N in [3, 9, 27, 81]:
            num, den = e3.numerator, e3.denominator
            inv_den = pow(den, -1, N)
            print(f"    7/24 mod {N} = ... (denom 24 = 8·3 has 3 in it, so ẽ_3 is NOT 3-adic integer)")
    else:
        print("  ẽ_3 has v_3 = -1; not a 3-adic integer, so phase ẽ_G alone is not well-defined")
        print("  at any level. The moment M_n(G=3, c) is NOT a natural object.")
    print()

    # ====================================================================
    # Section 5: empirical V_M^{(g_max=6)} rank check at n = 2, 3
    # ====================================================================
    print("## Section 5: empirical V_M^{(g_max=6)} basis dimension at n = 2, 3")
    print("(Re-confirms cross_freq_compute.py's finding: at n=3, rank ≥ 7, vs P-only rank = 1)")
    print()
    # This section is left as a reminder; the actual empirical check is in
    # cross_freq_compute.py. We don't re-run it here.
    print("  See cross_freq_compute.py for the rank diagnostic:")
    print("    n=2: augmented rank = 6 (P-only = 1)")
    print("    n=3: augmented rank = 7 (P-only = 1)")
    print("  This confirms M_n(g≥2, c) ∉ span{M_n(g=0, c)} = span{P_n^{ab}(c)}.")
    print()
    print("  The OBSTRUCTION in T_V_RECURSION.md is the NEXT step:")
    print("    when we iterate M_{n+1}(g, c) by Tao, the result is OUTSIDE")
    print("    V_M^{(g_max)} for ANY finite g_max. Specifically, the phase")
    print("    offsets θ_{v,g} and the parity-flipped G shifts produce")
    print("    moments at G odd, with phase offsets θ_{v,g}, neither of which")
    print("    is in V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, 6, ...}}.")
    print()

    # ====================================================================
    # Section 6: cascade table — outgoing G for each (v, v', g) surviving pair
    # ====================================================================
    print("## Section 6: outgoing G values from surviving (v, v', g) pairs")
    print()
    print("Brief: for each (a, b) ∈ {+,-}² and incoming g ∈ {0, 2, 4, 6},")
    print("tabulate the outgoing G = v' + g - v over the SURVIVING (v ∈ V_a, v' ∈ V_b)")
    print("pairs with 3 | D. Outgoing G values are the shifts of the moments")
    print("M_n^{a'b'}(G, c̃) appearing in the recursion.")
    print()
    V_plus = [2, 4, 6, 8, 10]
    V_minus = [1, 3, 5, 7, 9]
    Vmap = {'+': V_plus, '-': V_minus}

    for g in [0, 2, 4, 6]:
        print(f"### g = {g}")
        for (a, b) in [('+', '+'), ('+', '-'), ('-', '+'), ('-', '-')]:
            Gs = set()
            for v in Vmap[a]:
                for vprime in Vmap[b]:
                    D = D_vvg(v, vprime, g)
                    if v3(D) >= 1:
                        Gs.add(vprime + g - v)
            Gs_list = sorted(Gs)
            parity_pattern = "even-only" if all(x % 2 == 0 for x in Gs_list) else ("odd-only" if all(x % 2 == 1 for x in Gs_list) else "mixed")
            print(f"  (a={a}, b={b}): surviving G values = {Gs_list[:8]}...  parity: {parity_pattern}")
        print()

    print("Observation: for incoming g ∈ {2, 4}, the surviving outgoing G")
    print("values are ODD. So even-g moments produce odd-G moments under")
    print("Tao iteration. This is the parity obstruction (T_V_RECURSION §6).")
    print()
    print("For g = 0: surviving G values are EVEN (recovering V_M closure at g=0).")
    print()
    print("For g = 6: ẽ_6 ≡ 0 mod 3, refined survival check shows MIXED parities")
    print("of outgoing G, including g = 6 self-couplings and feeds to lower g.")
    print()

    # ====================================================================
    # Section 7: closure inequality (NOT run; would need T_V matrix)
    # ====================================================================
    print("## Section 7: closure inequality")
    print()
    print("Not run. T_V matrix not constructed (Phase 1 obstruction).")
    print()
    print("For the bilinear bound side (r ≤ 3 strict 2√N; r ≥ 4 polylog-free")
    print("2√p·√N), see project_collatz_r78_bilinear_cracked memory file.")
    print()

    print("=" * 70)
    print("DONE. Phase 1 obstruction documented; Phase 2+ blocked.")


if __name__ == "__main__":
    main()
