"""
PROBE QDIFF-2 (Wilson spec, Rung 3) -- level-1 difference-Galois group of A_1 = D_1(I+N).

Pinned (Wilson pen): diagonal d = PURE z-power (z^4), NOT 3*z^k. The 3^i in a_i=3^i R_e is
substitution-NORMALIZATION (tied to d_i=2*3^{i-1}, the z->z^3 orbit), hence a sigma-coboundary,
hence gauged away -- importing it as a free constant would be the coincident-3 trap.

Difference field (K, sigma):  K = Q(z),  sigma: z -> z^3  (Mahler substitution).
Constants C = {f: sigma(f)=f} = Q. Entries of A_1 live in Q(z) (no cyclotomic extension needed).

A_1 = D_1*(I+N),  D_1 = diag(z^4, 1),  N = [[1,1],[0,0]] (N^2=0),  I+N = [[1,1],[0,1]]
   => A_1(z) = [[z^4, z^4],[0, 1]]   (every entry a named element of Q(z))
System sigma(Y) = A_1 Y.

The group collapses to TWO decidable sigma-coboundary tests (this is the Lean-ready core):
  (C1) MULTIPLICATIVE, diagonal d=z^4:  does sigma(y)=d*y have a solution y in Q(z)*?
        YES (y=z^2) => d is a coboundary => G_m COLLAPSES.
  (C2) ADDITIVE, reduced off-diagonal e=z^{-2}: does sigma(psi)-psi=e have a solution psi in Q(z)?
        NO (forced infinite principal part) => e is a non-coboundary => G_a SURVIVES.
  => G_1 = G_a (unipotent only), a PROPER subgroup of Aff_1 = G_a |x G_m.

The G_m is a PRO-object phenomenon: coboundary at every finite level, non-gaugeable only in the
inverse limit (the accumulated prod 3^i = MAHLER's doubly-exp denominators). Fifth confirmation of
the structure(graded/G_m)-vs-value(unipotent) split. Everything exact; goal is a proof not a number.
"""
from fractions import Fraction as F

# Laurent polynomials over Q as {exponent: coeff}; sigma multiplies exponents by 3.
def sigma(f):        return {3 * k: c for k, c in f.items()}
def sub(f, g):
    out = dict(f)
    for k, c in g.items():
        out[k] = out.get(k, 0) - c
        if out[k] == 0: del out[k]
    return out
def mul_mono(f, k0, c0):  # multiply by c0*z^k0
    return {k + k0: c * c0 for k, c in f.items()}


def mult_coboundary_monomial(m):
    """sigma(y)=z^m * y, y=z^a: 3a=m+a => a=m/2. Rational-fn solution in Q(z) iff a is an integer."""
    if m % 2 == 0:
        a = m // 2
        # verify exactly: sigma(z^a) == z^m * z^a
        lhs = sigma({a: F(1)}); rhs = {m + a: F(1)}
        assert lhs == rhs, "solution check failed"
        return True, a
    return False, None


def additive_coboundary_cascade(p, J=8):
    """sigma(psi)-psi = z^{-p}. Match Laurent coeffs at z=0. Coeff of z^n in sigma(psi)-psi is
    c_{n/3}[3|n] - c_n. Solve forward => forces c_{-p*3^j} = -1 for all j>=0 (infinite principal
    part) => psi not rational. Return the forced coefficients (all should be -1 = nonzero)."""
    forced = {}
    # n = -p : (no c_{-p/3} unless 3|p) - c_{-p} = 1  => c_{-p} = -1
    forced[-p] = F(-1)
    n = -p
    for j in range(1, J + 1):
        n3 = n * 3                      # next forced negative exponent (= -p*3^j)
        # coeff of z^{n3}:  c_{n3/3}=c_n  minus c_{n3}  = 0 (RHS has nothing there)
        forced[n3] = forced[n]          # c_{n3} = c_n = -1
        n = n3
    return forced


def main():
    print("# PROBE QDIFF-2 -- level-1 difference-Galois group of A_1 = D_1(I+N)  (EXACT)\n")
    print("## QD2-A  the object, pinned")
    print("   K = Q(z),  sigma: z -> z^3,  constants C = Q.   A_1 = D_1(I+N):")
    print("   D_1 = diag(z^4, 1),  I+N = [[1,1],[0,1]]   =>   A_1(z) = [[z^4, z^4],[0, 1]]")
    print("   (diagonal = PURE z-power; the 3^i is substitution-normalization, gauged away -- pinned, not chosen)\n")

    # ---- (C1) multiplicative coboundary: the G_m factor ----
    print("## QD2-C.1  MULTIPLICATIVE coboundary test on diagonal d = z^4  (the G_m factor)")
    ok, a = mult_coboundary_monomial(4)
    print(f"   solve sigma(y) = z^4 * y  in Q(z):  y = z^{a}  (since 3a = 4 + a => a = 2)")
    print(f"   check sigma(z^2) = z^6 = z^4 * z^2 :  EXACT")
    print(f"   => d = z^4 IS a sigma-coboundary (solution lies in the base field Q(z))")
    print(f"   => the diagonal solution is FIXED by every K-automorphism => G_m COLLAPSES (trivial).\n")

    # ---- reduce the off-diagonal (gauge by y_1=z^2, y_2=1) ----
    print("## QD2-C.2  reduce the off-diagonal, then ADDITIVE coboundary test  (the G_a factor)")
    print("   Y = [[z^2, w],[0, 1]];  sigma(w) = z^4 w + z^4.  Set w = z^2 * psi:")
    print("   z^6 psi(z^3) = z^6 psi(z) + z^4  =>  sigma(psi) - psi = z^{-2}   (reduced off-diagonal e = z^{-2})")
    forced = additive_coboundary_cascade(2, J=7)
    exps = sorted(forced.keys(), reverse=True)
    print(f"   Laurent-coeff cascade forces  c_(-2*3^j) = -1  for all j>=0:")
    print("     " + ",  ".join(f"c[{k}]={forced[k]}" for k in exps[:7]) + ",  ...")
    allneg1 = all(v == F(-1) for v in forced.values())
    print(f"   all forced coeffs = -1 (nonzero), exponents -2,-6,-18,-54,... never terminate:  {allneg1}")
    print("   => psi has an INFINITE principal part at z=0 => psi is NOT a rational function.")
    print("   (clean proof: a rational psi with a pole of order m at 0 gives sigma(psi)-psi pole order 3m,")
    print("    which cannot equal 2 (3m=2 impossible); m=0 gives LHS regular != z^{-2}. No rational solution.)")
    print("   => e = z^{-2} is a sigma-NON-coboundary => the unipotent G_a SURVIVES (genuine Mahler-log psi).\n")

    # ---- verdict ----
    print("## QD2-C  VERDICT -- the level-1 difference-Galois group")
    print("   diagonal coboundary (G_m collapses) + off-diagonal non-coboundary (G_a genuine)")
    print("   =>  G_1 = G_a  (one-dimensional unipotent), a PROPER subgroup of Aff_1 = G_a |x G_m.")
    print("   As a matrix group: G_1 = { [[1, u],[0, 1]] : u in G_a }  (NOT the full [[t,u],[0,1]]).\n")

    print("## pro-object emergence -- where the G_m actually lives")
    print("   The diagonal is a coboundary at EVERY finite level r (gaugeable z-power). But the tower of")
    print("   normalizations prod_i 3^i does NOT converge to a coboundary: its non-gaugeable growth IS")
    print("   MAHLER's doubly-exponential denominator rate. So the reductive G_m is a PURE PRO-OBJECT")
    print("   phenomenon -- invisible at every finite stage, emerging only in the inverse limit. Sharper")
    print("   than 'Aff_1 at every level': the reductive part is purely pro; the unipotent is present finitely.\n")

    print("## QD2-D  Lean-ready core (build-by-hand; mathlib has field-Galois but NO Picard-Vessiot/difference-Galois)")
    print("   The whole group computation is TWO decidable predicates about Q(z):")
    print("     P1 (coboundary):    exists y in Q(z)* : sigma(y) = z^4 * y            -- TRUE  (y=z^2)")
    print("     P2 (non-coboundary): not exists psi in Q(z) : sigma(psi)-psi = z^{-2} -- TRUE  (pole-order proof)")
    print("   Then G_1 = G_a is: (P1 => no G_m factor) and (P2 => a G_a factor). No general PV theory needed;")
    print("   the two predicates are finite/decidable statements about rational functions -- the Rung-3 Lean target.")
    print("\n   [structure at level 1, not value; 7/15 excluded regardless, floor 0.473177.]")


if __name__ == "__main__":
    main()
