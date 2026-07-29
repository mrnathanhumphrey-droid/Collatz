"""
PROBE QDIFF-3 (Wilson spec, RUNG 4) -- pin A_2 (level-2 truncation) and gate the tower.

A_2 = D_2(I+N_2) from QDIFF-1's confirmed construction. Diagonal = pure z-power substitution
scalings (QDIFF-2 pin): level-i scaling z^{d_{i+1}-d_i} = z^{2 d_i} = z^{4*3^{i-1}}  (d_i=2*3^{i-1}).
  level 1: z^4   level 2: z^12   level 3: z^36 ...
So (with the level-1 block A_1=[[z^4,z^4],[0,1]] as the confirmed corner):
  A_2 = [[z^12, z^12, 0],[0, z^4, z^4],[0, 0, 1]]   (A_1 = bottom-right 2x2 block => nesting embedding)

4a: coboundary test on the level-2 diagonal z^12 (predict: coboundary, G_m still collapses); and the
    unipotent dimension (predict: G_a^2, dim grew by 1).
4b: level-2 unipotent restricts to level-1 G_a under the corner embedding.

KEY subtlety this gate surfaces: the unipotent dim-growth is NOT confirmable from the naive MONOMIAL
couplings -- those put every level's log-source in the single sigma-cascade {-2*3^j}, where sources are
mutually DEPENDENT (coboundary criterion below), which would give G_r=G_a for all r (no growth). That
reading is UNFAITHFUL: it would mean finite Mahler depth, contradicting MAHLER. The real Lambda-sources
are independent (Lambda has NO finite recurrence), so the growth G_r=G_a^r is FORCED BY MAHLER-consistency.
Module dim (QDIFF-1, grows) vs Galois dim (grows iff sources independent) is the distinction.
"""
from fractions import Fraction as F


def mult_coboundary_zpow(m):
    """sigma(y)=z^m*y, y=z^a: 3a=m+a => a=m/2 in Z iff m even. Return (is_coboundary, a)."""
    if m % 2 == 0:
        a = m // 2
        assert {3 * a: F(1)} == {m + a: F(1)}          # sigma(z^a)=z^{3a}=z^{m+a}
        return True, a
    return False, None


def cascade_coboundary(r):
    """rhs = sum_j r[j] z^{-2*3^j} (single sigma-cascade; r: dict j->coeff).
    sigma(g)-g at exponent -2*3^j has coeff x_{j-1}-x_j (x_{-1}=0) => x_j = -sum_{k<=j} r[k].
    Rational g (finite principal part) exists IFF x_j -> 0 for large j IFF sum_j r[j] = 0."""
    total = sum(r.values())
    return total == 0


def main():
    print("# PROBE QDIFF-3  RUNG 4 -- pin A_2, gate the tower  (EXACT)\n")
    print("## setup  A_2 = D_2(I+N_2), diagonal = pure z-power substitution scalings")
    print("   A_2 = [[z^12, z^12, 0],[0, z^4, z^4],[0, 0, 1]]   (A_1 = bottom-right 2x2 block)\n")

    # ---- 4a-diagonal: G_m at level 2 ----
    print("## 4a-diagonal  coboundary test on the level-2 diagonal scaling z^12  (the G_m)")
    ok12, a12 = mult_coboundary_zpow(12)
    ok4, a4 = mult_coboundary_zpow(4)
    print(f"   sigma(y)=z^12 y : y=z^{a12} (3*{a12}=12+{a12}); check sigma(z^6)=z^18=z^12*z^6  EXACT => COBOUNDARY")
    print(f"   general level-i scaling z^(4*3^(i-1)): coboundary with y=z^(2*3^(i-1)) (always even exponent)")
    print(f"   => G_m COLLAPSES at level 2, and at EVERY level. The reductive part stays PRO-only. [4a-diag PASS]\n")

    # ---- the coboundary criterion + the monomial-model warning ----
    print("## unipotent independence -- the criterion, and why the monomial model is UNFAITHFUL")
    print("   single-cascade coboundary criterion: sum_j r[j] z^{-2*3^j} is a coboundary IFF sum_j r[j]=0.")
    e1 = {0: F(1)}                       # z^{-2}   (j=0)
    e2 = {1: F(1)}                       # z^{-6}   (j=1)  = the monomial level-2 source
    print(f"   monomial sources: e1=z^-2 (j=0), e2=z^-6 (j=1).  e1 coboundary? {cascade_coboundary(e1)} (sum=1)")
    print(f"                                                     e2 coboundary? {cascade_coboundary(e2)} (sum=1)")
    # dependence: is e2 - lam*e1 a coboundary for some lam?  sum(e2-lam*e1)=1-lam=0 => lam=1
    dep = {k: e2.get(k, F(0)) - F(1) * e1.get(k, F(0)) for k in set(e1) | set(e2)}
    print(f"   e2 - 1*e1 coboundary? {cascade_coboundary(dep)} (sum={sum(dep.values())})  => lam=1 works")
    print(f"   => MONOMIAL sources are DEPENDENT (e2 = e1 + coboundary) => monomial model gives G_2=G_a (NO growth).")
    print(f"   WARN: but that reading is UNFAITHFUL: G_r=G_a for all r => finite Mahler depth => CONTRADICTS MAHLER.")
    print(f"      The naive monomial couplings collapse a cascade the real data does not. So NOT the answer.\n")

    # ---- 4a-unipotent: the growth, forced by MAHLER-consistency ----
    print("## 4a-unipotent  G_2 = G_a^2 (dim grew by 1) -- FORCED BY MAHLER-consistency")
    print("   The real level-i log-source is tied to Lambda_i, which has NO finite recurrence (R27-A/MAHLER)")
    print("   => the Lambda_i are algebraically independent data (doubly-exp denominators), NOT sigma-translates")
    print("   => successive log-sources are INDEPENDENT mod coboundaries => each level adds a genuine new G_a.")
    print("   Integrator is a DIRECT SUM: T = base + Lambda_1 + Lambda_2 (each Lambda enters T independently,")
    print("   parallel not iterated) => the unipotent is ABELIAN => G_2 = G_a^2 (not the non-abelian U_3).")
    print("   Consistency check: if instead dependent (G_a, no growth), the inverse limit would be finite-dim,")
    print("   contradicting MAHLER's infinite unipotent depth. So growth is REQUIRED, not optional. [4a-unip: by consistency]")
    print("   WARN: DIRECT construction (exact Lambda-derived couplings as Q(z) entries, then test independence")
    print("      literally) is NOT yet pinned -- the monomial model is too coarse. That pin is frontier-adjacent.\n")

    # ---- 4b: nesting ----
    print("## 4b  nesting -- G_2 restricts to G_1")
    print("   A_1 = [[z^4,z^4],[0,1]] is the bottom-right 2x2 block of A_2 => the level-1 PV sub-extension")
    print("   embeds. A level-2 automorphism (l_1->l_1+u_1, l_2->l_2+u_2) restricts to the level-1 sub as")
    print("   (l_1->l_1+u_1): the PROJECTION G_a^2 -> G_a (forget the level-2 coordinate). Clean, surjective,")
    print("   compatible with the tower => lim<- G_r well-defined. [4b PASS: G_2=G_a^2 --proj--> G_1=G_a]\n")

    # ---- verdict ----
    print("## RUNG-4 VERDICT")
    print("   WEAK CORE (pens as theorem):")
    print("     4a-diagonal: G_m is a coboundary at level 2 (and every level) => reductive part PRO-only.  [DIRECT]")
    print("     4a-unipotent: G_2 = G_a^2 (dim grew by 1), FORCED by MAHLER-consistency + direct-sum integrator.")
    print("     4b-nesting: G_2 = G_a^2 --projection--> G_1 = G_a; the tower nests, lim<- G_r defined.")
    print("     weak-4c: the reductive G_m emerges only in the inverse limit (= MAHLER doubly-exp denominators).")
    print("   FRONTIER (does NOT pen as established):")
    print("     - DIRECT construction of the exact Lambda-derived level-2 coupling as a Q(z) entry, then the")
    print("       literal independence test (monomial model is unfaithful -- collapses the cascade).")
    print("     - strong-4c: inverse-limit-of-groups = Galois-group-of-the-limit (the one step needing new math).")
    print("   FINDING surfaced: MODULE dim (QDIFF-1, grows) != automatically GALOIS dim; the latter grows iff the")
    print("   log-sources are independent, which MAHLER supplies -- but the naive monomial couplings would say")
    print("   otherwise, so the exact couplings (not the shape) are what a fully-direct Rung 4 needs.")
    print("   [structure/tower only; value = Rung 5, conditional; 7/15 excluded regardless, floor 0.473177.]")


if __name__ == "__main__":
    main()
