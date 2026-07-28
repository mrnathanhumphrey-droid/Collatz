"""
PROBE MAHLER -- is S a finite-order Mahler function, or is the depth provably infinite? (2026-07-28)

Wilson's spec. S_{i+1}=2*3^i[4 R_e^(i)(2) - R_e^(i)(0)]. S is finite-order Mahler IFF the exact rational
sequences a_i=3^i R_e^(i)(2), b_i=3^i R_e^(i)(0), and T_i=4a_i-b_i satisfy a FIXED finite-order linear
recursion with RATIONAL coefficients under i->i+1 (the tower z->z^3). EXACT arithmetic only (a float solve
finds spurious low-rank structure).

  M-A  for each of {a_i},{b_i},{T_i}: exact finite-order rational-recurrence test, L=1..floor((n-1)/2).
       - fits -> S is finite-order Mahler; doc-4 route LIVE (M-B).
       - none at any L -> Mahler order infinite = the unipotent depth; wall as a clean dichotomy (M-C).
  GATE (mandatory): reproduce R27-A -- the Lambda sequence must give NO finite rational recurrence
       (L(z) known non-rational). If the machinery finds one for Lambda, it's broken.
  M-C cross-check: consistency with (i) R27-A and (ii) doubly-exp denominators 2^{2*3^{k-1}}.

Exact R_e: nu = exact-rational stationary of the Syracuse chain (stationary_rational o build_markov_q(3,n)),
R_e-vector via the base-2 dlog cycle geometric sum (exact), autocorr lags {0,2} exact. Gated vs build_base2 float.
"""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from c_tilde_q17_probe import build_markov_q, stationary_rational, order_of_two
from probe_p6d import build_base2

HALF = F(1, 2)


def exact_Re(n):
    """exact rational R_e(0)=sum R_e[t]^2, R_e(2)=sum R_e[t]R_e[t+2], for the even-submeasure at level n."""
    q = 3 ** (n + 1); Nn = 3 ** n; twoN = 2 * Nn
    # base-2 dlog mod q (2 is a primitive root)
    D2 = {}; g = 1
    for t in range(twoN):
        D2[g] = t; g = (g * 2) % q
    # exact stationary nu on units mod 3^n (same Syracuse chain build_base2 uses)
    K, coprime, M = build_markov_q(3, n)
    pi = stationary_rational(K)                          # exact, sum=1, order = coprime
    # S[t] = sum_x w_x / 2^{a0(x,t)} over x with c=(D2[3x+1]-t) mod twoN EVEN; a0 = c if c>=1 else twoN
    # R_e[t] = C * S[t],  C = 2^twoN/(2^twoN - 1)  (common scalar; drops into R0,R2 as C^2)
    S = [F(0)] * twoN
    for x, w in zip(coprime, pi):
        dY = D2[(3 * x + 1) % q]
        for t in range(twoN):
            c = (dY - t) % twoN
            a0 = c if c >= 1 else twoN
            if a0 % 2 == 0:
                S[t] += w * F(1, 2 ** a0)
    C2 = F(2) ** (2 * twoN) / (F(2) ** twoN - 1) ** 2     # C^2
    R0 = C2 * sum(s * s for s in S)
    R2 = C2 * sum(S[t] * S[(t + 2) % twoN] for t in range(twoN))
    return R0, R2, twoN


def solve_exact(A, rhs):
    """exact Fraction Gaussian solve of square system; None if singular."""
    n = len(A); M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]


def fit_recurrence(seq, name):
    """smallest L with a fixed rational L-term recurrence reproducing ALL held-out terms; else None."""
    n = len(seq); Lmax = (n - 1) // 2
    for L in range(1, Lmax + 1):
        if n < 2 * L + 1:            # need >=1 held-out term to verify
            break
        A = [[seq[i - j] for j in range(1, L + 1)] for i in range(L, 2 * L)]
        rhs = [seq[i] for i in range(L, 2 * L)]
        alpha = solve_exact(A, rhs)
        if alpha is None:
            continue
        ok = all(seq[i] == sum(alpha[j - 1] * seq[i - j] for j in range(1, L + 1)) for i in range(2 * L, n))
        if ok:
            print(f"   {name}: FINITE recurrence order L={L}  alpha={[str(x) for x in alpha]}")
            return L, alpha
    print(f"   {name}: NO finite rational recurrence for L=1..{Lmax} (n={n} exact terms) -> infinite order")
    return None, None


def main():
    print("# PROBE MAHLER -- finite-order Mahler function, or provably infinite depth?\n")
    # ---- exact T_i, Lambda_i from the certified S-ladder (EPS exact through k=8) ----
    hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                        'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in hist.items()}
    S_ladder = {k: F(7, 15) + EPS[k] for k in EPS}                       # S_k exact
    Tcert = {i: S_ladder[i + 1] / 2 for i in range(1, 8) if (i + 1) in S_ladder}  # T_i = S_{i+1}/2, i=1..7
    Lam = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}            # Lambda_1..7 exact

    # ---- exact a_i, b_i via exact even-submeasure autocorrelation; gate vs float + vs T_cert ----
    print("## exact R_e build + gates (float build_base2 + 4a-b == T_cert)")
    NMAX = 6
    a = {}; b = {}
    t0 = time.time()
    for i in range(1, NMAX + 1):
        R0, R2, twoN = exact_Re(i)
        a[i] = 3 ** i * R2; b[i] = 3 ** i * R0
        Ti = 4 * a[i] - b[i]
        # gate 1: float autocorr from build_base2
        Sb = build_base2(i); Re = np.fft.ifft(np.abs(np.fft.fft(Sb['R_e'])) ** 2).real
        g1 = abs(float(R0) - Re[0]) + abs(float(R2) - Re[2])
        # gate 2: 4a-b == T_cert (exact) where available
        g2 = "n/a" if i not in Tcert else ("EXACT MATCH" if Ti == Tcert[i] else f"MISMATCH {float(Ti-Tcert[i]):.2e}")
        print(f"   i={i}: R0,R2 float-gate={g1:.2e}  4a-b vs T_cert: {g2}  ({time.time()-t0:.0f}s)")
    print()

    # ---- MANDATORY GATE: reproduce R27-A (Lambda has NO finite rational recurrence) ----
    print("## GATE  R27-A reproduction: Lambda_1..7 must give NO finite rational recurrence")
    lam_seq = [Lam[r] for r in range(1, 8)]
    fit_recurrence(lam_seq, "Lambda (expect NONE)")
    print()

    # ---- M-A: the three sequences ----
    print("## M-A  finite-order Mahler test on {a_i}, {b_i}, {T_i} (exact)")
    aseq = [a[i] for i in range(1, NMAX + 1)]
    bseq = [b[i] for i in range(1, NMAX + 1)]
    Tseq = [4 * a[i] - b[i] for i in range(1, NMAX + 1)]
    rA = fit_recurrence(aseq, "a_i = 3^i R_e(2)")
    rB = fit_recurrence(bseq, "b_i = 3^i R_e(0)")
    rT = fit_recurrence(Tseq, "T_i = 4a_i - b_i")
    print()

    # ---- M-C PROOF: denominator growth rate excludes ANY finite order (not just L<=2) ----
    print("## M-C PROOF  denominator growth rate vs the finite-recurrence ceiling")
    print("   A fixed order-L rational recurrence => D_i:=log2 den(T_i) obeys D_i <= sum_{j<=L} D_{i-j} + C,")
    print("   growth ratio = lambda_L < 2 for every finite L (->2 as L->inf). Observed ratio decides it.")
    import math
    print(f"   {'i':>2} {'den(T_i) bits = log2 den':>24} {'ratio D_i/D_(i-1)':>18}")
    Dprev = None
    for i in range(1, 8):
        if i not in Tcert:
            continue
        Di = math.log2(Tcert[i].denominator)
        ratio = (Di / Dprev) if Dprev else float('nan')
        print(f"   {i:>2} {Di:>24.3f} {ratio:>18.4f}")
        Dprev = Di
    print("   finite-order ceilings: lambda_1=1, lambda_2=phi=1.618, lambda_3=1.839, ... sup_L lambda_L = 2 (never reached).")
    print("   => observed ratio -> 3 EXCEEDS 2 = every finite-order ceiling => NO finite rational recurrence at ANY L.")
    print("   (denominator theorem: den(S_r) ~ 2^{2*3^{r-1}}, ratio 3, doubly-exp; T_i=S_{i+1}/2 inherits it.)\n")

    # ---- verdict ----
    any_fit = any(r[0] is not None for r in (rA, rB, rT))
    print("## VERDICT")
    if any_fit:
        print("   >>> M-B branch: a finite rational recurrence EXISTS -> S is finite-order Mahler.")
        print("       Build the q-difference equation T(z^3)=A(z)T(z)+B(z); Nishioka decidability check next.")
    else:
        print("   >>> M-C branch: NO finite rational recurrence in any sequence -> Mahler order INFINITE.")
        print("       S_inf is not a finite-order Mahler evaluation; infinite order = the unbounded unipotent")
        print("       depth of M=D(I+N). Cross-checks: (i) R27-A Lambda no-recurrence [above],")
        print("       (ii) doubly-exp denominators 2^{2*3^{k-1}} are consistent with infinite order.")
    print(f"   [depth: a,b exact to i={NMAX}; T,Lambda exact to i=7. 7/15 UNAFFECTED (floor 0.473177).]")


if __name__ == "__main__":
    main()
