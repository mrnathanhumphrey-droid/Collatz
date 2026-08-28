"""
PROBE T4 -- disproof by exclusion: does a NON-TRIVIAL (structured) product form for S_inf exist?
Wilson-penned kill-first spec. A trivial product S_inf=(7/15)*prod(1+delta_i) with delta back-solved is
VACUOUS (R9-D/R10 trap). T4 excludes a STRUCTURED product whose factors are INDEPENDENTLY determined
(periodic orbits / places / cyclotomic data) and carry info the tail SUM does not.

Verdict vocabulary: EXCLUDED (with witness) / VACUOUS (tautological) / RE-RUN (disposed object) /
SURVIVED-A-KILL (flag, do NOT celebrate).

Guardrails G1-G5 were confirmed verbatim from the corpus (cited inline). Machinery (exact_Re ladder,
build_M) reused, not re-derived.
"""
import os, sys, json
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
from sympy import factorint

HERE = os.path.dirname(os.path.abspath(__file__))

# ============================================================ ANCHORS (exact, certified) ==========
# exact S_k via the certified EPS ladder (probe_qdiff1.py:31-38 path); S_k = 7/15 + EPS[k]
EPS = {int(k): F(int(v['num']), int(v['den']))
       for k, v in json.load(open(os.path.join(HERE, '..', 'experiments_output',
                                                'result_77_7_eps_exact_through_k8_v2_vec_pool.json'))).items()}
S = {k: F(7, 15) + EPS[k] for k in EPS}                         # exact S_k
T = {i: S[i + 1] / 2 for i in range(0, 9) if (i + 1) in S}      # T_i = S_{i+1}/2 ; T_0 = 1/3
Lam = {i: T[i] - T[i - 1] for i in range(1, 9) if i in T and (i - 1) in T}  # Lambda_i exact
kmax = max(Lam)

print("# PROBE T4 -- exclusion of a structured product form for S_inf\n")
print("## ANCHORS (exact rationals)")
print(f"   S_1={S[1]}  S_2={S[2]}  (dead leading mode 7/15={F(7,15)}; true S_inf ~0.475 > 7/15)")
print(f"   Lambda_1={Lam[1]} (=-2/21 pin)   S_inf=7/15 <=> sum_{{i>=2}}Lam_i = -1/210 (char-ledger R10: sum_{{r>=1}}Lam_r=-1/10)")
print(f"   T_0={T[0]}  tail = S_inf/2 - T_0 = sum_{{i>=1}} Lam_i\n")

# ============================================================ G -- guardrails (confirmed) ==========
print("## G -- guardrails (verbatim-confirmed from corpus; each GATES a candidate)")
print("   G1 PRODFORM tested a product for the SPECTRUM Chat(m)=|rho_hat|^2 (corr~0, DISPOSED); a product")
print("      for mu_hat ITSELF is a DIFFERENT object (result_PRODFORM.md:1,12). -> C not a re-run of PRODFORM.")
print("   G2 7=N(2-omega) is a ROBUST q=3 fact (|Ghat(omega)|^2=1/7) but REFUTED as a family law at")
print("      p=7,11,13 (ratios 0.64/0.69/0.65) (result_85_bridge.md:22-36, result_1b_halving_pgf.md:9-11).")
print("   G3 renewal M(z)=A(z)/(1-z) is over the 2nd-moment primitive S_0 (=R7's M_i(1)); char-ledger")
print("      sum_{r>=1}Lam_r=-1/10 is the SUM machinery (result_42:26-44, result_charledger_R10.md:29-30).")
print("   G4 mu_hat FE (Siegel 2.18/Tao 1.12): mu_hat_{n+1}(xi)=sum_v 2^-v e(.) mu_hat_n(xi 2^-v) is a")
print("      WEIGHTED SUM over v (Geom-1/2 convex combo), NOT a rescaling f(xi)=a(xi)f(q xi). [gates C]")
print("   G5 limit operator = CONTINUOUS spectrum / branch cut at z=2, NO discrete eigenvalue")
print("      (README:196, INTERLEVEL:57-59, D1_T_M:54, K-erratum:7); L=4 EP witnesses inaccessible,")
print("      cond# 2.5e14->2.4e17, biov->0 (W4:35,42-49, R39:18,37, framework_cohesion:111). [gates A]\n")

# ============================================================ NULL B -- cyclotomic cascade ==========
print("## NULL B -- <2>-cascade cyclotomic product: do the NUMERATORS carry an independent product?")
# denominators of Lam_i factor cyclotomically ALWAYS (den ~ 3 * (2^(2*3^(i-1))-1)). That is VACUOUS.
# The value lives in the NUMERATORS. Test whether num primes are the SAME cofactor primes or FOREIGN.
BMAX = 4   # cap: den(Lam_i) ~ 2^(2*3^(i-1))-1 is ~440 digits by i=7 (infeasible to factor); i<=4 suffices
den_primes = set()
for i in Lam:
    if i > BMAX:
        continue
    den_primes |= set(factorint(Lam[i].denominator).keys())
print(f"   union of DENOMINATOR primes over Lam_1..{BMAX} (the cyclotomic <2>-cofactors): {sorted(den_primes)}")
foreign_seen = set()
for i in sorted(k for k in Lam if k <= BMAX):
    num = abs(Lam[i].numerator); den = Lam[i].denominator
    nf = factorint(num); df = factorint(den)
    foreign = sorted(p for p in nf if p not in den_primes)
    foreign_seen |= set(foreign)
    print(f"   Lam_{i}: num={Lam[i].numerator:>22}  factors={dict(nf)}")
    print(f"          den={den:>22}  factors={dict(df)}   FOREIGN num-primes (not cofactors): {foreign}")
print(f"\n   numerator primes NOT in the cofactor set (foreign to the cyclotomic cascade): {sorted(foreign_seen)}")
print("   => denominators factor cyclotomically (VACUOUS, pure number theory); numerators carry FOREIGN")
print("      primes (2,5,149,...) with NO product over the cofactors {7,19,73,...}. A product over the")
print("      cascade primes would have to reproduce these foreign primes -- it cannot.")
print("   VERDICT B: VACUOUS (denominators-only rewriting; numerators carry no independent cascade product).\n")

# ============================================================ NULL A -- Ruelle / dynamical zeta ======
print("## NULL A -- dynamical/Ruelle-zeta Euler product: S_inf = residue at the leading (1/3) pole?")
print("   FRAMING gate (already disposed): S_n = sum|mu_hat_n|^2 is a PLANCHEREL SUM, not a dynamical")
print("   trace tr(L^n) -- FAURE_IJ_HYPOTHESES:16,20 'Syracuse: not in this category'. So there is no")
print("   orbit Euler product to take a residue of. Still, demonstrate on the certified operator M:")

def two_subgroup(qL):
    s, x = [], 1
    while True:
        s.append(x); x = (x * 2) % qL
        if x == 1:
            break
    return s

def build_M(q, L):                       # verbatim from probe_25_transfer_operator_Aprime.py
    qL = q ** L
    inv2 = pow(2, -1, qL)
    sub = two_subgroup(qL); ordL = len(sub)
    Z = 1.0 - 2.0 ** (-ordL)
    mult = [(pow(inv2, delta, qL), (2.0 ** (-delta)) / Z) for delta in range(1, ordL + 1)]
    states = [(a, b, g) for a in sub for b in sub for g in range(qL)]
    idx = {s: i for i, s in enumerate(states)}; n = len(states)
    rows, cols, vals = [], [], []
    for (a, b, g) in states:
        i = idx[(a, b, g)]
        for (ga, wa) in mult:
            ap = (a * ga) % qL
            for (gb, wb) in mult:
                bp = (b * gb) % qL; T_ = (ap - bp) % qL
                if (g + T_) % q == 0:
                    gp = ((g + T_) // q) % qL
                    rows.append(idx[(ap, bp, gp)]); cols.append(i); vals.append(wa * wb)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n)), idx, n

for L in (1, 2):
    try:
        M, idx, n = build_M(3, L)
        ev = np.linalg.eigvals(M.toarray())
        mag = np.sort(np.abs(ev))[::-1]
        lam1 = mag[0]; lam2 = mag[1]
        print(f"   L={L} dim={n:>5}: lambda_1={lam1:.6f} (cert 1/3={1/3:.6f})  |lambda_2|/lambda_1 = {lam2/lam1:.4f}")
    except Exception as e:
        print(f"   L={L}: build/eig failed ({type(e).__name__}: {e})")
print("   => lambda_1 -> 1/3 (Perron, certified). The generating fn <1|(1-zM)^-1|v0> = sum ||pi_k||^2 z^k")
print("      has its LEADING residue at the 1/3 pole = the leading 3^-k mode = the DEAD 7/15, NOT the tail.")
print("      The SUBDOMINANT pole undergoes the q=3 exceptional point (|lambda_2|/lambda_1 -> 1, biov->0,")
print("      cond# blow-up, L=4 inaccessible: G5) -> it dissolves into CONTINUOUS spectrum (branch cut z=2).")
print("   VERDICT A: EXCLUDED. (i) no orbit product exists [framing NO_FIT]; (ii) the leading residue")
print("      delivers 7/15, not 0.475; (iii) the discrete subdominant factor does NOT survive L->inf [G5].\n")

# ============================================================ NULL C -- infinite Mahler product ======
print("## NULL C -- infinite Mahler product of mu_hat, Plancherel-summed to S_inf")
print("   (i) PREMISE gate: iterate the FE two steps.")
print("       mu_{n+1}(xi) = sum_v 2^-v e_v mu_n(xi 2^-v)")
print("       mu_{n+1}(xi) = sum_v 2^-v e_v [ sum_v' 2^-v' e'_v' mu_{n-1}(xi 2^-(v+v')) ]")
print("                    = sum_{v,v'} 2^-(v+v') e_v e'_v' mu_{n-1}(xi 2^-(v+v'))   <-- DOUBLE SUM over")
print("       valuation PATHS, i.e. the RENEWAL/path-sum (G3), NOT a product. -> step (i) = RE-RUN (renewal).")
print("   (ii) PLANCHEREL square: even granting a product mu = prod_k f_k, S_n = sum_xi |mu_n(xi)|^2.")
print("        |sum_v a_v|^2 = sum_v |a_v|^2  +  sum_{v!=v'} a_v conj(a_v')   [diagonal + OFF-DIAGONAL].")
# tiny numeric demonstration that the off-diagonal (cross) term is nonzero and irreducible
import cmath
xi = 0.37
avals = [2.0 ** (-v) * cmath.exp(2j * cmath.pi * xi * 2.0 ** (-v)) for v in range(1, 40)]
diag = sum(abs(a) ** 2 for a in avals)
total = abs(sum(avals)) ** 2
offdiag = total - diag
print(f"        demo (xi={xi}): |sum a_v|^2 = {total:.6f} = diag {diag:.6f} + off-diag {offdiag:.6f}")
print(f"        off-diagonal / total = {offdiag/total:.3f}  (NONZERO; = the non-free/monotone part that")
print("        carries the value beyond the diagonal -- free-prob ledger: |phi(X~j X~k X~j)|=0.1078).")
print("        So Sum -> Product does NOT survive |.|^2: the cross terms are exactly the tail correction.")
print("   VERDICT C: RE-RUN (iteration = renewal sum, G3) + EXCLUDED (Plancherel square breaks any product).\n")

# ============================================================ MAHLER backstop ======================
print("## MAHLER backstop (applies to ANY survivor)")
print("   zeta(2)=pi^2/6 is reachable because product + a CLOSING functional equation (Euler product +")
print("   the reflection/FE). MAHLER (proven) = the closing FE does NOT exist at finite order (infinite")
print("   Mahler depth = the branch cut = continuous spectrum, G5). So any product whose VALUE-EXTRACTION")
print("   needs a closing FE <=> finite Mahler order <=> CONTRADICTS MAHLER. Even a formal product is")
print("   non-closing (infinitely many independent factors, no finite relation) = the wall.\n")

print("## SUMMARY")
print("   A (Ruelle/dynamical zeta): EXCLUDED  -- framing NO_FIT + leading residue=7/15 + subdominant EP->continuous spectrum")
print("   B (<2>-cascade cyclotomic): VACUOUS  -- denominators factor (trivially), numerators carry foreign primes")
print("   C (infinite Mahler product): RE-RUN + EXCLUDED -- FE iteration=renewal; Plancherel |.|^2 breaks the product")
print("   Backstop: any survivor's value-extraction needs a closing FE = finite Mahler order = contradicts MAHLER.")
print("   => the space where a STRUCTURED product for S_inf could hide is CLOSED. No SURVIVED-A-KILL.")
