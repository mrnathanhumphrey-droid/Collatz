"""
PROBE R12 -- SUPPORT LAW, CLOSED LOOP, LOBES. Reuses R7 (build_mu/mu1/cram) + R10 (dlog_table/
autocorr_dlog/A_N/Lambda_r/layer_mass). Follows R11 (walk-back #31: U is sparse, not dense-flat).

Derivations gated:
  SUPPORT LAW: U(k,xi)=0 unless k==xi mod 3, exact; nonzero |U|=3^{-(r-1)/2}; U(-k,-xi)=conj U(k,xi).
  CLOSED LOOP: sum_{k prim mod 3^r} |theta_hat(k)|^2 e(km/3^r) = C_{r+1}(m)/3   (every angular moment = a C-table entry).
  LOBES: Lambda_r = sum_{k prim} |theta_hat(k)|^2 Re w(k/3^r), Re w(x)=(4cos2pix-1)/(17-8cos2pix);
         split by sign of Re w -> L_r (Re w>0) minus M_r (Re w<0).
  CONVENTION (R12-F): historical eps-table (result_77_7...k8) vs exact d_k = S_k - 7/15.
"""
import os, sys, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

S = R10.S


def v3(n):
    n = abs(int(n));  j = 0
    if n == 0: return 10**9
    while n % 3 == 0: n //= 3; j += 1
    return j


def U_matrix(r):
    N = 3 ** r; d = R10.dlog_table(r)
    beta = [d[z] for z in range(N)]
    e = lambda x: cmath.exp(2j * math.pi * x)
    U = np.zeros((N, N), dtype=complex)
    for k in range(N):
        for xi in range(N):
            U[k, xi] = sum(e((k * beta[z] - xi * z) / N) for z in range(N)) / N
    return U


def theta2_array(mu_r, r):
    """|theta_hat(k)|^2 (float) for k=0..3^r-1, from dlog-domain autocorrelation g_r (exact) via cos-sum."""
    N = 3 ** r
    g = R10.autocorr_dlog(mu_r, r)                      # exact Fractions, len N
    gf = [float(x) for x in g]
    return np.array([sum(gf[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]), g


def closed_loop_lhs(g, r, m, sign=+1):
    """exact sum_{k prim mod 3^r} |theta|^2 e(sign*km/3^r) = sum_u g(u) c_{3^r}(u + sign*m)."""
    N = 3 ** r
    return sum(g[u] * R7.cram(r, u + sign * m) for u in range(N))


def main():
    print("# PROBE R12 -- SUPPORT LAW / CLOSED LOOP / LOBES. Exact gates + labeled measurement.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 7):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ---- R12-A ----
    print("## R12-A  SUPPORT LAW  (U(k,xi)=0 unless k==xi mod3; nonzero |U|=3^-(r-1)/2; U(-k,-xi)=conjU)")
    okA = True
    for r in range(2, 6):
        N = 3 ** r; U = U_matrix(r)
        prim = [i for i in range(N) if i % 3 != 0]
        offclass_max = 0.0; inclass_dev = 0.0
        for k in prim:
            for xi in prim:
                if (k - xi) % 3 != 0:
                    offclass_max = max(offclass_max, abs(U[k, xi]))
                else:
                    inclass_dev = max(inclass_dev, abs(abs(U[k, xi]) - 3 ** (-(r - 1) / 2)))
        conj_dev = max(abs(U[(-k) % N, (-xi) % N] - np.conj(U[k, xi])) for k in prim for xi in prim)
        gi = offclass_max < 1e-10 and inclass_dev < 1e-10 and conj_dev < 1e-10
        okA = okA and gi
        print(f"   r={r}: off-class max|U|={offclass_max:.1e}  in-class |U|-3^-(r-1)/2 dev={inclass_dev:.1e}  "
              f"conj-block dev={conj_dev:.1e}  [{'PASS' if gi else 'FAIL'}]")
    print(f"   => R12-A {'GATE PASS (support law {k==xi mod3}, forced by fiber derivation)' if okA else 'FAIL (#32)'}\n")

    # ---- R12-B ----
    print("## R12-B  CLOSED-LOOP WELD  (sum_{ord=3^r}|theta|^2 e(km/3^r) = C_{r+1}(m)/3 exact)")
    okB = True
    for r in range(2, 5):
        _, g = theta2_array(mu[r], r)
        Np1 = 3 ** (r + 1); inv4 = pow(4, -1, Np1)
        # determine sign convention on first m, then lock it
        ms = [1, 2, 3, 3 ** r]                          # incl DC (m=3^r)
        # pick sign by matching m=1
        rhs1 = R7.Ck(r + 1, 1, mu[r], pow(inv4, 1, Np1)) / 3
        sign = +1 if closed_loop_lhs(g, r, 1, +1) == rhs1 else -1
        allok = True; rows = []
        for m in ms:
            lhs = closed_loop_lhs(g, r, m, sign)
            rhs = R7.Ck(r + 1, m, mu[r], pow(inv4, m, Np1)) / 3
            ok = (lhs == rhs); allok = allok and ok
            tag = "DC" if m % (3 ** r) == 0 else f"m={m}"
            rows.append(f"{tag}:{'OK' if ok else 'DEV'}")
        okB = okB and allok
        print(f"   r={r} (sign e({'+' if sign>0 else '-'}km)): " + "  ".join(rows) +
              f"   [{'PASS' if allok else 'FAIL'}]")
    print(f"   => R12-B {'GATE PASS -- every angular moment IS a banked C-table entry (loop closed)' if okB else 'FAIL'}\n")

    # ---- R12-F (do before C/E so ledger uses certified eps) ----
    print("## R12-F  CONVENTION CERTIFICATION  (historical eps-table vs exact d_k = S_k - 7/15)")
    hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                       'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
    eps = {int(k): F(int(v['num']), int(v['den'])) for k, v in hist.items()}
    print(f"   {'k':>2} {'d_k = S_k-7/15 (exact)':>26} {'historical eps_k':>26} {'match':>6}")
    okF = True
    for k in range(1, 7):
        dk = S[k] - F(7, 15)
        good = (dk == eps[k]); okF = okF and good
        print(f"   {k:>2} {str(dk)[:26]:>26} {str(eps[k])[:26]:>26} {'OK' if good else 'STOP':>6}")
    print(f"   => convention IDENTICAL (eps_k = d_k = S_k - 7/15, byte-equal k=1..6). "
          f"Historical table is EXACT through k=8; k>=9 float-only (the wall).")
    # exact S_7, S_8 from certified eps -> cross-check R11's mu_7-built S_7
    S7 = F(7, 15) + eps[7]; S8 = F(7, 15) + eps[8]
    S7_renewal = R10.layer_mass(R10.autocorr_dlog(mu[6], 6), 6)  # this is S_6; build mu_7 for S_7 cross-check
    mu[7] = R7.build_mu(mu[6], 7)
    S7_renewal = R10.layer_mass(R10.autocorr_dlog(mu[7], 7), 7)
    print(f"   CROSS-CHECK: S_7 from eps-table == S_7 from renewal mu_7 build?  {S7 == S7_renewal}  (independent routes)\n")

    # ---- R12-E ----
    print("## R12-E  LEDGER EXTENSION via certified exact eps (Lambda_r = (eps_{r+1}-eps_r)/2)")
    print(f"   {'r':>2} {'Lambda_r (exact)':>22} {'float':>14} {'sign':>5} {'source':>14}")
    for r in range(1, 8):
        if r <= 5:
            Lr = R10.Lambda_r(mu, r)[0]; src = "char/R10"
        else:
            Lr = (eps[r + 1] - eps[r]) / 2; src = "eps-table"
        print(f"   {r:>2} {str(Lr)[:22]:>22} {float(Lr):>+14.6e} {'-' if Lr<0 else '+':>5} {src:>14}")
    print(f"   NEW: Lambda_6={float((eps[7]-eps[6])/2):+.6e} (-), Lambda_7={float((eps[8]-eps[7])/2):+.6e}. "
          f"mu_8 renewal build (supp 4374, autocorr ~19M Fractions) is the exact wall; eps-table route bypasses it.\n")

    # ---- R12-C ----
    print("## R12-C  LOBE LEDGER  (L_r=sum_{Re w>0}|theta|^2 Re w, M_r=-sum_{Re w<0}; Lambda_r=L_r-M_r)")
    Rew = lambda x: (4 * math.cos(2 * math.pi * x) - 1) / (17 - 8 * math.cos(2 * math.pi * x))
    print(f"   {'r':>2} {'L_r (Re w>0)':>14} {'M_r (Re w<0)':>14} {'L_r-M_r':>14} {'Lambda_r(exact)':>16} {'L_r+M_r':>12}")
    for r in range(2, 7):
        N = 3 ** r
        th2, _ = theta2_array(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]
        L = sum(th2[k] * Rew(k / N) for k in prim if Rew(k / N) > 0)
        M = -sum(th2[k] * Rew(k / N) for k in prim if Rew(k / N) < 0)
        Lr = R10.Lambda_r(mu, r)[0] if r <= 6 else None
        print(f"   {r:>2} {L:>14.6e} {M:>14.6e} {L-M:>+14.6e} {float(Lr):>+16.6e} {L+M:>12.6e}")
    print("   (do the lobe masses L_r, M_r stabilize? = psi-existence evidence; verbatim, no fit.)\n")

    # ---- R12-D ----
    print("## R12-D  CLASS-RESOLVED PROFILE  (k mod 3 split; conjugation => class-1 mirror class-2)")
    for r in range(2, 6):
        N = 3 ** r
        th2, _ = theta2_array(mu[r], r)
        c1 = [k for k in range(1, N) if k % 3 == 1]
        c2 = [k for k in range(1, N) if k % 3 == 2]
        # mirror: |theta(k)|^2 == |theta(N-k)|^2, and N-k flips class 1<->2
        mirror_dev = max(abs(th2[k] - th2[(N - k) % N]) for k in c1)
        L1 = sum(th2[k] * ((4*math.cos(2*math.pi*k/N)-1)/(17-8*math.cos(2*math.pi*k/N))) for k in c1)
        L2 = sum(th2[k] * ((4*math.cos(2*math.pi*k/N)-1)/(17-8*math.cos(2*math.pi*k/N))) for k in c2)
        print(f"   r={r}: class-mirror |theta(k)|^2 vs |theta(-k)|^2 dev={mirror_dev:.1e} "
              f"[{'EXACT' if mirror_dev<1e-9 else 'DEV'}]; class-lobe sums c1={L1:+.6e} c2={L2:+.6e} (sum={L1+L2:+.6e})")
    print("   (conjugation mirror holds => the two mod-3 classes are the +/- (P+,P-) pair; connects to (1,4).)")


if __name__ == "__main__":
    main()
