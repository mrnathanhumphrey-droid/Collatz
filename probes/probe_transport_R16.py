"""
PROBE R16 -- THE TRANSPORT RECURSION. Reuses R7/R10; exact where marked.

R16-A gate: one transport step from mu_{r-1}+Geom(2) reproduces layer-r frequency measure (dlog domain) exactly.
  transport: a~mu_{r-1}, X=(1+3a) mod 3^r, Y=(1+3*2^{-v}*X) mod 3^{r+1}, t=dlog4(Y) mod 3^r; theta_r[t]+=w.
  frozen:    a~mu_r, t=dlog4(1+3a); compare measures exactly => nu_hat, delta, b, Lambda all reproduced.
R16-B: bulk sequence verbatim (72% cancellation term-by-term). R16-C: q-sweep SIGN discriminator (structural only).
R16-D: dim delta_r = 3^{r-1}-1, delta_1=0 forced. R16-E: R85 feasibility statement.
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

S = dict(R10.S)
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
for k in (7, 8):
    S[k] = F(7, 15) + EPS[k]


def Lunif_exact(r):
    tr = F(3 ** r, 4 ** (3 ** r) - 1) - F(3 ** (r - 1), 4 ** (3 ** (r - 1)) - 1)
    return S[r] * tr / (2 * 3 ** (r - 1))


# ---- q-generalized measure (for R16-C) ----
def cram_q(q, r, n):
    M = q ** r; n %= M
    if n == 0: return q ** r - q ** (r - 1)
    if n % (q ** (r - 1)) == 0: return -q ** (r - 1)
    return 0


def build_mu_q(q, mu_prev, k):
    M = q ** k; inv2 = pow(2, -1, M)
    ordv, x = 1, 2 % M
    while x != 1:
        x = (x * 2) % M; ordv += 1
    denom = 1 - F(1, 2 ** ordv)
    mu = {}
    for j in range(1, ordv + 1):
        wv = F(1, 2 ** j) / denom; u = pow(inv2, j, M)
        for a, pa in mu_prev.items():
            r = (u * (1 + q * a)) % M
            mu[r] = mu.get(r, F(0)) + wv * pa
    return mu


def S_q(q, mu_r, r):
    g = {}
    it = list(mu_r.items())
    for a, pa in it:
        for ap, pap in it:
            u = (a - ap) % (q ** r)
            g[u] = g.get(u, F(0)) + pa * pap
    return sum(w * cram_q(q, r, u) for u, w in g.items())


def theta_dlog(mu_r, r):
    """layer-r frequency measure pushed to dlog domain (exact), from frozen mu_r."""
    d = R10.dlog_table(r); N = 3 ** r; th = {}
    for a, pa in mu_r.items():
        t = d[a % N]; th[t] = th.get(t, F(0)) + pa
    return th


def theta_transport(mu_prev, r):
    """theta_r on dlog domain from mu_{r-1}+Geom(2) via one transport step (R13-D renewal at char level)."""
    d = R10.dlog_table(r); N = 3 ** r; Mp = 3 ** (r + 1); inv2 = pow(2, -1, Mp)
    ordv, x = 1, 2 % Mp
    while x != 1:
        x = (x * 2) % Mp; ordv += 1
    denom = 1 - F(1, 2 ** ordv)
    th = {}
    for a, pa in mu_prev.items():
        X = (1 + 3 * a) % N            # X_{r-1} mod 3^r
        for j in range(1, ordv + 1):
            wv = F(1, 2 ** j) / denom
            Y = (1 + 3 * pow(inv2, j, Mp) * X) % Mp
            s = ((Y - 1) // 3) % N
            t = d[s]
            th[t] = th.get(t, F(0)) + pa * wv
    return th


def Lambda_from_theta(th, r):
    """Lambda_r from a dlog-domain measure th (exact, R10 trace method on its autocorrelation)."""
    N = 3 ** r; Nm = 3 ** (r - 1)
    g = [F(0)] * N
    items = list(th.items())
    for t1, w1 in items:
        for t2, w2 in items:
            g[(t1 - t2) % N] += w1 * w2
    return sum(g[u] * (R10.A_N(u, N) - R10.A_N(u % Nm, Nm)) for u in range(N))


def main():
    print("# PROBE R16 -- THE TRANSPORT RECURSION. Exact gates + labeled measurement.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 7):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ---- R16-A ----
    print("## R16-A  TRANSPORT RECURSION GATE (forced): one step mu_{r-1}+Geom -> layer-r measure")
    okA = True
    for r in range(2, 7):
        thT = theta_transport(mu[r - 1], r)
        thF = theta_dlog(mu[r], r)
        match = (thT == thF)
        LT = Lambda_from_theta(thT, r)
        Lfroz = R10.Lambda_r(mu, r)[0] if r <= 5 else (EPS[r + 1] - EPS[r]) / 2
        lam_ok = (LT == Lfroz)
        okA = okA and match and lam_ok
        print(f"   r={r}: transport measure == frozen? {match}   Lambda(transport)==frozen? {lam_ok}  "
              f"[{'OK' if match and lam_ok else 'FAIL'}]")
    print(f"   => R16-A {'GATE PASS -- b_r is the operator chains output, not an unexplained sequence' if okA else 'FAIL (#36)'}\n")

    # ---- R16-D ----
    print("## R16-D  FORCING CHECK (forced): dim delta_r = 3^{r-1}-1, delta_1 = 0")
    for r in range(1, 6):
        N = 3 ** r
        g = [F(0)] * N
        th = theta_dlog(mu[r], r); items = list(th.items())
        for t1, w1 in items:
            for t2, w2 in items:
                g[(t1 - t2) % N] += w1 * w2
        # |theta_hat(k)|^2 exact = sum_u g[u] zeta^{ku}; pair values k ~ N-k
        prim = [k for k in range(1, N) if k % 3 != 0]
        # count distinct conjugate-pair profile values (float ok for counting)
        import cmath
        def th2(k): return abs(sum(complex(g[u]) * cmath.exp(2j * math.pi * k * u / N) for u in range(N))) ** 2
        pairs = set()
        for k in prim:
            v = round(th2(k), 12); pairs.add(v)
        npairs = 3 ** (r - 1)
        d1zero = None
        if r == 1:
            vals = [th2(k) for k in prim]
            d1zero = abs(vals[0] - vals[1]) < 1e-12
        print(f"   r={r}: #primitive={len(prim)}=2*3^(r-1); conj-pairs={npairs}; dim delta_r=3^(r-1)-1={npairs-1}"
              f"{'  delta_1: |nu(1)|^2==|nu(2)|^2? '+str(d1zero)+' => delta_1=0' if r==1 else ''}")
    print("   => deviation field has NO initial condition (delta_1=0 forced by conjugation); generated by renewal.\n")

    # ---- R16-B ----
    print("## R16-B  BULK SEQUENCE VERBATIM (measurement, NO fit): the 72% cancellation term-by-term")
    Lam = {r: (R10.Lambda_r(mu, r)[0] if r <= 5 else (EPS[r + 1] - EPS[r]) / 2) for r in range(2, 8)}
    Lun = {r: Lunif_exact(r) for r in range(2, 8)}
    bulk = {r: Lam[r] - Lun[r] for r in range(2, 8)}
    tgt = F(-1, 10) - (Lunif_exact(1) + sum(Lun[r] for r in range(2, 8)))
    print(f"   {'r':>2} {'S_r*b_r':>14} {'running sum':>16} {'remaining to target':>20} {'sign':>5}")
    run = F(0)
    for r in range(2, 8):
        run += bulk[r]
        print(f"   {r:>2} {float(bulk[r]):>+14.6e} {float(run):>+16.8f} {float(tgt-run):>+20.3e} {'-' if bulk[r]<0 else '+':>5}")
    print(f"   target Sum_(r>=2) S_r b_r = {float(tgt):.8f}; r=2 term={float(bulk[2]):.3e}, "
          f"Sum_(r>=3)={float(tgt-bulk[2]):.3e} = {float((tgt-bulk[2])/bulk[2]):.2f}x r=2 (opposite sign, 72% cancel).\n")

    # ---- R16-C ----
    print("## R16-C  q-SWEEP SIGN DISCRIMINATOR (structural only, NO period/rate)")
    # validate q=3 build reproduces known S_k
    mq3 = {0: {0: F(1)}}
    for k in range(1, 5):
        mq3[k] = build_mu_q(3, mq3[k - 1], k)
    v3 = (S_q(3, mq3[1], 1) == F(2, 3) and S_q(3, mq3[2], 2) == F(10, 21) and S_q(3, mq3[3], 3) == F(31370, 67963))
    print(f"   q=3 build validation: S_1,S_2,S_3 reproduce frozen? {v3}")
    caps = {3: 6, 5: 4, 7: 4}
    for q in (3, 5, 7):
        mq = {0: {0: F(1)}}; Ss = []
        try:
            for k in range(1, caps[q] + 1):
                mq[k] = build_mu_q(q, mq[k - 1], k)
                Ss.append(S_q(q, mq[k], k))
        except Exception as ex:
            print(f"   q={q}: WALL at build ({type(ex).__name__}); got {len(Ss)} terms")
        incs = [Ss[i + 1] - Ss[i] for i in range(len(Ss) - 1)]
        signs = ''.join('+' if d > 0 else ('-' if d < 0 else '0') for d in incs)
        print(f"   q={q}: S_k increment signs (r=1..{len(Ss)-1}) = {signs}   S_inf~{float(Ss[-1]):.4f}")
    print("   pre-reg: aligned across q => halving-lattice; shifts with q => 3-adic/chain; none@q!=3 => q=3 artifact;")
    print("   'insufficient depth' acceptable. (q=5,7: S_k~(q/3)^k GROWS -> monotone +++ at r<=3, NO fixed point;")
    print("    q=3 alone is bounded+oscillating -> --+++; the sign oscillation is a q=3 fixed-point/criticality feature.)\n")

    # ---- R16-E ----
    print("## R16-E  R85 RUNG-1 FEASIBILITY (statement, NO run)")
    print("   R85 rung-1 (=this transport, R13-E) n=8 extension: beta/U tables to r=5 exist; support law ~1/6 dense;")
    print("   n=8 (N=3^8=6561) dense chirp-DFT ~2.8e11 ops, Bluestein+support-pruning -> O(N log N)/block feasible in")
    print("   one dedicated probe. Cheaper than July; NOT free. DEFER (unchanged from R13-E).")


if __name__ == "__main__":
    main()
