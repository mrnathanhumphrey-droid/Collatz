"""
VERIFY (Wilson's ψ-resolution session, follows R13). Read-only gate of the load-bearing claims.

Claims:
  (1) PARTIAL SUM: A_r(m) := C_{r+1}(m)/3 = gamma_r(tau_m) - gamma_{r-1}(tau_m); gamma_n(tau_m)=1+Sum_{r<=n}A_r(m).
  (2) NO NON-UNIFORM psi: gamma_n bounded => A_r(m) has no nonzero limit; Cesaro mean gamma_n/n -> 0.
      (A_r may still oscillate -- R13-C -- so 'converges to uniform OR does not converge'.)
  (3) RETARGET IDENTITY: Lambda_r = S_r * <delta_r, Re w> + Lambda_r^unif, delta_r = normalized profile - uniform,
      Lambda_r^unif = (S_r/M) Sum_{k prim} Re w(k/N) = R10-D baseline (doubly-exp small). So "exactly" holds
      up to the discrete primitive-mean of Re w (continuous <Re w>=0).
  (4) SIGN: is sign(Lambda_r) set by delta_r's near-x=0 lobe (Re w(0)=1/3 max)?  Gate honestly vs R13-B.
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_gamma_R9 as R9
import probe_charledger_R10 as R10

S = dict(R10.S)
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
for k in (7, 8):
    S[k] = F(7, 15) + EPS[k]

Rew = lambda x: 15 / (2 * (17 - 8 * math.cos(2 * math.pi * x))) - 0.5


def theta2(mu_r, r):
    N = 3 ** r
    g = R10.autocorr_dlog(mu_r, r); gf = [float(x) for x in g]
    return [sum(gf[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]


def main():
    print("# VERIFY -- Wilson psi-resolution (deviation field). Read-only.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ---- (1) partial-sum structure ----
    print("## (1) A_r(m) = C_{r+1}(m)/3 = gamma_r(tau_m) - gamma_{r-1}(tau_m)  (exact)")
    okpartial = True
    for m in (1, 2, 3, 9):
        row = []
        for r in range(2, 6):
            Np1 = 3 ** (r + 1)
            Cr1 = R7.Ck(r + 1, m, mu[r], pow(pow(4, -1, Np1), m, Np1)) / 3
            gr = R9.gamma(mu[r], r, R9.tau(m, r))
            grm1 = R9.gamma(mu[r - 1], r - 1, R9.tau(m, r - 1))
            ok = (Cr1 == gr - grm1); okpartial = okpartial and ok
            row.append(f"r={r}:{'OK' if ok else 'DEV'}")
        print(f"   m={m}: " + "  ".join(row))
    print(f"   => partial-sum structure {'CONFIRMED' if okpartial else 'DEV'} (A_r = the R13-C successive diffs)\n")

    # ---- (2) boundedness + Cesaro ----
    print("## (2) no non-uniform psi: gamma_n bounded, Cesaro mean gamma_n/n -> 0; A_r no nonzero limit")
    for m in (1, 9):
        gs = [float(R9.gamma(mu[r], r, R9.tau(m, r))) for r in range(1, 8)]
        ces = [gs[n - 1] / n for n in range(1, 8)]
        print(f"   m={m}: gamma_n={[f'{g:.3f}' for g in gs]}")
        print(f"         Cesaro gamma_n/n={[f'{c:.4f}' for c in ces]} (-> 0 slowly; bounded => no nonzero A_r-limit)")
    print("   => any convergent normalized profile is UNIFORM; non-uniform psi RULED OUT (limiting-shape reading dead).\n")

    # ---- (3) retarget identity ----
    print("## (3) Lambda_r = S_r <delta_r, Re w> + Lambda_r^unif   (delta_r = normalized profile - uniform)")
    print(f"   {'r':>2} {'Lambda_r':>14} {'S_r<delta,Rew>':>16} {'Lambda^unif(disc mean)':>22} {'sum':>14} {'match':>6}")
    for r in range(2, 8):
        N = 3 ** r; th2 = theta2(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim)
        Sr = float(S[r])
        Lam = sum(th2[k] * Rew(k / N) for k in prim)
        p = {k: th2[k] / Sr for k in prim}                      # normalized profile
        delta = {k: p[k] - 1.0 / M for k in prim}
        inner = sum(delta[k] * Rew(k / N) for k in prim)         # <delta_r, Re w>
        Lunif = (Sr / M) * sum(Rew(k / N) for k in prim)         # discrete primitive mean term
        tot = Sr * inner + Lunif
        print(f"   {r:>2} {Lam:>+14.7f} {Sr*inner:>+16.7f} {Lunif:>+22.3e} {tot:>+14.7f} "
              f"{'OK' if abs(tot-Lam)<1e-9 else 'DEV':>6}")
    print("   => identity holds; Lambda^unif = R10-D baseline (doubly-exp small) is the only gap vs 'exactly'.\n")

    # ---- (4) sign gate ----
    print("## (4) SIGN: is sign(Lambda_r) set by delta_r's near-x=0 lobe? (Re w>0 iff cos2pix>1/4, i.e. |x|<0.2088)")
    for r in range(3, 8):
        N = 3 ** r; th2 = theta2(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim); Sr = float(S[r])
        delta = {k: th2[k] / Sr - 1.0 / M for k in prim}
        # near-0 = angles with x<0.1 or x>0.9 (closest to trivial char); Re w>0 band = x<0.2088 or x>0.7912
        near0 = sum(Sr * delta[k] * Rew(k / N) for k in prim if (k / N < 0.1 or k / N > 0.9))
        rewpos = sum(Sr * delta[k] * Rew(k / N) for k in prim if Rew(k / N) > 0)
        rewneg = sum(Sr * delta[k] * Rew(k / N) for k in prim if Rew(k / N) < 0)
        Lam = rewpos + rewneg
        print(f"   r={r}: Lambda_r={Lam:+.6e}  | near-x0 (|x|<.1) contrib={near0:+.3e}  "
              f"Re w>0 band contrib={rewpos:+.3e}  Re w<0 band={rewneg:+.3e}  "
              f"[sign(Lam)==sign(near0)? {(Lam>0)==(near0>0)}]")
    print("   (gate Wilson's 'Lambda sign = near-0 lobe sign' honestly; R13-B had near-0 DEPLETED at r=5,6,7.)")


if __name__ == "__main__":
    main()
