"""
PROBE P6K -- MICROCOSM: the lambda-deformation (2026-07-26).

GUARDRAIL FOR THE WRITE-UP: at lambda != 1/2 the map is NOT Syracuse. This is a MODEL FAMILY; every result is about the
family's behaviour at and near the physical point lambda=1/2, not about Collatz directly.

Keep q=3, b=2, c=1. Deform the valuation law to P(v) prop lambda^v (v>=1), normalized; map still x -> (3x+1)/2^v.
lambda=1/2 => P(v)=2^-v = Syracuse. Off-critical => geometric convergence, T_i(lambda) converges in a few levels.
Drift E[3*2^-v] = 3(1-lambda)/(2-lambda) = 1 at lambda=1/2 (critical); <1 for lambda>1/2 (converge), >1 for lambda<1/2.

Construction: T_i(lambda) = 3^i Sum_{k>=1} 4^-k <rho_lambda, shift_k rho_lambda>, rho_lambda = base-4 numerator profile
from the lambda-stationary measure. (4^-k is the b=2 branch weight, hypothesised lambda-independent; M-A + M-D test it.)

M-A gate: lambda=1/2 reproduces the certified S-ladder to machine precision (T_i(0.5) == P6H T_i). If not, the
          deformation isn't passing through the real object.
M-B: lambda in {0.40,0.45,0.48,0.52,0.55,0.60}: T_i to convergence, lim, rate, full Lambda_i(lambda).
M-C: (1) Lambda_i(lambda) > 0 all i, all lambda? (2) does lim Lambda change sign in the family -> boundary in lambda?
     (3) as lambda->1/2 both sides, does lim T_i approach 0.238 continuously or with a KINK at the critical point?
M-D sanity: convergence rate should improve markedly off-critical (gap opening); compare to 3(1-lambda)/(2-lambda).
            If it doesn't speed up, the criticality id is wrong and M-C is meaningless.

stationary_trunc_lam = gated refactor of probe_27.stationary_trunc (weight lambda^v). No new transport.
"""
import os, sys, time
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
from probe_6_conservation_generalize import order_of_two
from probe_27_high_k_rho_q5 import stationary_trunc


def stationary_trunc_lam(q, k, lam, vmax=64):
    """Stationary pi with valuation law P(v) prop lam^v (v>=1). lam=0.5 == probe_27.stationary_trunc."""
    N = q ** k
    inv2 = pow(2, -1, N)
    M = order_of_two(N)
    vm = min(vmax, M)
    cp = np.array([r for r in range(N) if gcd(r, q) == 1], dtype=np.int64)
    n = len(cp)
    Z = sum(lam ** v for v in range(1, vm + 1))          # normalization Sum_{v=1}^{vm} lam^v
    inv_idx = np.full(N, -1, dtype=np.int64)
    inv_idx[cp] = np.arange(n)
    base_t = (q * cp + 1) % N
    rows_l, cols_l, vals_l = [], [], []
    inv2v = 1
    src = np.arange(n)
    for v in range(1, vm + 1):
        inv2v = (inv2v * inv2) % N
        t = (base_t * inv2v) % N
        rows_l.append(src); cols_l.append(inv_idx[t]); vals_l.append(np.full(n, (lam ** v) / Z))
    K = sp.csr_matrix((np.concatenate(vals_l), (np.concatenate(rows_l), np.concatenate(cols_l))), shape=(n, n))
    Kt = K.T.tocsr()
    pi = np.full(n, 1.0 / n)
    for _ in range(6000):
        nxt = Kt.dot(pi); s = nxt.sum()
        if s == 0:
            break
        nxt /= s
        if np.max(np.abs(nxt - pi)) < 1e-15:
            pi = nxt; break
        pi = nxt
    full = np.zeros(N); full[cp] = pi
    return full, cp


def rho4_lam(n, lam):
    q = 3 ** (n + 1); Nn = 3 ** n
    DL = np.full(q, -1, dtype=np.int64); g = 1
    for s in range(Nn):
        DL[g] = s; g = (g * 4) % q
    piW, _ = stationary_trunc_lam(3, n, lam)
    r = np.arange(Nn); cp = r[r % 3 != 0]
    nu = piW[cp]; nu = nu / nu.sum()
    Y = (3 * cp + 1) % q
    return np.bincount(DL[Y], weights=nu, minlength=Nn)


def autocorr(f):
    F = np.fft.fft(f); return np.fft.ifft(F * np.conj(F)).real


def T_lam(n, lam, KMAX=64):
    C = autocorr(rho4_lam(n, lam)); Nn = 3 ** n
    s = sum((4.0 ** -k) * C[k % Nn] for k in range(1, KMAX + 1))
    return 3 ** n * s


def main():
    t0 = time.time()
    print("# PROBE P6K -- MICROCOSM: lambda-deformation (MODEL FAMILY near lambda=1/2; NOT Collatz at lambda!=1/2)\n")

    # ---- M-A gate: lam=1/2 reproduces certified S-ladder ----
    print("## M-A gate: T_i(0.5) == certified S-ladder (P6H T_i)")
    P6H_T = {1: 0.23809524, 2: 0.23078734, 3: 0.23210720, 4: 0.23275746, 5: 0.23308438, 6: 0.23274571}
    for i in range(1, 7):
        t = T_lam(i, 0.5)
        print(f"   i={i}: T_i(0.5)={t:.8f}  P6H={P6H_T[i]:.8f}  diff={t-P6H_T[i]:.1e}")
    print()

    # ---- M-B / M-D: T_i(lambda) to convergence ----
    IMAX = 10
    lams = [0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60]
    print("## M-B/M-D: T_i(lambda) trajectory, rate, lim  (drift 3(1-l)/(2-l): <1 conv, =1 crit@0.5, >1 div)")
    lim = {}; sign_all_pos = {}; conv_rate = {}
    for lam in lams:
        T = {i: T_lam(i, lam) for i in range(1, IMAX + 1)}
        Lam = {i: T[i] - T[i - 1] for i in range(2, IMAX + 1)}
        drift = 3 * (1 - lam) / (2 - lam)
        # convergence rate = geometric ratio of |Lambda| over last levels (deparitied two-step)
        rates = [(abs(Lam[i]) / abs(Lam[i - 2])) ** 0.5 for i in range(8, IMAX + 1) if abs(Lam[i - 2]) > 1e-16]
        rate = np.mean(rates) if rates else float('nan')
        conv_rate[lam] = rate; lim[lam] = T[IMAX]
        sign_all_pos[lam] = all(Lam[i] > 0 for i in range(4, IMAX + 1))
        tail = "" if lam == 0.5 else (" DIVERGES" if drift > 1.0001 else "")
        print(f"   lam={lam:.2f} drift={drift:.4f}: T_10={T[IMAX]:.6f}  rate~{rate:.4f}  "
              f"Lam signs {'ALL+' if sign_all_pos[lam] else 'MIXED'}  Lam_10={Lam[IMAX]:+.2e}{tail}")
        print(f"        T: " + " ".join(f"{T[i]:.5f}" for i in range(1, IMAX + 1)))
    print()

    # ---- M-C: the three readings ----
    print("## M-C readings")
    print(f"   (1) Lambda_i > 0 for all i, all lambda?  " +
          ", ".join(f"{l}:{'+' if sign_all_pos[l] else 'MIX'}" for l in lams))
    print(f"   (2) sign of lim Lambda across family (does it flip -> boundary in lambda):")
    for lam in lams:
        # last Lambda sign as proxy for lim Lambda direction
        T9 = T_lam(9, lam); T10 = T_lam(10, lam)
        print(f"        lam={lam:.2f}: Lambda_10 = {T10-T9:+.3e}  ({'increasing' if T10>T9 else 'decreasing'})")
    print(f"   (3) lim T_i(lambda) vs lambda -- continuity / kink at 0.5 (critical value = limit of off-critical?):")
    for lam in lams:
        print(f"        lam={lam:.2f}: T_10={lim[lam]:.6f}  (2*T=S~{2*lim[lam]:.5f})")
    print(f"\n   M-D sanity: rate should DROP (converge faster) as |lambda-0.5| grows:")
    print("        " + " ".join(f"{l}:{conv_rate[l]:.3f}" for l in lams))
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
