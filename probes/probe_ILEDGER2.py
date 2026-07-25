"""
PROBE ILEDGER2 -- interference ledger, UN-NORMALIZED (mass-contracting) diagonal, extended to late j (G1-G4).

Corrections to ILEDGER: (1) diagonal is the UN-NORMALIZED |D|^2 kernel (mass-contracting; DC/avg weight
<|D|^2> = sum_v P(v)^2 = 1/3 = R17-B/C). Drop the /Z from push AND pullback (kept as an exact adjoint pair).
(2) delta_j via build_nu -> dlog pushforward -> |FFT|^2, so j can reach ~12 (mu-ladder walls at 8).

Gates:
  G0  cross-check: <delta_j, Re w> (FFT path) == known g_j for j<=7.  (else the profile path is wrong; STOP.)
  P4  telescoping gate sum_j A[j,r-j] = g_r (holds for any adjoint pair; sanity).
  G1  source rate rho_src = column rate A[j+1,k]/A[j,k] on the LATE window j=8..12. climbs 0.66 -> ~0.9?
  G2  closure: g_r ratios g_{r+1}/g_r vs exact Lambda ladder (0.950,0.897,0.913,0.893 at r=12..16) -- consistency
      of the ledger's overall rate (=max(rho_src,rho_prop)) with the data it must reproduce.
  G3  margin constant sign(c0 + c1 * rho_prop/(rho_src - rho_prop)) at the late rho_src; reconcile vs empirical.
  G4  rho_prop = row rate A[j,k+1]/A[j,k]: stable in j? (eigenvalue of a fixed kernel -- must be.)
"""
import os, sys, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

Rew = lambda x: 15.0 / (2 * (17 - 8 * np.cos(2 * np.pi * x))) - 0.5
JMAX = 12
LMAX = 12

_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): float(F(int(v['num']), int(v['den']))) for k, v in _hist.items()}
EPS.update({9: -7.520257156400000e-6, 10: 7.207509171100000e-4, 11: 1.501967012082273e-3, 12: 2.274713720558208e-3})
LAM_NU = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}
for r in range(12, 17):
    EPS[r + 1] = EPS[r] + 2 * LAM_NU[r]
LAM = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 17)}
Sval = {r: 7.0 / 15 + EPS[r] for r in range(1, 17)}


def push(s, L):                                  # UN-normalized forward T~_diag: level L -> L+1
    N = 3 ** L; Np = 3 * N
    kk = np.arange(Np)
    wD = 1.0 / (5 - 4 * np.cos(2 * np.pi * kk / Np))
    return wD * s[kk % N]                         # NO /Z


def delta_from_nu(nuj, j):
    """delta_j (dlog Plancherel profile deviation) from build_nu measure nu_j."""
    N = 3 ** j
    mu = np.zeros(N)
    for X, w in nuj.items():
        mu[(X - 1) // 3 % N] += float(w)          # a = (X-1)/3 mod 3^j
    d = R10.dlog_table(j)                          # beta(a)
    rho = np.zeros(N)
    np.add.at(rho, np.array([d[a] for a in range(N)]), mu)
    th = np.fft.fft(rho)
    prof = np.abs(th) ** 2                          # |theta_hat(k)|^2
    prim = np.array([k for k in range(1, N) if k % 3 != 0]); M = len(prim)
    S = float(prof[prim].sum())
    dd = np.zeros(N); dd[prim] = prof[prim] / S - 1.0 / M
    return dd, S


def main():
    t0 = time.time()
    print(f"# PROBE ILEDGER2 -- un-normalized (mass-contracting) diagonal, j to {JMAX}.\n")
    print("building nu (float) ...")
    nus = build_nu(0.5, JMAX)
    print(f"  nu built ({time.time()-t0:.1f}s); computing delta_j ...")
    dlt = {1: np.zeros(3)}; Sv = {1: 2.0 / 3}
    for j in range(2, JMAX + 1):
        dlt[j], Sv[j] = delta_from_nu(nus[j], j)
    print(f"  delta_j done ({time.time()-t0:.1f}s)\n")

    # ---- G0 cross-check ----
    print("## G0  cross-check <delta_j, Re w> (FFT path) vs known g_j = (Lambda_j - Lambda^unif)/S_j, j<=7")
    okG0 = True
    for j in range(2, 8):
        N = 3 ** j
        g_fft = float(np.sum(dlt[j] * Rew(np.arange(N) / N)))
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim)
        Lunif_over_S = (1.0 / M) * sum(Rew(k / N) for k in prim)
        g_known = LAM[j] / Sval[j] - Lunif_over_S
        ok = abs(g_fft - g_known) < 1e-7
        okG0 = okG0 and ok
        print(f"   j={j}: <delta,Rew>={g_fft:+.6e}  known g_j={g_known:+.6e}  [{'OK' if ok else 'DEV'}]")
    print(f"   => G0 {'PASS -- profile path correct' if okG0 else 'FAIL -- STOP'}\n")
    if not okG0:
        return

    # ---- sources + array ----
    s = {}
    for j in range(2, JMAX + 1):
        s[j] = dlt[j] - push(dlt[j - 1], j - 1)
    A = {j: {} for j in range(2, JMAX + 1)}
    for j in range(2, JMAX + 1):
        cur = s[j].copy(); L = j
        while L <= LMAX:
            N = 3 ** L
            A[j][L - j] = float(np.sum(cur * Rew(np.arange(N) / N)))
            if L < LMAX:
                cur = push(cur, L)
            L += 1

    # ---- P4 gate + fresh-source signs ----
    print("## P4  gate sum_j A[j,r-j] = g_r=<delta_r,Rew>  &  fresh source <s_j,Rew>=A[j,0] sign")
    for r in range(2, JMAX + 1):
        tot = sum(A[j][r - j] for j in range(2, r + 1))
        gr = float(np.sum(dlt[r] * Rew(np.arange(3 ** r) / 3 ** r)))
        fresh = A[r][0]
        print(f"   r={r:>2}: sum={tot:+.6e}  g_r={gr:+.6e}  [{'OK' if abs(tot-gr)<1e-7 else 'DEV'}]  "
              f"fresh<s_r,Rew>={fresh:+.3e} ({'+' if fresh>0 else '-'})")
    print()

    # ---- G4 row rate (rho_prop) stability ----
    print("## G4  rho_prop = row rate A[j,k+1]/A[j,k] (should be stable in j = eigenvalue of fixed kernel)")
    for j in range(2, JMAX + 1):
        rr = [A[j][k + 1] / A[j][k] for k in range(1, min(6, LMAX - j)) if abs(A[j].get(k, 0)) > 1e-15]
        print(f"   j={j:>2}: " + " ".join(f"{x:+.3f}" for x in rr))
    print()

    # ---- G1 column rate (rho_src) late window ----
    print("## G1  rho_src = column rate A[j+1,k]/A[j,k] vs j (does it climb 0.66 -> ~0.9 on late window?)")
    for k in (0, 1, 2):
        row = [(j, A[j + 1].get(k, float('nan')) / A[j].get(k, float('nan')))
               for j in range(2, JMAX) if abs(A[j].get(k, 0)) > 1e-15]
        print(f"   k={k}: " + " ".join(f"j{j}:{x:+.3f}" for j, x in row))
    print()

    # ---- G2 closure: g_r decay vs exact Lambda ladder ----
    print("## G2  closure: g_r ratio g_r/g_{r-1} (ledger, =Lambda_r ratio up to S) vs exact Lambda_r/Lambda_{r-1}")
    for r in range(4, JMAX + 1):
        gr = float(np.sum(dlt[r] * Rew(np.arange(3 ** r) / 3 ** r)))
        grm = float(np.sum(dlt[r - 1] * Rew(np.arange(3 ** (r - 1)) / 3 ** (r - 1))))
        lam_ratio = float(LAM[r] / LAM[r - 1])
        print(f"   r={r:>2}: g_r/g_(r-1)={gr/grm:+.4f}   Lambda_r/Lambda_(r-1)={lam_ratio:+.4f}")
    print("   [ledger overall rate = max(rho_src,rho_prop); must match Lambda's ~0.89-0.91 at large r.]\n")

    # ---- G3 margin constant ----
    print("## G3  margin constant  sign(c0 + c1 * rho_prop/(rho_src - rho_prop))")
    # estimate c0,c1 from the late, cleanest source (largest j with clean columns); rho_prop from G4 late
    jl = JMAX - 2                                        # need k=0,1,2 present -> j+2 <= LMAX
    c0 = A[jl][0]; c1 = A[jl][1]
    rp = A[jl][2] / A[jl][1] if abs(A[jl].get(1, 0)) > 1e-15 else float('nan')       # rho_prop late
    # rho_src late from k=1 column
    rs_pairs = [A[j + 1][1] / A[j][1] for j in range(JMAX - 4, JMAX - 1) if abs(A[j].get(1, 0)) > 1e-15]
    rs = float(np.mean(rs_pairs)) if rs_pairs else float('nan')
    tail = rp / (rs - rp) if (rs - rp) != 0 else float('nan')
    margin = c0 + c1 * tail
    print(f"   late j={jl}: c0=A[j,0]={c0:+.3e}, c1=A[j,1]={c1:+.3e}")
    print(f"   rho_prop(late)={rp:.3f}, rho_src(late,k=1)={rs:.3f}, tail=rho_prop/(rho_src-rho_prop)={tail:+.3f}")
    print(f"   margin = c0 + c1*tail = {margin:+.3e}  => sign {'NEG (7/15/rollover)' if margin<0 else 'POS (0.477)'}")
    print(f"   (c0<0 fresh source, c1>0 propagated; if rho_src climbs to ~0.9 the tail shrinks and NEG side gains.)")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
