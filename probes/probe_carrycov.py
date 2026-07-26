"""
PROBE CARRYCOV -- gate Wilson's carry-covariance conjecture (2026-07-26).

CONJECTURE:  dpi_{4x}  ~  T_{-c(x)} dpi_x     (the profile is CARRIED ALONG by the carry, not randomized).
If exact:  T_c dpi_{4x} = dpi_x  =>  q_r(1)-1/3 = E_mu ||dpi||^2  > 0   (a VARIANCE -- right SHAPE, not just size).

EXACT SPLIT (no assumption).  Write  dpi_{4x} = T_{-c} dpi_x + rho_x.  Then, with the lemma coherence
t = <dpi_x, T_c dpi_{4x}> (= u/v/w by c mod 3), one has identically
    q_r(1)-1/3 = E_mu||dpi_x||^2  +  E_mu<dpi_x, T_c rho_x>   ==   VAR + RES.
Cauchy-Schwarz:  RES >= -sqrt(VAR * E||rho||^2) = -VAR * D,   D := sqrt(E_mu||rho||^2)/sqrt(E_mu||dpi_x||^2).
=> SUFFICIENT CONDITION  D < 1  gives  q-1/3 >= VAR*(1-D) > 0.

Machinery reused verbatim from probe_carrylemma.gates_float: M=3^r; x=class mod M; c=floor(4x/M) in {0,1,2,3};
L=4x mod M; dpi_x = nu_hi(x+dM)/nu_low(x) - 1/3 (3-vector); mu(x) prop nu_low(x) nu_low(L). Rotations mod 3:
(T_a f)(d)=f(d+a) => T_a f = roll(f,-a); T_{-c} dpi_x = roll(dpi_x, c%3). rho = dpL - roll(dpx, c%3).

WILSON'S THREE READINGS (r<=16, top of ladder):
  (i)   D stably < 1                       -> sign established conditional on a MEASURED bound; pen target = prove D<1
                                              (an R8 self-similarity statement about dpi under x4, in the carry-native coord).
  (ii)  D ~ 1 and drifting UP              -> covariance is a low-r accident, dies like the rest.
  (iii) D < 1 but margin VAR*(1-D) (or the slack) TRACKS q-1/3 itself -> circular, no gain.

GATES:
  G0  VAR + RES == q_r(1)-1/3   (exact identity; confirms the split code) -- exact Fractions r=4,5, float r=4..16.
  G1  Cauchy-Schwarz slack: is VAR*(1-D) > 0 (rigorous-style lower bound holds) and how it compares to actual q-1/3.

PLUS the corrected RATE check (Wilson): DEPARITY first -- s_r=(e_r+e_{r+1})/2 kills the period-2 term -- then read
the geometric rate at the TOP (r=12..16), NOT the bottom. (Exact r<=8 can't see a change first visible at r=12-16.)
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact

RTOP = 16


def dense(nu, mod):
    a = np.zeros(mod)
    for X, w in nu.items():
        a[X % mod] += float(w)
    return a / a.sum()


def cov_gates_float(hi, r):
    """VAR, RES, lemma(=q-1/3), defect D, pointwise mean ratio, VAR*(1-D)."""
    M = 3 ** r
    low = hi[:M] + hi[M:2 * M] + hi[2 * M:3 * M]
    xs = np.nonzero(low > 0)[0]
    c = (4 * xs) // M
    L = (4 * xs) % M
    wgt = low[xs] * low[L]
    keep = wgt > 0
    xs, c, L, wgt = xs[keep], c[keep], L[keep], wgt[keep]
    Z = wgt.sum()
    cm = c % 3
    dpx = np.stack([hi[xs + d * M] / low[xs] for d in range(3)], 1) - 1.0 / 3     # (n,3)
    dpL = np.stack([hi[L + d * M] / low[L] for d in range(3)], 1) - 1.0 / 3
    # T_{-c} dpx = roll(dpx, cm) per row:  rolled[i,d] = dpx[i,(d-cm[i])%3]
    dd = np.arange(3)[None, :]
    idx_m = (dd - cm[:, None]) % 3                      # roll by +cm  (= T_{-c})
    Tm_dpx = np.take_along_axis(dpx, idx_m, axis=1)
    rho = dpL - Tm_dpx                                  # defect vector
    # lemma coherence t = <dpx, T_c dpL>,  T_c dpL = roll(dpL,-cm): rolled[i,d]=dpL[i,(d+cm)%3]
    idx_p = (dd + cm[:, None]) % 3
    Tc_dpL = np.take_along_axis(dpL, idx_p, axis=1)
    t = (dpx * Tc_dpL).sum(1)
    # RES = <dpx, T_c rho>, T_c rho = roll(rho,-cm)
    Tc_rho = np.take_along_axis(rho, idx_p, axis=1)
    res_i = (dpx * Tc_rho).sum(1)
    nsq_dpx = (dpx * dpx).sum(1)
    nsq_rho = (rho * rho).sum(1)
    VAR = float((wgt * nsq_dpx).sum() / Z)
    RES = float((wgt * res_i).sum() / Z)
    LEM = float((wgt * t).sum() / Z)                   # = q_r(1)-1/3
    E_nsq_rho = float((wgt * nsq_rho).sum() / Z)
    D = (E_nsq_rho / VAR) ** 0.5                        # mu-weighted RMS defect ratio
    good = nsq_dpx > 1e-30
    ptw = float((wgt[good] * np.sqrt(nsq_rho[good] / nsq_dpx[good])).sum() / wgt[good].sum())
    return dict(VAR=VAR, RES=RES, LEM=LEM, D=D, ptw=ptw, slack=VAR * (1 - D), Z=float(Z))


def cov_gates_exact(nex, r):
    M = 3 ** r
    hi = {}
    for X, w in nex[r].items():
        hi[X % (3 * M)] = hi.get(X % (3 * M), F(0)) + w
    tot = sum(hi.values()); hi = {k: w / tot for k, w in hi.items()}
    low = {}
    for X, w in hi.items():
        low[X % M] = low.get(X % M, F(0)) + w
    Z = F(0); VAR = F(0); RES = F(0); LEM = F(0)
    for x, lw in low.items():
        L = (4 * x) % M; cm = ((4 * x) // M) % 3
        lwL = low.get(L, F(0)); w0 = lw * lwL
        if w0 == 0:
            continue
        dpx = [hi.get(x + d * M, F(0)) / lw - F(1, 3) for d in range(3)]
        dpL = [hi.get(L + d * M, F(0)) / lwL - F(1, 3) for d in range(3)]
        Tm = [dpx[(d - cm) % 3] for d in range(3)]      # T_{-c} dpx
        rho = [dpL[d] - Tm[d] for d in range(3)]
        t = sum(dpx[d] * dpL[(d + cm) % 3] for d in range(3))         # <dpx, T_c dpL>
        res = sum(dpx[d] * rho[(d + cm) % 3] for d in range(3))       # <dpx, T_c rho>
        var = sum(dpx[d] * dpx[d] for d in range(3))
        Z += w0; VAR += w0 * var; RES += w0 * res; LEM += w0 * t
    return VAR / Z, RES / Z, LEM / Z


def main():
    t0 = time.time()
    print("# PROBE CARRYCOV -- carry-covariance conjecture: q-1/3 = VAR + RES, is defect D<1?\n")

    print("## EXACT (Fractions) r=4,5:  VAR+RES == q-1/3  (split identity)")
    nex = build_nu_exact(5)
    for r in (4, 5):
        VAR, RES, LEM = cov_gates_exact(nex, r)
        print(f"   r={r}: VAR={float(VAR):+.6e}  RES={float(RES):+.6e}  VAR+RES={float(VAR+RES):+.6e}"
              f"  q-1/3(lemma)={float(LEM):+.6e}  EQUAL: {VAR + RES == LEM}")
    print()

    print(f"## FLOAT r=4..{RTOP}   (build_nu to {RTOP} ... ~7 min)")
    nus = build_nu(0.5, RTOP)
    print(f"   built ({time.time()-t0:.1f}s)\n")

    rows = {}
    print(f"   {'r':>2} {'q-1/3':>11} {'VAR':>11} {'RES':>11} | {'D(RMS)':>7} {'ptw':>6} "
          f"{'VAR(1-D)':>11} {'(q-1/3)/VAR':>11}")
    for r in range(4, RTOP + 1):
        hi = dense(nus[r], 3 ** (r + 1))
        g = cov_gates_float(hi, r)
        del hi
        rows[r] = g
        gate = abs((g['VAR'] + g['RES']) - g['LEM']) / abs(g['LEM'])
        flag = " !split" if gate > 1e-6 else ""
        print(f"   {r:>2} {g['LEM']:>+11.4e} {g['VAR']:>+11.4e} {g['RES']:>+11.4e} | {g['D']:>7.4f} {g['ptw']:>6.3f} "
              f"{g['slack']:>+11.4e} {g['LEM']/g['VAR']:>11.4f}{flag}")
    print()

    # ---- readings ----
    Ds = [rows[r]['D'] for r in range(4, RTOP + 1)]
    print("## READINGS (Wilson's three)")
    print(f"   D sequence r=4..16: " + " ".join(f"{d:.3f}" for d in Ds))
    d_top = np.mean(Ds[-4:]); d_trend = Ds[-1] - Ds[-5] if len(Ds) >= 5 else float('nan')
    print(f"   D top-4 mean = {d_top:.4f};  D[16]-D[12] = {d_trend:+.4f}  "
          f"({'DRIFTING UP -> accident (ii)' if d_trend > 0.02 else 'stable/falling'})")
    below = all(d < 1 for d in Ds[-5:])
    print(f"   D<1 at all r>=12: {below}   -> "
          f"{'reading (i): sign conditional on measured bound D<1' if below and d_top < 1 else 'D not safely <1'}")
    # circularity: does slack VAR(1-D) track q-1/3? (ratio ~const across r => the smallness IS q-1/3, no gain)
    ratio_slack = [rows[r]['slack'] / rows[r]['LEM'] for r in range(12, RTOP + 1)]
    print(f"   slack/(q-1/3) r=12..16: " + " ".join(f"{x:+.3f}" for x in ratio_slack)
          + "   (if ~const & >0: C-S bound holds w/ room; if it TRACKS/~ -> circular (iii); if <0: C-S too lossy, sign via signed RES)")
    print()

    # ---- deparitied excess rate (corrected check) ----
    print("## DEPARITIED excess rate  (s_r=(e_r+e_{r+1})/2 kills period-2; read rate at TOP)")
    e = {r: rows[r]['LEM'] for r in range(4, RTOP + 1)}
    s = {r: 0.5 * (e[r] + e[r + 1]) for r in range(4, RTOP)}
    print("   raw excess two-step (|e_r/e_{r-2}|)^.5:  " +
          " ".join(f"r{r}:{(abs(e[r]/e[r-2]))**0.5:.3f}" for r in range(12, RTOP + 1)))
    print("   deparitied s_r successive ratio s_{r+1}/s_r: " +
          " ".join(f"r{r}:{s[r+1]/s[r]:.3f}" for r in range(11, RTOP - 1)))
    print("   [asymptotic excess rate -> tail sum -> S_inf. raw two-step ~0.89 was the cancellation tail;")
    print("    deparitied top ratio is the honest read of where the rate is actually heading.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
