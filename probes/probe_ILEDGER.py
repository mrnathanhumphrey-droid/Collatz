"""
PROBE ILEDGER -- the interference ledger (Wilson's reformulation of the coupling).

Second-moment decomposition: delta_r = second moment; the map on it is NOT the pushforward of the measure map.
Diagonal (|A|^2+|B|^2) = |D|^2 kernel (positivity, magnitudes); the CROSS term 2Re(A Bbar) = the SOURCE (bilinear,
phases) -- one inter-block interference object. Recursion delta_r = T~_diag delta_{r-1} + s_r, s_r = the cross term.

Projection = iterated pullback of the fixed weight (NOT an eigenvector -- the eigenvector doesn't exist):
  W^(0) = Re w;  W^(k+1)(x) = |D|^2-weighted AVERAGE of W^(k) over the 3 preimages {x/3,(x+1)/3,(x+2)/3},
  wD(y) = 1/(5-4cos 2pi y) = 1/|2 e(y)-1|^2  (=1 at x=0, =1/9 at x=1/2).
Exact decomposition (telescopes):  g_r = sum_{j=2}^r <s_j, W^(r-j)>,  computed as A[j,k]=<(T~_diag)^k s_j, Re w>
via the FORWARD adjoint of the pullback (exact adjoint pair on the dense primitive dlog grid).

P1 pullback family W^(k): positive-set (Haar) measure + value at x=0.
P2 sources s_j = delta_j - T~_diag delta_{j-1}, j=2..JMAX (exact data; s_j = subtraction, watch precision).
P3 array A[j,k].
P4 GATE: sum_j A[j,r-j] = g_r = <delta_r, Re w> for r<=7. FAIL => pullback convention wrong (#43 5th). STOP.
P5 two rates: row-rate in k (transport contraction), column-rate in j (source rate). Dominant direction's SIGN
   = the coupling.  Pre-reg verdict: dominant A[j,k] eventually single-signed NEG => rollover, 7/15 live;
   POS => 0.477; sign-indefinite with rates in noise => undecided, banked plainly.
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

Rew = lambda x: 15.0 / (2 * (17 - 8 * np.cos(2 * np.pi * x))) - 0.5
JMAX = 7          # sources s_j for j=2..JMAX (mu built to JMAX; the mu_8 build is the wall)
LMAX = 13         # push forward up to level 3^13 (1.6M)


def theta2(muj, j):
    N = 3 ** j
    g = R10.autocorr_dlog(muj, j); gf = np.array([float(x) for x in g])
    u = np.arange(N)
    return np.array([float(np.sum(gf * np.cos(2 * np.pi * k * u / N))) for k in range(N)])


def delta(muj, j):
    N = 3 ** j; th2 = theta2(muj, j)
    prim = np.array([k for k in range(1, N) if k % 3 != 0]); M = len(prim)
    S = float(th2[prim].sum())
    d = np.zeros(N); d[prim] = th2[prim] / S - 1.0 / M
    return d, S


def push(s, L):
    """forward T~_diag: level L (dense 3^L) -> level L+1 (dense 3^{L+1}); adjoint of the |D|^2 pullback."""
    N = 3 ** L; Np = 3 * N
    kk = np.arange(Np)
    wD = 1.0 / (5 - 4 * np.cos(2 * np.pi * kk / Np))
    Z = wD[0:N] + wD[N:2 * N] + wD[2 * N:3 * N]          # per-parent normalization
    par = kk % N
    return wD * s[par] / Z[par]


def pullback(W, L):
    """adjoint pullback: level L (dense 3^L) -> level L-1 (dense 3^{L-1}); |D|^2-weighted AVERAGE over 3 preimages."""
    N = 3 ** L; Nm = N // 3
    kk = np.arange(N)
    wD = 1.0 / (5 - 4 * np.cos(2 * np.pi * kk / N))
    x = np.arange(Nm)
    c0, c1, c2 = wD[x], wD[x + Nm], wD[x + 2 * Nm]        # children of x at level L: x, x+Nm, x+2Nm
    Z = c0 + c1 + c2
    return (c0 * W[x] + c1 * W[x + Nm] + c2 * W[x + 2 * Nm]) / Z


def main():
    print("# PROBE ILEDGER -- the interference ledger (second-moment source decomposition).\n")
    mu = {1: R7.mu1()}
    for k in range(2, JMAX + 1):
        mu[k] = R7.build_mu(mu[k - 1], k)
    dlt = {1: np.zeros(3)}; Sv = {1: float(R10.S[1])}
    for j in range(2, JMAX + 1):
        dlt[j], Sv[j] = delta(mu[j], j)

    # ---- P1: pullback family W^(k) ----
    print("## P1  PULLBACK FAMILY W^(k): positive-set (Haar) measure + value at x=0  (start Rew at level 12)")
    L0 = 12; N0 = 3 ** L0
    W = Rew(np.arange(N0) / N0); Ls = L0
    print(f"   {'k':>3} {'level':>6} {'Haar(W>0)':>10} {'W^k(0)':>12} {'max|W|':>10}")
    for k in range(0, 11):
        primmask = (np.arange(len(W)) % 3) != 0
        pos = float(np.mean(W[primmask] > 0))
        print(f"   {k:>3} {Ls:>6} {pos:>10.4f} {W[0]:>12.5e} {np.max(np.abs(W)):>10.4e}")
        if Ls > 1:
            W = pullback(W, Ls); Ls -= 1
    print("   [Haar(W>0)->fragments; W^k(0) = value at DC where |D|^2 concentrates & Re w is max +1/3 -- watch it.]\n")

    # ---- P2: sources ----
    print("## P2  SOURCES s_j = delta_j - T~_diag delta_{j-1}  (interference cross-term)")
    s = {}
    print(f"   {'j':>2} {'<s_j,Rew>':>13} {'||s_j||':>11} {'||s_j||/||d_j||':>15} {'sign':>5}")
    for j in range(2, JMAX + 1):
        tdj = push(dlt[j - 1], j - 1)                    # T~_diag delta_{j-1}, at level j
        s[j] = dlt[j] - tdj
        N = 3 ** j
        sr = float(np.sum(s[j] * Rew(np.arange(N) / N)))
        print(f"   {j:>2} {sr:>+13.5e} {np.linalg.norm(s[j]):>11.4e} {np.linalg.norm(s[j])/np.linalg.norm(dlt[j]):>15.4f} "
              f"{'+' if sr > 0 else '-':>5}")
    print("   [||s_j||/||d_j|| small => |D|^2 diagonal captures the bulk, s_j is the pure interference; large => not.]\n")

    # ---- P3 + P4: array A[j,k] and the GATE ----
    print("## P3/P4  A[j,k] = <(T~_diag)^k s_j, Re w>  and GATE  sum_j A[j,r-j] = g_r  (r<=7, exact target)")
    A = {j: {} for j in range(2, JMAX + 1)}
    for j in range(2, JMAX + 1):
        cur = s[j].copy(); L = j
        while L <= LMAX:
            N = 3 ** L
            A[j][L - j] = float(np.sum(cur * Rew(np.arange(N) / N)))
            if L < LMAX:
                cur = push(cur, L)
            L += 1
    print(f"   {'r':>2} {'sum_j A[j,r-j]':>16} {'g_r=<d_r,Rew>':>16} {'gate':>6}")
    okgate = True
    for r in range(2, JMAX + 1):
        tot = sum(A[j][r - j] for j in range(2, r + 1))
        gr = float(np.sum(dlt[r] * Rew(np.arange(3 ** r) / 3 ** r)))
        ok = abs(tot - gr) < 1e-9
        okgate = okgate and ok
        print(f"   {r:>2} {tot:>+16.9e} {gr:>+16.9e} {'OK' if ok else 'FAIL':>6}")
    print(f"   => GATE {'PASS -- forward/pullback adjoint pair consistent; decomposition valid' if okgate else 'FAIL -- pullback convention wrong (#43 5th). STOP, do not patch.'}\n")
    if not okgate:
        return

    # ---- P3 array print ----
    print("## P3  A[j,k] triangular array (rows j=2..%d, cols k=0..)" % JMAX)
    kmax = LMAX - 2
    print("   j\\k " + " ".join(f"{k:>10}" for k in range(0, min(9, kmax + 1))))
    for j in range(2, JMAX + 1):
        row = " ".join(f"{A[j].get(k, float('nan')):>+10.3e}" for k in range(0, min(9, kmax + 1)))
        print(f"   {j:>2}  {row}")
    print()

    # ---- P5: two rates + dominant-direction sign ----
    print("## P5  TWO RATES + dominant-direction sign")
    print("   row-rate A[j,k+1]/A[j,k] (transport contraction, per j):")
    for j in range(2, JMAX + 1):
        ks = [k for k in range(0, LMAX - j) if abs(A[j].get(k, 0)) > 1e-14 and abs(A[j].get(k + 1, 0)) > 1e-14]
        rr = [A[j][k + 1] / A[j][k] for k in ks[:8]]
        print(f"     j={j}: " + " ".join(f"{x:+.3f}" for x in rr))
    print("   column-rate A[j+1,k]/A[j,k] (source rate, per k):")
    for k in range(0, 6):
        cc = [A[j + 1].get(k, float('nan')) / A[j].get(k, float('nan')) for j in range(2, JMAX)
              if abs(A[j].get(k, 0)) > 1e-14]
        print(f"     k={k}: " + " ".join(f"{x:+.3f}" for x in cc))
    print("   diagonal g_r = sum_j A[j,r-j]; dominant term (largest |A[j,r-j]|) and its sign, per r:")
    for r in range(4, min(JMAX + 6, LMAX)):
        terms = [(j, A[j].get(r - j)) for j in range(2, min(r, JMAX) + 1) if (r - j) in A[j]]
        if not terms:
            continue
        jd, ad = max(terms, key=lambda t: abs(t[1]))
        tot = sum(a for _, a in terms)
        print(f"     r={r}: dominant j={jd} A={ad:+.3e} ({'+' if ad>0 else '-'}); partial sum_j(j<=min(r,{JMAX}))={tot:+.3e}")
    print("   [dominant-direction sign eventually NEG => rollover/7-15 live; POS => 0.477; indefinite => undecided.]")


if __name__ == "__main__":
    main()
