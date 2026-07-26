"""
PROBE CONE -- is the positive sign self-generated? A FALSIFIER for the flubber/self-generation frame.

K := { mu : Phi(mu) >= 0 },  Phi(mu) = C_q(1)/C_q(0),  q = P rho_mu  (fiber-fluctuation, kill 3|m modes).
Phi(nu_r) = d1_r > 0.  If the sign is manufactured by the transport, K is forward-invariant & attracting under T,
so adversarial starts with Phi<0 must ENTER K and stay. ONE surviving start (Phi<=0 for all 20 steps) kills the frame.

T = the CERTIFIED R16-A renewal (build_nu's step), reused, NOT reconstructed:
   on the a-index (X=1+3a, a in Z/3^r):  a -> a' = (2^{-v} (1+3a)) mod 3^r,  weight 2^{-v}, v~Geom(1/2).
   fixed point = nu_r (R8 self-similarity mu = sum_v p_v (f_v)_* mu).

P (kill 3|m frequency modes) in REAL space = subtract the (N/3)-periodic part (EXACT, no FFT):
   q(s) = rho(s) - (1/3)[rho(s)+rho(s+N/3)+rho(s+2N/3)].

Phi is a ratio of QUADRATICS in q => EVEN in q. So family-5 "flip the fluctuation sign" (rho_nu - 2 q_nu) and
"reflect" both give q -> +-q and Phi = +d1 UNCHANGED -- they start INSIDE K. Reported as the lead finding.

EXACT rational at r=4,5 (gate + constructions + first steps); float64 for r=4..8 trajectories (20 exact iterations
blow up denominators; sign-crossings -1 -> +0.002 are robust to 1e-12; float validated vs exact at r=5).

Non-convex K: NO convex-combination arguments. One survivor kills it; no success rates, no averaging.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact

D1_BANK = {4: 8.643329308003e-3, 5: 8.039124362558e-3, 6: 7.742271269236e-3,
           7: 5.786625297669e-3, 8: 4.665965059268e-3}
VMAX = 54


def dlog_array(r):
    d = R10.dlog_table(r); N = 3 ** r
    return np.array([d[a] for a in range(N)], dtype=np.int64)


# ---------------- float measures on the a-index ----------------
def mu_nu_float(r):
    N = 3 ** r; mu = np.zeros(N)
    for X, w in build_nu(0.5, r)[r].items():
        mu[(X - 1) // 3 % N] += float(w)
    return mu / mu.sum()

def rho_of(mu, dlog):
    rho = np.zeros_like(mu); rho[dlog] = mu       # rho[dlog[a]] = mu[a]
    return rho

def q_of(rho, N):
    n3 = N // 3
    fm = (rho + np.roll(rho, -n3) + np.roll(rho, -2 * n3)) / 3.0
    return rho - fm

def phi_of(mu, dlog, N):
    q = q_of(rho_of(mu, dlog), N)
    c0 = float(q @ q)
    if c0 <= 0:
        return float('nan'), 0.0
    return float(q @ np.roll(q, -1)) / c0, c0

def T_float(mu, r, inv2pows):
    N = 3 ** r; X = 1 + 3 * np.arange(N)
    new = np.zeros(N)
    for v in range(1, VMAX + 1):
        ar = (inv2pows[v] * X) % N
        np.add.at(new, ar, (2.0 ** -v) * mu)
    return new / new.sum()


# ---------------- exact (Fraction) for the gate + validation ----------------
def T_exact(mu, r, inv2pow_ex):
    N = 3 ** r; new = {}
    for a, w in mu.items():
        X = 1 + 3 * a
        for v in range(1, VMAX + 1):
            ar = (inv2pow_ex[v] * X) % N
            new[ar] = new.get(ar, F(0)) + F(1, 2 ** v) * w
    tot = sum(new.values())
    return {a: w / tot for a, w in new.items()}

def phi_exact(mu, r, d):
    N = 3 ** r; n3 = N // 3
    rho = [F(0)] * N
    for a, w in mu.items():
        rho[d[a]] += w
    fm = [(rho[s] + rho[(s + n3) % N] + rho[(s + 2 * n3) % N]) / 3 for s in range(N)]
    q = [rho[s] - fm[s] for s in range(N)]
    c0 = sum(x * x for x in q)
    if c0 == 0:
        return None
    c1 = sum(q[s] * q[(s + 1) % N] for s in range(N))
    return c1 / c0


def build_starts(r, dlog, mu_nu):
    """adversarial starts as float a-index measures; returns dict name->mu."""
    N = 3 ** r; starts = {}
    # F1 extremal high-freq mode (irrational; Phi ~ cos(2pi k/N) ~ -1)
    k = min((kk for kk in range(1, N) if kk % 3 != 0), key=lambda kk: np.cos(2 * np.pi * kk / N))
    rho_t = 1.0 + 0.99 * np.cos(2 * np.pi * k * np.arange(N) / N)
    m1 = rho_t[dlog].copy(); starts['F1_extremal'] = m1 / m1.sum()   # mu[a]=rho_t[dlog[a]]
    # F2 alternating dlog comb
    rho_t = (np.arange(N) % 2 == 0).astype(float)
    m2 = rho_t[dlog].copy(); starts['F2_altcomb'] = m2 / m2.sum()
    # F3 x4-anticorrelated support (greedy in X-space: 4S cap S = empty)
    inS = np.zeros(N, bool); taken4 = np.zeros(N, bool)
    for a in range(N):
        X = 1 + 3 * a; x4 = (4 * X) % (3 ** (r + 1)); a4 = (x4 - 1) // 3 % N
        if not taken4[a] and not inS[a4]:
            inS[a] = True; taken4[a4] = True
    m3 = np.zeros(N); m3[inS] = 1.0; starts['F3_x4anti'] = m3 / m3.sum()
    # F4 point masses
    m4a = np.zeros(N); m4a[N // 2 + 1] = 1.0; starts['F4a_pt_generic'] = m4a
    m4b = np.zeros(N); m4b[0] = 1.0; starts['F4b_pt_X1'] = m4b
    # F5 (should start INSIDE K, Phi=+d1, by evenness in q)
    rho_nu = rho_of(mu_nu, dlog)
    rho_ref = rho_nu[(-np.arange(N)) % N]
    m5a = rho_ref[dlog].copy(); starts['F5a_reflected'] = m5a / m5a.sum()
    qn = q_of(rho_nu, N); rho_flip = rho_nu - 2 * qn
    rho_flip = rho_flip - min(0.0, rho_flip.min())        # +const (kills nothing in q) to stay >=0
    m5b = rho_flip[dlog].copy(); starts['F5b_flip'] = m5b / m5b.sum()
    return starts


def main():
    t0 = time.time()
    print("# PROBE CONE -- is the positive sign self-generated? (falsifier)\n")

    # =============== CONE-A GATE (exact r=4,5 ; float r=4..8) ===============
    print("## CONE-A  GATE: Phi(nu_r)=d1_r, T[nu]=nu (fixed point), Phi(uniform): q==0")
    for r in (4, 5):
        N = 3 ** r; d = R10.dlog_table(r)
        nu = {(X - 1) // 3 % N: w for X, w in build_nu_exact(r)[r].items()}
        s = sum(nu.values()); nu = {a: w / s for a, w in nu.items()}
        inv2 = pow(2, -1, N); ip = [F(1)] + [F(pow(inv2, v, N)) for v in range(1, VMAX + 1)]
        phi_nu = float(phi_exact(nu, r, d))
        Tnu = T_exact(nu, r, ip)
        l1 = float(sum(abs(Tnu.get(a, F(0)) - nu.get(a, F(0))) for a in set(nu) | set(Tnu)))
        uni = {a: F(1, N) for a in range(N)}
        phu = phi_exact(uni, r, d)
        print(f"   r={r} (EXACT): Phi(nu)={phi_nu:+.12e} vs d1={D1_BANK[r]:+.12e} rel={abs(phi_nu-D1_BANK[r])/D1_BANK[r]:.1e}"
              f" | ||T nu - nu||1={l1:.1e} | Phi(uniform): {'q==0 OK' if phu is None else phu}")
    dlogs = {}; inv2f = {}
    for r in range(4, 9):
        N = 3 ** r; dlogs[r] = dlog_array(r)
        inv2 = pow(2, -1, N); inv2f[r] = [1] + [pow(inv2, v, N) for v in range(1, VMAX + 1)]
        mu_nu = mu_nu_float(r)
        phi_nu, _ = phi_of(mu_nu, dlogs[r], N)
        Tnu = T_float(mu_nu, r, inv2f[r]); l1 = float(np.abs(Tnu - mu_nu).sum())
        print(f"   r={r} (float): Phi(nu)={phi_nu:+.9e} vs d1={D1_BANK[r]:+.9e} rel={abs(phi_nu-D1_BANK[r])/D1_BANK[r]:.1e}"
              f" | ||T nu - nu||1={l1:.1e}")
    print()

    # exact-vs-float validation of the iteration at r=5 (family F4b, 4 steps)
    r = 5; N = 3 ** r; d = R10.dlog_table(r); inv2 = pow(2, -1, N)
    ip = [F(1)] + [F(pow(inv2, v, N)) for v in range(1, VMAX + 1)]
    me = {0: F(1)}; mf = np.zeros(N); mf[0] = 1.0
    print("## validation: exact vs float Phi trajectory, r=5, start=delta_{X=1}, 4 steps")
    for n in range(5):
        pe = phi_exact(me, r, d); pf = phi_of(mf, dlogs[5], N)[0]
        print(f"   step {n}: exact={float(pe) if pe is not None else 'q==0':>+.9e}  float={pf:+.9e}"
              if pe is not None else f"   step {n}: exact=q==0  float={pf:+.9e}")
        me = T_exact(me, r, ip); mf = T_float(mf, r, inv2f[5])
    print()

    # =============== CONE-B/C/D/E  float trajectories, r=4..8 ===============
    for r in (4, 6, 8):
        N = 3 ** r; dlog = dlogs[r]; mu_nu = mu_nu_float(r)
        _, c0_nu = phi_of(mu_nu, dlog, N)
        starts = build_starts(r, dlog, mu_nu)
        print(f"## CONE-B/C/D/E  r={r} (N={N})  [Phi(nu)=d1={D1_BANK[r]:+.4e}]")
        print(f"   {'family':<16}{'Phi_0':>11}{'entry':>7}{'stay+?':>7}{'settle':>11}{'L1@entry':>10}{'Phi@entry':>11}")
        for name, mu0 in starts.items():
            _, c0 = phi_of(mu0, dlog, N)
            if c0 < 1e-6 * c0_nu:
                print(f"   {name:<16}  excluded: C_q(0)={c0:.1e} < 1e-6*nu ({c0_nu:.1e}) -- q~0 ill-conditioned")
                continue
            mu = mu0.copy(); phis = []; l1s = []
            for n in range(21):
                ph, _ = phi_of(mu, dlog, N); phis.append(ph)
                l1s.append(float(np.abs(mu - mu_nu).sum()))
                mu = T_float(mu, r, inv2f[r])
            phis = np.array(phis)
            entry = next((n for n in range(len(phis)) if phis[n] > 0), None)
            stay = (entry is not None) and all(phis[n] > 0 for n in range(entry, len(phis)))
            settle = phis[-1]
            l1e = l1s[entry] if entry is not None else float('nan')
            phe = phis[entry] if entry is not None else float('nan')
            tag = "" if entry is not None else "  *** SURVIVOR (Phi<=0 all 20) -> FRAME DIES ***"
            print(f"   {name:<16}{phis[0]:>+11.4e}{('--' if entry is None else entry):>7}{str(stay):>7}"
                  f"{settle:>+11.4e}{l1e:>10.3f}{phe:>+11.4e}{tag}")
        print()
    print(f"# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
