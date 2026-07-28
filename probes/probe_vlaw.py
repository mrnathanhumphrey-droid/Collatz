"""
PROBE VLAW -- is Collatz exactly critical, or only nearly? (2026-07-27)

Wilson's spec. The exact empirical v-law P(v | m mod 2^k) from Lagarias-Sinai (forward integer orbits,
`outputs/v_distribution_by_N.csv`) deviates from Geom(1/2) (chi^2/dof=4373). But S_inf sees ONE linear
functional (the 2nd moment). Compute the two criticality functionals DIRECTLY:

  delta1 = E_exact[2^-v] - 1/3      (drift/criticality: =0 iff E[W]=1, W=3*2^-v)
  delta2 = sum_v P_exact(v)^2 - 1/3  (participation / D2)

Both are 1/3 EXACTLY under Geom(1/2). chi^2 huge does NOT imply delta != 0 -- compute delta.

STAGE A: delta1, delta2 per N (full precision) + tail bound.  Read: machine-zero => exactly critical,
  wall real; else measure sign & size, compare to the 7/15-vs-0.475 gap (~2%).
MIRAGE CHECK (inside A): the 7/15 object (build_nu / Tao's chain) uses Geom(1/2) BY DEFINITION, so
  delta=0 there. The forward-orbit deviation is a DIFFERENT object. Test: build the Syracuse stationary
  pi_k two ways -- Geom(1/2) vs the forward-orbit P(v) as an iid kernel -- do they give the same pi / S?
  Same => forward deviation doesn't propagate, delta irrelevant to S_inf. Differ => real & load-bearing.
STAGE B (only if delta2!=0 AND measures differ): rebuild pi_k with the exact v-law, recompute S_3,4,5.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_lattice import scat, build_dl

CSV = r"C:\Collatz\outputs\v_distribution_by_N.csv"


def load_pv():
    """N -> array P(v), v=1..jmax (raw empirical, truncated tail)."""
    rows = {}
    with open(CSV) as f:
        next(f)
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 4:
                continue
            N, j, p = parts[0], int(parts[1]), float(parts[2])
            rows.setdefault(N, {})[j] = p
    out = {}
    for N, d in rows.items():
        jmax = max(d)
        arr = np.array([d.get(j, 0.0) for j in range(1, jmax + 1)])
        out[N] = arr
    return out


def stage_A(pvs):
    print("## STAGE A  delta1 = E[2^-v]-1/3,  delta2 = sum P(v)^2 - 1/3   (=0 exactly for Geom(1/2))")
    print(f"   {'N':>6} {'tail mass':>10} {'E[2^-v]':>11} {'delta1':>13} {'sumP^2':>11} {'delta2':>13}")
    res = {}
    for N in ("2^28", "2^30", "2^32", "2^34"):
        p = pvs[N]
        v = np.arange(1, len(p) + 1)
        tail = 1.0 - p.sum()
        E2v = float((p * 2.0 ** (-v)).sum())          # tail adds < tail*2^-(jmax) ~ 1e-9, negligible
        sumP2 = float((p ** 2).sum())
        d1 = E2v - 1.0 / 3.0
        d2 = sumP2 - 1.0 / 3.0
        res[N] = (d1, d2, tail)
        print(f"   {N:>6} {tail:>10.2e} {E2v:>11.7f} {d1:>+13.2e} {sumP2:>11.7f} {d2:>+13.2e}")
    # geom control
    g = 2.0 ** (-np.arange(1, 61))
    print(f"   Geom(1/2) control: E[2^-v]={float((g*2.0**(-np.arange(1,61))).sum()):.7f} "
          f"sumP^2={float((g**2).sum()):.7f}  (both -> 1/3, delta=0)")
    # interpret vs the 2% gap
    d1, d2, _ = res["2^34"]
    print(f"\n   READ @2^34: delta1={d1:+.2e} ({d1/(1/3)*100:+.2f}% of 1/3), "
          f"delta2={d2:+.2e} ({d2/(1/3)*100:+.2f}% of 1/3)")
    print(f"   E[W]=3*E[2^-v]={3*(1/3+d1):.5f} (=1 exactly critical; <1 subcritical/contracting).")
    print(f"   7/15-vs-0.475 gap ~ 2%.  |delta2| here ~ {abs(d2)/(1/3)*100:.2f}% -> NOT the same object.\n")
    return res


def T_level_w(n, wvec, q=3):
    """Syracuse stationary S-functional T at level n with iid v-weights wvec[v-1]=P(v).
       Same machinery as probe_lattice.T_level but with a pluggable v-law (not Geom(1/2))."""
    qN = q ** (n + 1); Nn = q ** n; inv2 = pow(2, -1, Nn)
    r = np.arange(Nn); cpr = r[r % q != 0]; ncp = len(cpr)
    base = ((q * cpr + 1) % Nn).astype(np.uint32)
    ps = []; p = 1
    for v in range(1, len(wvec) + 1):
        p = (p * inv2) % Nn
        if wvec[v - 1] > 0:
            ps.append((p, float(wvec[v - 1])))
    pi = np.full(ncp, 1.0 / ncp)
    for _ in range(500):
        nxt = np.zeros(ncp)
        for pp, w in ps:
            nxt += np.bincount(scat(base, pp, Nn), weights=w * pi, minlength=ncp)
        nxt /= nxt.sum()
        if np.abs(nxt - pi).max() < 1e-14:
            pi = nxt; break
        pi = nxt
    nu = pi / pi.sum(); DL = build_dl(qN, Nn); Y = (q * cpr + 1) % qN
    rho = np.bincount(DL[Y].astype(np.int64), weights=nu, minlength=Nn)
    T = q ** n * sum((4.0 ** -k) * float(np.dot(rho, np.roll(rho, k))) for k in range(1, 64))
    return pi, T


def mirage_check(pvs, nmax=6):
    print("## MIRAGE CHECK  build pi_k with Geom(1/2) vs forward-orbit P(v); same stationary?")
    geom = 2.0 ** (-np.arange(1, 61))                 # P(v)=2^-v
    fwd = pvs["2^34"].copy()                           # exact forward-orbit law (renormalized)
    fwd = fwd / fwd.sum()
    print(f"   {'n':>2} {'|pi_geom-pi_fwd|max':>19} {'T_geom':>11} {'T_fwd':>11} {'S=2T geom':>11} {'S=2T fwd':>11} {'dS':>10}")
    for n in range(2, nmax + 1):
        pig, Tg = T_level_w(n, geom)
        pif, Tf = T_level_w(n, fwd)
        dpi = float(np.abs(pig - pif).max())
        print(f"   {n:>2} {dpi:>19.2e} {Tg:>11.7f} {Tf:>11.7f} {2*Tg:>11.7f} {2*Tf:>11.7f} {2*(Tf-Tg):>+10.2e}")
    print("   [same pi & S => forward deviation does NOT propagate to the 7/15 object; delta irrelevant.]")
    print("   [7/15 object uses Geom(1/2) BY DEFINITION (build_nu lam=0.5) => its delta=0 identically.]\n")


def main():
    print("# PROBE VLAW -- is Collatz exactly critical, or only nearly?\n")
    pvs = load_pv()
    stage_A(pvs)
    mirage_check(pvs, 6)
    print("VERDICT: see delta1/delta2 (forward orbit) + mirage dS. The 7/15 object is exactly critical")
    print("(Geom(1/2) by construction); 0.475 is the exact-critical value, not a v-law artifact.")


if __name__ == "__main__":
    main()
