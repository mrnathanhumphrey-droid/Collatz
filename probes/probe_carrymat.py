"""
PROBE CARRYMAT -- Wilson's amended Option 1: the per-cell angular matrix, cone-positivity not PSD (2026-07-26).

The sign of q_r(1)-1/3 is the sign of cos(angle) between dpi_x and T_c dpi_{4x} -- an ANGLE-sign question, the one
class norm tools (C-S/Bochner/Poincare/D<1) are structurally blind to. Move it into the class that CAN see it:
a quadratic form. In the 2-D mean-zero plane of Z/3 (orthonormal basis e1=(1,-1,0)/sqrt2, e2=(1,1,-2)/sqrt6),
project a_x = P dpi_x, and per cell c (rotation cm=c%3, T_c = cyclic shift) form the mu-averaged 2x2 matrix

    S_c = sym( E_mu[ P(T_c dpi_{4x}) (x) P(dpi_x)^T | cell c ] ),   contribution U_c = p_c * tr(S_c)  (GATE: sum = q-1/3).

sym(T_c) = cos(2pi cm/3) I  =>  the ISOTROPIC drag is cos(120 cm): -1/2 in cells c=1,2 (rot 1,2). So PSD there needs
M far from I -- PROBABLY FALSE, and that is NOT a kill. We don't need PSD; we need positivity on the OBSERVED CONE
of a_x directions. So report, per cell (Wilson's amended 1/2/3):
  (1) eigenvalues + eigenvectors of S_c  (both, not trace/det);
  (2) distribution of observed directions phi_x = arg a_x under mu (weighted histogram);
  (3) overlap: observed direction-energy fraction on S_c's NEGATIVE eigendirection (f_neg; <0.5 = avoidance = mechanism).
Plus: isotropic tr(S_c)/2 vs cos(120 cm); the -1/2-vs-twist story for V1 (c=1).

ALSO (Wilson's cos-theta fix + a flag):
  kappa^2 = VAR_L/VAR, cos_theta_true = (1+kappa^2-D^2)/(2 kappa)  [geometric "1-D^2/2" = the kappa=1 shortcut];
  does kappa^2 cross 1 at the SAME r (~10) as U0+U3 crosses 0 (~12)?  "might be one event."

Machinery reused from probe_carrylemma/probe_fourcell (dense hi, cells, dpx/dpL, mu). One build_nu(0.5,16) ~7 min.
"""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_gapop_R28 import build_nu

RTOP = 16
SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"
# orthonormal basis of the mean-zero plane in R^3
E = np.array([[1, 1], [-1, 1], [0, -2]], float)
E[:, 0] /= np.sqrt(2.0); E[:, 1] /= np.sqrt(6.0)          # columns e1,e2 orthonormal


def dense(nu, mod):
    a = np.zeros(mod)
    for X, w in nu.items():
        a[X % mod] += float(w)
    return a / a.sum()


def cell_data(hi, r):
    """per x (nonzero low): cm, weight, projected a_x=P dpx, projected bC=P(T_c dpL)."""
    M = 3 ** r
    low = hi[:M] + hi[M:2 * M] + hi[2 * M:3 * M]
    xs = np.nonzero(low > 0)[0]
    c = (4 * xs) // M
    L = (4 * xs) % M
    wgt = low[xs] * low[L]
    keep = wgt > 0
    xs, c, L, wgt = xs[keep], c[keep], L[keep], wgt[keep]
    cm = c % 3
    dpx = np.stack([hi[xs + d * M] / low[xs] for d in range(3)], 1) - 1.0 / 3     # (n,3)
    dpL = np.stack([hi[L + d * M] / low[L] for d in range(3)], 1) - 1.0 / 3
    dd = np.arange(3)[None, :]
    idx = (dd + cm[:, None]) % 3                          # T_c dpL: (T_c f)(d)=f(d+cm) => roll(dpL,-cm)[d]=dpL[d+cm]
    TcdpL = np.take_along_axis(dpL, idx, axis=1)
    aX = dpx @ E                                          # (n,2) observed direction
    bC = TcdpL @ E                                        # (n,2) rotated-image direction; <aX,bC> = dpx.(T_c dpL)
    aL = dpL @ E                                          # for kappa^2 = VAR_L/VAR
    idxm = (dd - cm[:, None]) % 3                          # T_{-c} dpx = roll(dpx, cm): rolled[d]=dpx[(d-cm)%3]
    rho3 = dpL - np.take_along_axis(dpx, idxm, axis=1)     # exact CARRYCOV defect vector (3-vec)
    return c, cm, wgt, aX, bC, aL, rho3


def eig2(S):
    """eigen-decomp of symmetric 2x2, return (l_hi, l_lo, v_hi(angle deg), v_lo(angle deg))."""
    w, V = np.linalg.eigh(S)                              # ascending
    lo, hi = w[0], w[1]
    ang = lambda v: float(np.degrees(np.arctan2(v[1], v[0])))
    return hi, lo, ang(V[:, 1]), ang(V[:, 0])


def main():
    t0 = time.time()
    print("# PROBE CARRYMAT -- per-cell angular matrix S_c: cone-positivity, eigenvectors, direction histograms\n")
    print(f"building build_nu to {RTOP} ... (~7 min)")
    nus = build_nu(0.5, RTOP)
    print(f"  built ({time.time()-t0:.1f}s)\n")

    # ---- trend arrays ----
    trend = {}
    for r in range(4, RTOP + 1):
        hi = dense(nus[r], 3 ** (r + 1))
        c, cm, wgt, aX, bC, aL, rho3 = cell_data(hi, r)
        del hi
        Z = wgt.sum()
        VAR = float((wgt * (aX * aX).sum(1)).sum() / Z)
        VARL = float((wgt * (aL * aL).sum(1)).sum() / Z)
        Erho = float((wgt * (rho3 * rho3).sum(1)).sum() / Z)   # exact 3-vec defect (= CARRYCOV D)
        D = (Erho / VAR) ** 0.5
        kap2 = VARL / VAR
        cth = (1 + kap2 - D * D) / (2 * np.sqrt(kap2))
        # per-cell contributions U_c and matrices
        Uc = {}; Smat = {}
        for cell in range(4):
            sel = (c == cell)
            w = wgt[sel]; a = aX[sel]; b = bC[sel]
            G = (b * w[:, None]).T @ a / Z               # 2x2, sum_cell w b a^T / Z
            S = 0.5 * (G + G.T)
            Smat[cell] = S; Uc[cell] = float(np.trace(S))
        trend[r] = dict(VAR=VAR, kap2=kap2, D=D, cth=cth, Uc=Uc,
                        exc=sum(Uc.values()))
    print("## cos-theta FIX + kappa crossing  (Wilson: geometric 1-D^2/2 was the kappa=1 shortcut)")
    print(f"   {'r':>2} {'VAR':>10} {'kap2':>8} {'D':>7} {'1-D^2/2':>8} {'cos_th_true':>11} {'exc=q-1/3':>11} {'U0+U3':>11}")
    for r in range(4, RTOP + 1):
        t = trend[r]
        geo = 1 - t['D'] ** 2 / 2
        u03 = t['Uc'][0] + t['Uc'][3]
        print(f"   {r:>2} {t['VAR']:>10.4e} {t['kap2']:>8.4f} {t['D']:>7.4f} {geo:>8.4f} {t['cth']:>11.4f} "
              f"{t['exc']:>+11.4e} {u03:>+11.4e}")
    kcross = next((r for r in range(4, RTOP + 1) if trend[r]['kap2'] > 1), None)
    u03cross = next((r for r in range(5, RTOP + 1)
                     if (trend[r]['Uc'][0] + trend[r]['Uc'][3]) > 0 and
                        (trend[r - 1]['Uc'][0] + trend[r - 1]['Uc'][3]) <= 0), None)
    print(f"   => kappa^2 crosses 1 at r={kcross};  U0+U3 crosses 0 at r={u03cross}  "
          f"({'SAME window -- possibly one event' if kcross and u03cross and abs(kcross-u03cross) <= 2 else 'different windows'})")
    print()

    # ---- eigen-decomposition + histograms at r=16 ----
    r = RTOP
    hi = dense(nus[r], 3 ** (r + 1))
    c, cm, wgt, aX, bC, aL, _ = cell_data(hi, r)
    del hi
    Z = wgt.sum()
    print(f"## PER-CELL ANGULAR MATRIX S_c at r={r}  (eigenvalues+eigenvectors; contribution = p_c*tr = U_c)")
    print("   sym(T_c) isotropic drag = cos(120*cm): rot0=+1, rot1=rot2=-0.5\n")
    recon = 0.0
    NB = 12
    edges = np.linspace(-np.pi, np.pi, NB + 1)
    for cell in range(4):
        sel = (c == cell); rot = cell % 3
        w = wgt[sel]; a = aX[sel]; b = bC[sel]
        pc = float(w.sum() / Z)
        G = (b * w[:, None]).T @ a / w.sum()             # normalized (per unit mass): tr = U_c/p_c
        S = 0.5 * (G + G.T)
        recon += pc * float(np.trace(S))
        lh, ll, vh_ang, vl_ang = eig2(S)
        iso = float(np.trace(S)) / 2
        # observed direction covariance + its top axis
        Sig = (a * w[:, None]).T @ a / w.sum()
        sw, sV = np.linalg.eigh(Sig)
        obs_ang = float(np.degrees(np.arctan2(sV[1, 1], sV[0, 1])))
        # f_neg: observed energy fraction on S's negative eigendirection (v_lo)
        wS, VS = np.linalg.eigh(S)
        vneg = VS[:, 0]                                  # smaller eigenvalue's eigenvector
        f_neg = float((vneg @ Sig @ vneg) / np.trace(Sig))
        # phi histogram (weighted), folded to [0,180) since direction is a line (sign of a_x arbitrary? no -- keep full)
        phi = np.arctan2(a[:, 1], a[:, 0])
        h, _ = np.histogram(phi, bins=edges, weights=w)
        h = h / h.sum()
        bar = "".join("#" if v > 0.12 else (":" if v > 0.06 else ".") for v in h)
        print(f"   cell c={cell} (rot {rot}): p_c={pc:.4f}  U_c=p_c*tr={pc*np.trace(S):+.4e}  iso={iso:+.4f} "
              f"(drag {np.cos(2*np.pi*rot/3):+.2f})")
        print(f"       eig S_c: lambda+ = {lh:+.4f} @ {vh_ang:+6.1f}deg ,  lambda- = {ll:+.4f} @ {vl_ang:+6.1f}deg   "
              f"(anisotropy {(lh-ll)/2:.4f})")
        print(f"       observed dir: top axis @ {obs_ang:+6.1f}deg ,  f_neg(mass on lambda- dir) = {f_neg:.3f} "
              f"({'AVOIDS neg dir (mechanism)' if f_neg < 0.45 else ('on neg dir!' if f_neg > 0.55 else 'diffuse')})")
        print(f"       arg a_x hist [-pi,pi) x{NB}: {bar}")
    print(f"\n   GATE: sum p_c tr(S_c) = {recon:+.6e}  vs banked q_16(1)-1/3 = +4.1789e-04  "
          f"rel {abs(recon-4.1789e-4)/4.1789e-4:.1e}")

    # save trend
    json.dump({str(r): {'VAR': trend[r]['VAR'], 'kap2': trend[r]['kap2'], 'D': trend[r]['D'],
                        'cth': trend[r]['cth'], 'exc': trend[r]['exc'],
                        'Uc': {str(k): v for k, v in trend[r]['Uc'].items()}} for r in trend},
              open(os.path.join(SCRATCH, 'carrymat.json'), 'w'))
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
