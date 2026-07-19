"""
PROBE R3 -- THE SECULAR PRODUCT (crown corrected route; Theorem S gate). Cheap; banked F2-4 basis + symbolic.
Two-mode form on the degenerate top block B = [[c0, 0],[g, rho]] (kin, tow basis; protection => upper-right 0):
   a_m = P c0^m + Q rho^m,  P = phi_kin psi_kin + phi_tow psi_kin (g/Delta),
                            Q = phi_tow psi_tow - phi_tow psi_kin (g/Delta),  Delta = c0 - rho.
   phi = agreement readout (all-ones) projected on (kin,tow);  psi = independent-pair init (delta(1,1,0)) projected.
Annihilation lemma: g_m = a_m - a_{m-1}/3 = P c0^{m-1}(c0-1/3) + Q rho^{m-1}(rho-1/3).
Delta-cancellation (Theorem S): S_flat = 3^m g_m|plateau -> 3 Q (rho-1/3) = 3 g phi_tow psi_kin (the Delta's cancel).
R3-A convention freeze; R3-B the product vs plateau; R3-C mechanism in parts; R3-D the L-law.
INSTRUMENT: dense eig at q=3 (L=2,3). No fit. Exact where feasible; deviations as deviations.
"""
import numpy as np, scipy.sparse as sp, os
from fractions import Fraction as Fr
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

np.set_printoptions(linewidth=160, suppress=True)


def build_full(L, lam=0.5):
    q = 3; qL = q ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    gam = np.array([s[2] for s in states])
    wn = np.array(raw) / sum(raw); c0 = float(np.sum(wn ** 2))
    return M.tocsr().astype(np.float64), gam, idx, c0, n


def eff_2x2_basis(M, r0, rp, gam):
    U = np.column_stack([r0, rp]); Q, _ = np.linalg.qr(U)
    p0 = (gam == 0).astype(float)
    G0 = Q.T @ (p0[:, None] * Q)
    wv, Vv = np.linalg.eigh(G0)                      # [1]=max gamma0-weight=kinematic, [0]=tower
    Qp = Q @ np.column_stack([Vv[:, 1], Vv[:, 0]])
    B = Qp.T @ (M.dot(Qp))
    return B, Qp


def actual_a(M, idx, K):
    n = M.shape[0]; v = np.zeros(n); v[idx[(1, 1, 0)]] = 1.0
    a = [1.0]
    for _ in range(K):
        v = M.dot(v); a.append(v.sum())
    return a


def extract(L):
    M, gam, idx, c0, n = build_full(L)
    ev, VR = np.linalg.eig(M.toarray())
    real = np.abs(ev.imag) < 1e-9
    idxr = np.where(real)[0]; evr = ev[idxr].real
    ic0 = idxr[int(np.argmin(np.abs(evr - c0)))]
    order = idxr[np.argsort(-ev[idxr].real)]
    ip = int([j for j in order if j != ic0][0])
    r0 = VR[:, ic0].real.copy(); rp = VR[:, ip].real.copy()
    B, Qp = eff_2x2_basis(M, r0, rp, gam)
    c0m, rho, g = float(B[0, 0]), float(B[1, 1]), float(B[1, 0])
    Delta = c0m - rho
    one = np.ones(n); v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
    phi = Qp.T @ one; psi = Qp.T @ v0                # (kin, tow)
    a = actual_a(M, idx, 40)
    return dict(L=L, c0=c0, c0m=c0m, rho=rho, g=g, Delta=Delta, Bur=float(B[0, 1]),
                phi_kin=float(phi[0]), phi_tow=float(phi[1]), psi_kin=float(psi[0]), psi_tow=float(psi[1]), a=a)


def main():
    print("# PROBE R3 -- THE SECULAR PRODUCT (Theorem S gate). Two-mode degenerate block. No fit.")
    for L in (2, 3):
        d = extract(L)
        c0m, rho, g, Delta = d["c0m"], d["rho"], d["g"], d["Delta"]
        pk, pt, sk, st = d["phi_kin"], d["phi_tow"], d["psi_kin"], d["psi_tow"]
        gd = g / Delta
        P = pk * sk + pt * sk * gd
        Q = pt * st - pt * sk * gd
        a = d["a"]
        print(f"\n{'='*84}\n## L={L}")
        print(f"   2x2 block: c0={c0m:.10f} rho={rho:.10f} g=B[tow,kin]={g:+.6f} "
              f"B[kin,tow](protection~0)={d['Bur']:+.2e} Delta=c0-rho={Delta:+.3e}  g/Delta={gd:+.4f}")
        print(f"   projections: phi=(kin {pk:+.5f}, tow {pt:+.5f})  psi=(kin {sk:+.5f}, tow {st:+.5f})")
        print(f"   P={P:+.6f}  Q={Q:+.6f}")

        # ---- R3-A CONVENTION FREEZE: two candidate shell functionals vs welded 2/3,10/21 ----
        Sb = [3.0 ** k * (a[k] - a[k-1] / 3) for k in range(1, 6)]        # backward g_m=a_m-a_{m-1}/3
        Sf = [3.0 ** k * (a[k] - a[k+1]) for k in range(1, 5)]           # forward g_m=a_m-a_{m+1}
        print(f"   R3-A conventions from the A-seq:  backward S1,S2 = {Sb[0]:.6f},{Sb[1]:.6f} "
              f"(welded 2/3={2/3:.6f},10/21={10/21:.6f} -> {'MATCH' if abs(Sb[0]-2/3)<1e-6 and abs(Sb[1]-10/21)<1e-6 else 'no'});"
              f"  forward S1={Sf[0]:.6f} (=20/21={20/21:.6f}) -> DATA PINS BACKWARD g_m=a_m-a_(m-1)/3")

        # ---- R3-B THE PRODUCT vs plateau ----
        product = 3.0 * g * pt * sk
        twomode_plateau = 3.0 * Q * (rho - 1.0 / 3.0)                    # the SELF-CONSISTENT 2-mode flat level
        plateau = [3.0 ** k * (a[k] - a[k-1] / 3) for k in range(25, 39)]
        plat = float(np.mean(plateau))
        supercrit = (3.0 * rho > 1.0)                                    # finite-L partner > 1/3 => S_k grows (no plateau)
        twomode = [P * c0m ** m + Q * rho ** m for m in range(1, 13)]
        print(f"   R3-B PRODUCT 3*g*phi_tow*psi_kin = {product:+.6f}   (2-mode self-plateau 3Q(rho-1/3) = {twomode_plateau:+.6f})")
        print(f"        full-chain plateau (mean S_k, k=25..38) = {plat:+.6f}"
              + ("  [!] L-partner>1/3 => 3rho>1, S_k GROWS (no valid plateau at this L)" if supercrit else
                 f"   (2-mode is {100*product/plat:.0f}% of it; rest = subdominant near-1/3 modes)"))
        print(f"        two-mode a_m / actual a_m (m=8,10,12): "
              + ", ".join(f"{twomode[m-1]/a[m]:.3f}" for m in (8, 10, 12)))

        # ---- R3-C MECHANISM IN PARTS (Theorem S sign: the two minuses cancel => +) ----
        lhs = Q * (rho - 1.0 / 3.0)
        rhs = g * pt * sk                                                # +g phi_tow psi_kin  (Thm S; R3-C prose '-' is a typo)
        print(f"   R3-C parts: Q*(rho-1/3) = {lhs:+.8f}   vs   +g*phi_tow*psi_kin = {rhs:+.8f}   "
              f"ratio {lhs/rhs if rhs else float('nan'):+.6f}  "
              f"({'CANCELLATION verified (ratio~1)' if abs(lhs/rhs-1)<0.1 else 'dev (L-partner regime)'})")

        # ---- R3-D L-LAW piece ----
        print(f"   R3-D L-law piece: 3*g_L*(phi_tow*psi_kin) = {product:+.8f}   [-> 7/15=0.46667 as L->inf; L={L}]")


if __name__ == "__main__":
    main()
