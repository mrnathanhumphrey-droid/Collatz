"""
W4c -- validate the DEFLATION method against DENSE at L=3 (do we trust the L=4 overlap?).
At L=3: (a) dense eig -> true r_p, r0, overlap_dense (=0.99999 banked); (b) SAME deflation pipeline as L=4
(tower r_p/l_p embedded -> deflate -> c0 r0) -> overlap_defl. If they agree, L=4's 0.034 is real; if
deflation gives ~0.03 while dense gives 0.99999, the method can't measure a near-1 overlap (ill-conditioned).
"""
import numpy as np, scipy.sparse as sp, time, json
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

def build_full(L, lam=0.5):
    q = 3; qL = q ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    M, idx, n = build_M_gen(q, L, 2, [lam ** d for d in range(1, D + 1)])
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states])
    wn = np.array([lam ** d for d in range(1, D + 1)]); wn /= wn.sum(); c0 = float(np.sum(wn ** 2))
    return M.tocsr().astype(np.float64), gam, c0

def pit(A, tol=1e-11, maxit=20000, seed=1):
    n = A.shape[0]; rng = np.random.default_rng(seed)
    v = np.abs(rng.standard_normal(n)); v /= np.linalg.norm(v); rp = 0.0
    for it in range(1, maxit + 1):
        w = A.dot(v); rho = float(v @ w); res = float(np.linalg.norm(w - rho * v))
        v = w / np.linalg.norm(w)
        if res <= tol and abs(rho - rp) < 1e-15: break
        rp = rho
    return rho, v, it, res

def ov(a, b): return abs(float(a @ b)) / (np.linalg.norm(a) * np.linalg.norm(b))

def main():
    L = 3
    M, gam, c0 = build_full(L)
    n = M.shape[0]; tw = np.where(gam != 0)[0]
    print(f"L={L} full {n} tower {len(tw)} c0={c0:.10f}", flush=True)
    # (a) DENSE truth
    ev, VR = np.linalg.eig(M.toarray())
    real = np.abs(ev.imag) < 1e-9; idxr = np.where(real)[0]; evr = ev[idxr].real
    ic0 = idxr[int(np.argmin(np.abs(evr - c0)))]
    order = idxr[np.argsort(-ev[idxr].real)]; ip = int([j for j in order if j != ic0][0])
    r0_d = VR[:, ic0].real.copy(); rp_d = VR[:, ip].real.copy()
    lam_p = float(ev[ip].real)
    ov_dense = ov(r0_d, rp_d)
    print(f"  DENSE: partner={lam_p:.10f} c0_mode={ev[ic0].real:.10f} overlap_dense={ov_dense:.8f}", flush=True)
    # how much of the dense partner eigvec lives on gamma=0?
    p0 = (gam == 0); frac_p_g0 = float(np.sum(rp_d[p0] ** 2) / np.sum(rp_d ** 2))
    frac_0_g0 = float(np.sum(r0_d[p0] ** 2) / np.sum(r0_d ** 2))
    print(f"  dense partner gamma=0 mass frac={frac_p_g0:.4e}; c0-mode gamma=0 frac={frac_0_g0:.4e}", flush=True)

    # (b) DEFLATION pipeline (same as L=4)
    Mt = M[tw][:, tw].tocsr()
    lp_r, rv, _, _ = pit(Mt, seed=1); lp_l, lv, _, _ = pit(Mt.T.tocsr(), seed=2)
    r_p = np.zeros(n); r_p[tw] = rv; r_p /= np.linalg.norm(r_p)
    l_p = np.zeros(n); l_p[tw] = lv; l_p /= np.linalg.norm(l_p)
    denom = float(l_p @ r_p); lamp2 = 0.5 * (lp_r + lp_l)
    print(f"  DEFL tower partner={lamp2:.10f} l_p.r_p={denom:.4e}", flush=True)
    # overlap of the tower-embedded partner vs DENSE partner (are they the same vector?)
    print(f"  <r_p_tower, r_p_dense> = {ov(r_p, rp_d):.8f}   (1.0 => embedded partner == dense partner)", flush=True)
    class Defl:
        shape = (n, n)
        def dot(self, v): return M.dot(v) - lamp2 * r_p * (float(l_p @ v) / denom)
    c0v, r0_defl, it3, r3 = pit(Defl(), tol=1e-10, maxit=30000, seed=5)
    print(f"  DEFL c0={c0v:.10f} ({it3}it res{r3:.1e})", flush=True)
    ov_defl_towerp = ov(r0_defl, r_p)          # overlap using tower-embedded partner (as L=4 did)
    ov_defl_densep = ov(r0_defl, rp_d)         # overlap using dense partner
    ov_r0 = ov(r0_defl, r0_d)                  # is the deflated c0 vec == dense c0 vec?
    print(f"  overlap(r0_defl, r_p_tower) = {ov_defl_towerp:.8f}   [what L=4 computed]", flush=True)
    print(f"  overlap(r0_defl, r_p_dense) = {ov_defl_densep:.8f}", flush=True)
    print(f"  <r0_defl, r0_dense> = {ov_r0:.8f}   (1.0 => deflation found the true c0 eigenvector)", flush=True)
    print(json.dumps(dict(L=3, ov_dense=ov_dense, ov_defl_towerp=ov_defl_towerp,
                          ov_defl_densep=ov_defl_densep, r0_match=ov_r0,
                          rp_match=ov(r_p, rp_d), frac_p_g0=frac_p_g0, frac_0_g0=frac_0_g0)), flush=True)

if __name__ == "__main__":
    main()
