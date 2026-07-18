"""
W4b -- L=4 c0-mode via OBLIQUE DEFLATION (block-2 orthogonal split failed: partner & c0 right eigenvectors
are numerically parallel at overlap->1, so an orthogonal frame collapses onto the complex pair, not c0).
Instrument-legal: get partner right r_p (power iter on M_tower) + left l_p (power iter on M_tower^T), embed
to full (0 on gamma=0), deflate M1 = M - lam_p r_p l_p^T/(l_p^T r_p), power-iterate M1 -> c0 dominant.
No shift-invert, no ARPACK. VALIDATION GATE: recovered c0 eigenvalue == 1/3+(2/3)2^-54 to >=6 digits.
Then overlap4 = |<r0,r_p>|, g4 via F2-4 2x2, c0 gamma-profile.
"""
import numpy as np, scipy.sparse as sp, time, json, os
LADDER = None

def v3(n):
    if n == 0: return -1
    k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def pit(A, tol=1e-7, maxit=4000, seed=1, tag=""):
    n = A.shape[0]; rng = np.random.default_rng(seed)
    v = np.abs(rng.standard_normal(n)); v /= np.linalg.norm(v)
    rp = 0.0
    for it in range(1, maxit + 1):
        w = A.dot(v); rho = float(v @ w); res = float(np.linalg.norm(w - rho * v))
        v = w / np.linalg.norm(w)
        if it % 100 == 0: print(f"    {tag} it{it} rho={rho:.10f} res={res:.1e}", flush=True)
        if res <= tol and abs(rho - rp) < 1e-13: break
        rp = rho
    return rho, v, it, res

def eff_2x2(M, r0, rp, gam):
    U = np.column_stack([r0, rp]); Q, _ = np.linalg.qr(U)
    p0 = (gam == 0).astype(float)
    G0 = Q.T @ (p0[:, None] * Q)
    wv, Vv = np.linalg.eigh(G0); kin = Vv[:, 1]; tow = Vv[:, 0]
    Qp = Q @ np.column_stack([kin, tow]); MQp = M.dot(Qp); B = Qp.T @ MQp
    return B, float(wv[1]), float(wv[0])

def main():
    M = sp.load_npz(os.path.expanduser("~/w4_full_M.npz")).tocsr()
    gam = np.load(os.path.expanduser("~/w4_full_gam.npy"))
    wn = np.array([0.5 ** d for d in range(1, 55)]); wn /= wn.sum(); c0 = float(np.sum(wn ** 2))
    n = M.shape[0]; tw = np.where(gam != 0)[0]
    print(f"  full M {n} nnz={M.nnz}; tower {len(tw)}; c0={c0:.12f}", flush=True)
    Mt = M[tw][:, tw].tocsr(); MtT = Mt.T.tocsr()
    lam_p_true = 0.33349990132218854

    t0 = time.time()
    lp_r, rv, it1, rr = pit(Mt, tol=1e-8, maxit=4000, seed=1, tag="r_p")      # partner right (tower)
    lp_l, lv, it2, rl = pit(MtT, tol=1e-8, maxit=4000, seed=2, tag="l_p")     # partner left (tower)
    print(f"  partner tower: rho_right={lp_r:.10f} ({it1}it), rho_left={lp_l:.10f} ({it2}it) {time.time()-t0:.0f}s", flush=True)
    lam_p = 0.5 * (lp_r + lp_l)
    # embed to full
    r_p = np.zeros(n); r_p[tw] = rv; r_p /= np.linalg.norm(r_p)
    l_p = np.zeros(n); l_p[tw] = lv; l_p /= np.linalg.norm(l_p)
    denom = float(l_p @ r_p)
    print(f"  lam_p={lam_p:.10f} (vs G4 {lam_p_true:.10f}); l_p.r_p={denom:.6e}", flush=True)

    class Defl:
        shape = (n, n)
        def dot(self, v):
            return M.dot(v) - lam_p * r_p * (float(l_p @ v) / denom)
    M1 = Defl()
    t0 = time.time()
    c0_val, r0, it3, r3 = pit(M1, tol=1e-8, maxit=4000, seed=5, tag="c0")
    print(f"  deflated c0: {c0_val:.12f} ({it3}it res{r3:.1e}) {time.time()-t0:.0f}s", flush=True)

    gate_c0 = abs(c0_val - c0) < 1e-6
    print(f"\n  VALIDATION: c0_recovered={c0_val:.10f} vs closed-form {c0:.10f}  |diff|={abs(c0_val-c0):.2e}  "
          f"gate {'PASS' if gate_c0 else 'FAIL'}", flush=True)

    r0 = r0 / np.linalg.norm(r0)
    ov4 = abs(float(r0 @ r_p)) / (np.linalg.norm(r0) * np.linalg.norm(r_p))
    B, kw, tw_wt = eff_2x2(M, r0, r_p, gam)
    vp = np.abs(r0); levs = np.array([v3(int(g)) for g in gam]); tot = vp.sum()
    prof = {int(lv): float(vp[levs == lv].sum() / tot) for lv in sorted(set(levs.tolist()))}
    out = dict(L=4, method="oblique_deflation", gate_c0=bool(gate_c0), c0_recovered=c0_val, c0=c0,
               lam_p=lam_p, overlap4=ov4, B=B.tolist(), Bkt=float(B[0, 1]), Btk=float(B[1, 0]),
               kin_wt=kw, tow_wt=tw_wt, c0_gamma_profile=prof, deflate_it=it3, deflate_res=r3)
    print("RESULT_JSON_W4B " + json.dumps(out), flush=True)
    print(f"\n  === W4b L=4 (deflation) ===", flush=True)
    print(f"  W1 overlap4 = {ov4:.8f}   (0.99834 L2, 0.99999 L3, {ov4:.6f} L4)", flush=True)
    print(f"  W2 g4 = B_tk = {B[1,0]:+.6f}  (0.05053 L2, 0.01882 L3); B_kt={B[0,1]:+.2e}", flush=True)
    dlt = 1.666e-4
    print(f"     defectiveness g4/|Delta4| = {abs(B[1,0])/dlt:.1f}  (17 L2, 189 L3)", flush=True)
    print(f"  W4 c0-mode gamma-profile: {prof}", flush=True)
    with open("w4b_out.json", "w") as f: json.dump(out, f, indent=1)

if __name__ == "__main__":
    main()
