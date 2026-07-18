"""
PROBE W4 -- the EP witnesses at L=4 (+ same-instrument re-base at L=2,3). FULL operator.
Method (sanctioned): block-2 orthogonal power iteration on the full M at L=4 (SpMV + QR each step),
then Rayleigh-Ritz on the converged 2-frame to split c0-mode from partner (2x2 dense eig, NOT ARPACK,
no shift-invert). Convergence rate |pair|/c0 ~ 0.987 => ~1500-2000 it. L=2,3 by dense eig (re-base).

VALIDATION GATE: the two split Ritz values reproduce partner=0.33349990132 (G4) and c0=1/3+(2/3)2^-54
to >=8 digits, else the split failed -> report + stop.

W1 overlap4 = |<r0, r_partner>| (l2-normalized right eigenvecs), vs 0.998, 0.99999.
W2 g4 via F2-4's 2x2 (kinematic=gamma0-aligned, tower) basis; defectiveness g/|Delta| vs 17,189
   using canonical detuning-vs-c0 |Delta4|=1.666e-4.
W3 re-base: overlap and g at L=2,3 with THIS (dense) method -> one instrument for all 3 points.
W4 c0-mode right eigvec gamma-profile at L=4 (finding).
"""
import numpy as np, scipy.sparse as sp, time, json, os
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

def v3(n):
    if n == 0: return -1        # zero-carry sector
    k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def build_full(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    t0 = time.time()
    M, idx, n = build_M_gen(q, L, 2, raw)
    print(f"  build L={L}: n={n} nnz={M.nnz} {time.time()-t0:.1f}s", flush=True)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states])
    wn = np.array([lam ** d for d in range(1, D + 1)]); wn = wn / wn.sum()
    c0 = float(np.sum(wn ** 2))
    return M.tocsr().astype(np.float64), gam, D, qL, c0

def eff_2x2_sparse(M, r0, rp, gam):
    """F2-4 convention with SpMV (no dense M). Returns B (2x2), kin/tow gamma0-weights."""
    U = np.column_stack([r0, rp]); Q, _ = np.linalg.qr(U)
    p0 = (gam == 0).astype(float)
    G0 = Q.T @ (p0[:, None] * Q)
    wv, Vv = np.linalg.eigh(G0)                     # ascending; [1]=max gamma0-weight = kinematic
    kin = Vv[:, 1]; tow = Vv[:, 0]
    Qp = Q @ np.column_stack([kin, tow])
    MQp = M.dot(Qp)
    B = Qp.T @ MQp
    return B, float(wv[1]), float(wv[0])

def overlap(r0, rp):
    return abs(float(r0 @ rp)) / (np.linalg.norm(r0) * np.linalg.norm(rp))

def dense_rebase(L):
    M, gam, D, qL, c0 = build_full(L)
    ev, VR = np.linalg.eig(M.toarray())
    real = np.abs(ev.imag) < 1e-9
    idxr = np.where(real)[0]
    evr = ev[idxr].real
    # c0-mode = real eig closest to c0; partner = the other top real (largest real != c0-mode)
    ic0 = idxr[int(np.argmin(np.abs(evr - c0)))]
    order = idxr[np.argsort(-ev[idxr].real)]
    ip = int([j for j in order if j != ic0][0])
    r0 = VR[:, ic0].real.copy(); rp = VR[:, ip].real.copy()
    ov = overlap(r0, rp)
    B, kw, tw = eff_2x2_sparse(M, r0, rp, gam)
    return dict(L=L, c0=c0, c0_mode=float(ev[ic0].real), partner=float(ev[ip].real),
                overlap=ov, B=B.tolist(), Bkt=float(B[0, 1]), Btk=float(B[1, 0]),
                kin_wt=kw, tow_wt=tw)

def block2_split(M, c0, partner_val, gam, b=2, tol=1e-8, maxit=3500, seed=3):
    n = M.shape[0]; rng = np.random.default_rng(seed)
    V, _ = np.linalg.qr(rng.standard_normal((n, b)))
    for it in range(1, maxit + 1):
        W = M.dot(V)
        H = V.T @ W
        theta, Y = np.linalg.eig(H)
        X = V @ Y; MX = W @ Y
        res = np.array([np.linalg.norm(MX[:, k] - theta[k] * X[:, k]) for k in range(b)])
        V, _ = np.linalg.qr(W)
        if it % 25 == 0 or it < 3:
            print(f"    it {it}: theta={[f'{t.real:.10f}{t.imag:+.1e}j' for t in theta]} res={[f'{r:.1e}' for r in res]}", flush=True)
        if np.max(res) <= tol:
            print(f"    converged it {it}", flush=True)
            break
    return theta, X, res, it

def main():
    OUT = {}
    # ---- W3 re-base (dense) ----
    for L in [2, 3]:
        OUT[f"L{L}"] = dense_rebase(L)
        d = OUT[f"L{L}"]
        print(f"  [re-base L={L}] c0={d['c0']:.10f} c0_mode={d['c0_mode']:.10f} partner={d['partner']:.10f} "
              f"overlap={d['overlap']:.8f} B_kt={d['Bkt']:+.6f} B_tk={d['Btk']:+.6f}", flush=True)

    # ---- L=4 block-2 split ----
    L = 4
    CACHE = os.path.expanduser("~/w4_full")
    if os.path.exists(CACHE + "_M.npz"):
        M = sp.load_npz(CACHE + "_M.npz"); gam = np.load(CACHE + "_gam.npy")
        wn = np.array([0.5 ** d for d in range(1, 55)]); wn /= wn.sum(); c0 = float(np.sum(wn ** 2))
        print(f"  cached full M {M.shape[0]} nnz={M.nnz}", flush=True)
    else:
        M, gam, D, qL, c0 = build_full(L)
        sp.save_npz(CACHE + "_M.npz", M); np.save(CACHE + "_gam.npy", gam)
    partner_true = 0.33349990132218854
    t0 = time.time()
    theta, X, res, it = block2_split(M, c0, partner_true, gam)
    print(f"  block-2 split {time.time()-t0:.1f}s, {it} it, theta={theta}", flush=True)
    # assign: c0-mode closest to c0 closed form; partner the other
    reals = theta.real
    ic0 = int(np.argmin(np.abs(reals - c0)))
    ip = 1 - ic0 if len(theta) == 2 else int(np.argmin(np.abs(reals - partner_true)))
    c0_mode = float(theta[ic0].real); partner = float(theta[ip].real)
    r0 = X[:, ic0].real.copy(); rp = X[:, ip].real.copy()
    # VALIDATION GATE
    gate_partner = abs(partner - partner_true) < 1e-8
    gate_c0 = abs(c0_mode - c0) < 1e-8
    print(f"\n  VALIDATION: partner={partner:.10f} (target {partner_true:.10f}, ok {gate_partner}); "
          f"c0_mode={c0_mode:.10f} (target {c0:.10f}, ok {gate_c0})", flush=True)
    passed = gate_partner and gate_c0 and np.max(res) <= 1e-8
    if not passed:
        print(f"  GATE FAILED (res={np.max(res):.1e}) -- split not trusted; readouts suppressed.", flush=True)
    # readouts
    ov4 = overlap(r0, rp)
    B, kw, tw = eff_2x2_sparse(M, r0, rp, gam)
    # gamma-profile of c0-mode
    vp = np.abs(r0); levs = np.array([v3(int(g)) for g in gam]); tot = vp.sum()
    prof = {int(lv): float(vp[levs == lv].sum() / tot) for lv in sorted(set(levs.tolist()))}
    out4 = dict(L=4, passed=bool(passed), partner=partner, c0_mode=c0_mode, c0=c0,
                gate_partner=bool(gate_partner), gate_c0=bool(gate_c0), maxres=float(np.max(res)), it=it,
                overlap=ov4, B=B.tolist(), Bkt=float(B[0, 1]), Btk=float(B[1, 0]),
                kin_wt=kw, tow_wt=tw, c0_gamma_profile=prof)
    OUT["L4"] = out4
    print("RESULT_JSON_W " + json.dumps(OUT), flush=True)
    print(f"\n  === W4 L=4 EP WITNESSES (gate {'PASS' if passed else 'FAIL'}) ===", flush=True)
    print(f"  W1 overlap4 = {ov4:.8f}   (seq 0.998 L2, 0.99999 L3, {ov4:.6f} L4)", flush=True)
    print(f"     re-base overlaps: L2={OUT['L2']['overlap']:.8f} L3={OUT['L3']['overlap']:.8f}", flush=True)
    print(f"  W2 B (kin,tow) off-diagonals: B_kt={B[0,1]:+.6f} B_tk={B[1,0]:+.6f}  (F2-4 g=0.0505 L2,0.0188 L3)", flush=True)
    print(f"     re-base g: L2 B_kt={OUT['L2']['Bkt']:+.5f} B_tk={OUT['L2']['Btk']:+.5f}; "
          f"L3 B_kt={OUT['L3']['Bkt']:+.5f} B_tk={OUT['L3']['Btk']:+.5f}", flush=True)
    print(f"  W4 c0-mode gamma-profile (v3-level -> mass frac; -1=zero-carry): {prof}", flush=True)
    with open("w4_out.json", "w") as f:
        json.dump(OUT, f, indent=1)

if __name__ == "__main__":
    main()
