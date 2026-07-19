"""
R3-D L=4 point: the secular product 3*g*phi_tow*psi_kin at L=4 from the banked w4_full cache.
Block-2 subspace iteration (SpMV+QR, instrument-legal) -> converged c0-mode r0 + partner rp, then eff_2x2 +
project agreement readout (all-ones) and independent-pair init (delta(1,1,0)=index 0). No fit.
"""
import numpy as np, scipy.sparse as sp, os, time, json

CACHE = os.path.expanduser("~/w4_full")
PARTNER = 0.33349990132218854


def block2_split(M, b=2, tol=1e-8, maxit=3000, seed=3):
    n = M.shape[0]; rng = np.random.default_rng(seed)
    V, _ = np.linalg.qr(rng.standard_normal((n, b)))
    t0 = time.time()
    for it in range(1, maxit + 1):
        W = M.dot(V); H = V.T @ W
        theta, Y = np.linalg.eig(H)
        X = V @ Y; MX = W @ Y
        res = np.array([np.linalg.norm(MX[:, k] - theta[k] * X[:, k]) for k in range(b)])
        V, _ = np.linalg.qr(W)
        if it % 100 == 0:
            print(f"    it {it}: theta={[f'{t.real:.9f}' for t in theta]} res={[f'{r:.1e}' for r in res]} ({time.time()-t0:.0f}s)", flush=True)
        if np.max(res) <= tol:
            print(f"    converged it {it} ({time.time()-t0:.0f}s)", flush=True)
            break
    return theta, X, res, it


def eff_2x2_basis(M, r0, rp, gam):
    U = np.column_stack([r0, rp]); Q, _ = np.linalg.qr(U)
    p0 = (gam == 0).astype(float)
    G0 = Q.T @ (p0[:, None] * Q)
    wv, Vv = np.linalg.eigh(G0)
    Qp = Q @ np.column_stack([Vv[:, 1], Vv[:, 0]])
    B = Qp.T @ (M.dot(Qp))
    return B, Qp


def build_full_L4():
    from probe_phase2a_q2b_q6 import build_M_gen, subgroup
    q = 3; L = 4; qL = q ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [0.5 ** d for d in range(1, D + 1)]
    t0 = time.time()
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    gam = np.array([s[2] for s in states])
    M = M.tocsr().astype(np.float64)
    print(f"  built L=4: n={n} nnz={M.nnz} ({time.time()-t0:.0f}s)", flush=True)
    sp.save_npz(CACHE + "_M.npz", M); np.save(CACHE + "_gam.npy", gam)
    return M, gam


def main():
    if os.path.exists(CACHE + "_M.npz"):
        M = sp.load_npz(CACHE + "_M.npz").tocsr(); gam = np.load(CACHE + "_gam.npy")
    else:
        M, gam = build_full_L4()
    wn = np.array([0.5 ** d for d in range(1, 55)]); wn /= wn.sum(); c0 = float(np.sum(wn ** 2))
    print(f"# R3-D L=4: full M {M.shape[0]} nnz={M.nnz}, c0={c0:.12f}", flush=True)
    theta, X, res, it = block2_split(M)
    reals = theta.real
    ic0 = int(np.argmin(np.abs(reals - c0))); ip = 1 - ic0
    c0_mode = float(theta[ic0].real); partner = float(theta[ip].real)
    print(f"  split: c0_mode={c0_mode:.10f} (c0={c0:.10f}) partner={partner:.10f} (target {PARTNER:.10f}) maxres={np.max(res):.1e}", flush=True)
    r0 = X[:, ic0].real.copy(); rp = X[:, ip].real.copy()
    B, Qp = eff_2x2_basis(M, r0, rp, gam)
    c0m, rho, g = float(B[0, 0]), float(B[1, 1]), float(B[1, 0])
    Delta = c0m - rho
    n = M.shape[0]; one = np.ones(n); v0 = np.zeros(n); v0[0] = 1.0     # (1,1,0) = index 0
    phi = Qp.T @ one; psi = Qp.T @ v0
    pk, pt, sk, st = float(phi[0]), float(phi[1]), float(psi[0]), float(psi[1])
    gd = g / Delta
    Q = pt * st - pt * sk * gd
    product = 3.0 * g * pt * sk
    lhs = Q * (rho - 1.0 / 3.0); rhs = g * pt * sk
    out = dict(L=4, c0=c0, c0m=c0m, rho=rho, g=g, Delta=Delta, g_over_Delta=gd,
               phi_kin=pk, phi_tow=pt, psi_kin=sk, psi_tow=st, Q=Q, product=product,
               R3C_lhs=lhs, R3C_rhs=rhs, R3C_ratio=lhs / rhs, maxres=float(np.max(res)), it=it)
    print("RESULT_JSON " + json.dumps(out), flush=True)
    print(f"\n  === R3 L=4 ===\n  c0={c0m:.10f} rho={rho:.10f} g={g:+.6f} Delta={Delta:+.3e} g/Delta={gd:+.2f} "
          f"(protection B[kin,tow]={B[0,1]:+.1e})", flush=True)
    print(f"  phi=(kin {pk:+.5f}, tow {pt:+.5f})  psi=(kin {sk:+.5f}, tow {st:+.5f})", flush=True)
    print(f"  R3-C: Q*(rho-1/3)={lhs:+.8f} vs +g*phi_tow*psi_kin={rhs:+.8f} ratio {lhs/rhs:+.6f}", flush=True)
    print(f"  R3-D PRODUCT 3*g*phi_tow*psi_kin = {product:+.8f}   [L-law: 0.3066(L2), 0.3512(L3), {product:.4f}(L4) -> 7/15=0.46667]", flush=True)
    with open("outputs/r3_L4.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
