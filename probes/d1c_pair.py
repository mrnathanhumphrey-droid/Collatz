"""
PROBE D1-C -- the tower's leading complex PAIR at L=4 (233k states, block SpMM subspace iteration).
Sanctioned instrument at L=4: block subspace (real top + complex pair + buffer), doubles, residuals <=1e-8.
No ARPACK, no shift-invert. Reports: partner (real), pair (modulus AND phase), within-block ratio |pair|/|partner|,
arg(pair). (Sequence context: within-block ratio 0.87 (L2) -> 0.97 (L3) -> ? ; arg across L.)

D1-D (deflated c0-side / overlap / g4): SKIPPED -- reported. At L=4 partner (0.33350) > c0 (~1/3), so the c0-side
right eigenvector is subdominant to the partner in the FULL operator AND its exact overlap/g4 need F2's convention;
deferred to avoid a mislabeled EP number. A-C are the spine (per spec).
"""
import numpy as np, scipy.sparse as sp, time, json, os
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

def build_tower(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    t0 = time.time()
    M, idx, n = build_M_gen(q, L, 2, raw)
    print(f"  build L={L}: n={n} nnz={M.nnz} {time.time()-t0:.1f}s", flush=True)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states])
    tw = np.where(gam != 0)[0]
    return M[tw][:, tw].tocsr().astype(np.float64), gam[tw], D

def block_subspace(Mt, b=6, tol=1e-8, maxit=1500, seed=7):
    rng = np.random.default_rng(seed)
    n = Mt.shape[0]
    V, _ = np.linalg.qr(rng.standard_normal((n, b)))
    hist = []
    for it in range(1, maxit + 1):
        W = Mt.dot(V)                      # SpMM, b columns
        H = V.T @ W                        # b x b Rayleigh quotient
        theta, Y = np.linalg.eig(H)
        oo = np.argsort(-np.abs(theta)); theta = theta[oo]; Y = Y[:, oo]
        X = V @ Y; MX = W @ Y              # Ritz vectors & their images
        res = np.array([np.linalg.norm(MX[:, k] - theta[k] * X[:, k]) for k in range(b)])
        V, _ = np.linalg.qr(W)
        if it % 20 == 0 or it < 3:
            print(f"    it {it}: top4 theta={[f'{t.real:.6f}{t.imag:+.6f}j' for t in theta[:4]]} "
                  f"res123={res[0]:.1e},{res[1]:.1e},{res[2]:.1e}", flush=True)
        hist.append((it, theta.copy(), res.copy()))
        if np.max(res[:3]) <= tol:
            print(f"    converged it {it}", flush=True)
            return theta, X, res, it
    return theta, X, res, maxit

def main():
    L = 4
    CACHE = os.path.expanduser("~/d1_tower")
    if os.path.exists(CACHE + "_Mt.npz"):
        Mt = sp.load_npz(CACHE + "_Mt.npz"); gam = np.load(CACHE + "_gam.npy")
        print(f"  cached tower {Mt.shape[0]} nnz={Mt.nnz}", flush=True)
    else:
        Mt, gam, D = build_tower(L)
        sp.save_npz(CACHE + "_Mt.npz", Mt); np.save(CACHE + "_gam.npy", gam)
    n = Mt.shape[0]
    print(f"  M_tower L=4: {n} states nnz={Mt.nnz}", flush=True)
    t0 = time.time()
    theta, X, res, it = block_subspace(Mt, b=6, tol=1e-8)
    print(f"  block subspace {time.time()-t0:.1f}s, {it} it", flush=True)
    # identify partner (top real) and leading complex pair
    isreal = np.abs(theta.imag) < 1e-7
    partner = float(theta[np.where(isreal)[0][0]].real)
    cpl = [t for t in theta if t.imag > 1e-7]
    pair = complex(sorted(cpl, key=lambda z: -abs(z))[0]) if cpl else None
    ratio = abs(pair) / partner if pair else None
    out = dict(L=4, partner=partner, pair=[pair.real, pair.imag] if pair else None,
               pair_mod=abs(pair) if pair else None, pair_arg=float(np.angle(pair)) if pair else None,
               within_ratio=ratio, top6=[[t.real, t.imag] for t in theta], res=[float(r) for r in res[:4]], it=it)
    print("RESULT_JSON_C " + json.dumps(out), flush=True)
    print(f"\n  === D1-C L=4 PAIR ===", flush=True)
    print(f"  partner (real) = {partner:.10f}", flush=True)
    if pair:
        print(f"  leading pair = {pair.real:.8f} {pair.imag:+.8f}j  |.|={abs(pair):.8f}  arg={np.angle(pair):.8f} rad "
              f"({np.degrees(np.angle(pair)):.4f} deg)", flush=True)
        print(f"  within-block ratio |pair|/partner = {ratio:.6f}   (seq 0.87 L2, 0.97 L3, {ratio:.4f} L4)", flush=True)
    print(f"  top-6 Ritz: {[f'{t.real:.6f}{t.imag:+.6f}j' for t in theta]}", flush=True)
    print(f"  D1-D SKIPPED (c0-side overlap/g4: partner>c0 at L4 + F2 convention needed; reported).", flush=True)
    with open("d1c_pair_L4.json", "w") as f:
        json.dump(out, f, indent=1)

if __name__ == "__main__":
    main()
