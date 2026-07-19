"""
PROBE J2 RIDER -- the 12-digit L=4 doublet SPLIT (block-splitting ledger's judge).
Sanctioned within-block subspace iteration (block Rayleigh-Ritz) on the cached L=4 tower -- NO ARPACK,
NO shift-invert (INSTRUMENT LAW; the q=3 operator is defective near the EP). Real block (M is real);
Ritz values from the projected nonsymmetric H = Q^T M Q carry the conjugate doublet pair.
Report the two doublet members + splitting, to the precision ACTUALLY achieved (near the defective EP the
matvec floor may bite ~1e-10; no fabricated digits). Banked target: split L=3 2.644e-3 -> L=4 ~2.0e-4 (x0.076).
"""
import numpy as np, scipy.sparse as sp, os, time, json

CACHE = os.path.expanduser("~/j2_L4")


def subspace_iter(Mt, block=28, maxit=400, seed=7, tol=1e-13):
    n = Mt.shape[0]
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((n, block)))
    prev = None
    hist = []
    t0 = time.time()
    for it in range(1, maxit + 1):
        Z = Mt.dot(Q)
        Q, _ = np.linalg.qr(Z)
        H = Q.T @ (Mt.dot(Q))                                    # block x block projected operator
        ev = np.linalg.eig(H)[0]
        # the doublet = the complex conj pair of LARGEST modulus (partner is real & larger)
        cpx = sorted([z for z in ev if abs(z.imag) > 1e-9], key=lambda z: -abs(z))
        top = cpx[0] if cpx else 0
        hist.append(top)
        if it % 20 == 0 or it < 3:
            realmax = max((z.real for z in ev if abs(z.imag) < 1e-9), default=float('nan'))
            print(f"    it {it}: partner~{realmax:.12f}  doublet-top {top.real:.12f}{top.imag:+.12f}j  "
                  f"({time.time()-t0:.0f}s)", flush=True)
        if prev is not None and abs(top - prev) < tol:
            print(f"    converged at it {it} (d_top {abs(top-prev):.1e})", flush=True)
            break
        prev = top
    return ev, hist, it


def main():
    print("# PROBE J2 RIDER -- 12-digit L=4 doublet split. Within-block subspace iteration (instrument-legal).")
    Mt = sp.load_npz(CACHE + "_Mt.npz").tocsr()
    print(f"  tower {Mt.shape[0]} nnz {Mt.nnz}", flush=True)
    ev, hist, it = subspace_iter(Mt)
    cpx = sorted([z for z in ev if abs(z.imag) > 1e-9 and z.imag > 0], key=lambda z: -abs(z))
    reals = sorted([z.real for z in ev if abs(z.imag) < 1e-9], reverse=True)
    partner = reals[0] if reals else float('nan')
    print(f"\n  === RIDER: L=4 doublet split ===", flush=True)
    print(f"  partner (real, largest) = {partner:.12f}", flush=True)
    print(f"  top complex pairs (upper half), by modulus:", flush=True)
    for z in cpx[:4]:
        print(f"     {z.real:.12f} {z.imag:+.12f}j   |.|={abs(z):.12f}  arg={np.angle(z):.12f}", flush=True)
    if len(cpx) >= 2:
        m0, m1 = cpx[0], cpx[1]
        split = abs(m0 - m1)
        print(f"\n  doublet members: {m0.real:.12f}{m0.imag:+.12f}j  &  {m1.real:.12f}{m1.imag:+.12f}j", flush=True)
        print(f"  SPLITTING |p0-p1| = {split:.6e}", flush=True)
        print(f"  banked: L=3 2.644285e-3 -> L=4 (block-6) 2.002e-4 (ratio 0.0757);  this split = {split:.4e}  "
              f"ratio-to-L3 = {split/2.644285e-3:.4f}", flush=True)
        out = dict(partner=partner, m0=[m0.real, m0.imag], m1=[m1.real, m1.imag], split=split,
                   ratio_to_L3=split / 2.644285e-3, iters=it,
                   doublet_top_history=[[z.real, z.imag] for z in hist[-8:]])
        with open("outputs/judge_L4_doublet_split.json", "w") as f:
            json.dump(out, f, indent=1)
        print("RESULT_JSON " + json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
