"""
PROBE G4 -- THE PARTNER AT L=4 (pure CPU numpy/scipy; power iteration WITHIN the tower block).
Partner = rho(M_tower) = dominant eigenvalue of the gamma!=0 principal submatrix (G0). Within-block
runner-up ~0.29 (ratio ~0.87) => power iteration stable & sanctioned (instrument law; no ARPACK/shift-invert).

R1 rho(M_tower,4) >=8 sig digits, two independent starts, plateau agreement.
R2 Delta_4 = 1/3 - rho_4 (sign+magnitude). c0(4)=1/3+(2/3)2^-54, corr 3.7e-17 < double eps => 1/3 is ref.
   Sequence: Delta = -2.911e-3 (L2), +9.958e-5 (L3), ?   (braid point 3 = the SIGN, reported, no pre-pick)
R3 ON-RECORD SHOT (zero weight until DERIVED): 27^-L predicts |Delta_4| ~ 2*27^-4 = 3.76e-6. measured vs pred.
R4 gamma-level mass profile of converged partner eigenvector. PRE-REG shape: geometric cascade+truncation
   [2/3,2/9,2/27,1/27,0] vs L=3 anchor [0.67,0.22,0.11,0].
R5 SKIPPED (reported): F2-4 coupling needs the c0 RIGHT eigenvector r0; at L=4 r0 is near-degenerate with the
   partner in the FULL operator (ratio ~0.9997) -> power iteration cannot separate (F2-4's own caveat). Needs
   r0 closed-form or deflation beyond the sanctioned within-block power iteration. Skipped per spec.
"""
import numpy as np, scipy.sparse as sp, time, json
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

def v3(n):
    if n == 0: return 99
    k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def build_tower(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    t0 = time.time()
    M, idx, n = build_M_gen(q, L, 2, raw)
    print(f"  build_M_gen L={L}: n={n} nnz={M.nnz} {time.time()-t0:.1f}s", flush=True)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states])
    tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsr().astype(np.float64)
    return Mt, gam[tw], D, qL

def power_iter(Mt, seed, tol=1e-10, maxit=3000):
    rng = np.random.default_rng(seed)
    v = np.abs(rng.standard_normal(Mt.shape[0])); v /= np.linalg.norm(v)
    rho_prev = 0.0
    for it in range(1, maxit + 1):
        w = Mt.dot(v)
        rho = float(v @ w)                       # Rayleigh quotient (v unit)
        res = float(np.linalg.norm(w - rho * v))  # ||Mv - rho v|| (v unit => /||v||=res)
        v = w / np.linalg.norm(w)
        if it % 25 == 0 or it < 5:
            print(f"    [seed {seed}] it {it}: rho={rho:.12f} res={res:.2e} d_rho={abs(rho-rho_prev):.1e}", flush=True)
        if res <= tol and abs(rho - rho_prev) <= 1e-13:
            return rho, v, it, res
        rho_prev = rho
    return rho, v, it, res

def main():
    L = 4
    import os
    CACHE = os.path.expanduser("~/g4_tower")
    if os.path.exists(CACHE + "_Mt.npz"):
        Mt = sp.load_npz(CACHE + "_Mt.npz"); gam = np.load(CACHE + "_gam.npy")
        print(f"  loaded cached tower: {Mt.shape[0]} states nnz={Mt.nnz}", flush=True)
    else:
        Mt, gam, D, qL = build_tower(L)
        sp.save_npz(CACHE + "_Mt.npz", Mt); np.save(CACHE + "_gam.npy", gam)
        print(f"  cached tower to {CACHE}", flush=True)
    nt = Mt.shape[0]
    print(f"  M_tower L=4: {nt} states, nnz={Mt.nnz}", flush=True)

    # ---- R1: partner eigenvalue, two starts ----
    t0 = time.time()
    rho1, v1, it1, r1 = power_iter(Mt, seed=12345)
    rho2, v2, it2, r2 = power_iter(Mt, seed=98765)
    print(f"  power-iter: start1 rho={rho1:.12f} ({it1} it, res {r1:.1e}); "
          f"start2 rho={rho2:.12f} ({it2} it, res {r2:.1e}); {time.time()-t0:.1f}s", flush=True)
    agree = abs(rho1 - rho2)
    rho = 0.5 * (rho1 + rho2)

    # ---- R2/R3 ----
    third = 1.0 / 3.0
    Delta = third - rho
    pred27 = 2 * 27.0 ** (-4)
    seq = {"L2": -2.911e-3, "L3": +9.958e-5, "L4": Delta}

    # ---- R4: gamma-level mass profile (use start-1 eigenvector; nonneg Perron) ----
    vp = np.abs(v1)
    levs = np.array([v3(int(g)) for g in gam])
    prof = {}
    tot = vp.sum()
    for lv in sorted(set(levs.tolist())):
        prof[int(lv)] = float(vp[levs == lv].sum() / tot)
    # also start-2 for robustness
    vp2 = np.abs(v2); prof2 = {int(lv): float(vp2[levs == lv].sum() / vp2.sum()) for lv in sorted(set(levs.tolist()))}

    out = dict(L=4, tower=int(nt), rho1=rho1, rho2=rho2, rho=rho, agree=agree,
               it=[it1, it2], res=[r1, r2], Delta=Delta, sign=("+" if Delta > 0 else "-"),
               absDelta=abs(Delta), pred27=pred27, ratio_meas_pred=abs(Delta) / pred27,
               profile=prof, profile_start2=prof2, seq=seq)
    print("RESULT_JSON " + json.dumps(out), flush=True)
    print(f"\n  === G4 L=4 PARTNER ===", flush=True)
    print(f"  R1 rho(M_tower,4) = {rho1:.10f} / {rho2:.10f}  (two starts agree to {agree:.1e})", flush=True)
    print(f"  R2 Delta_4 = 1/3 - rho = {Delta:+.6e}  (sign {out['sign']}); "
          f"sequence -2.911e-3 (L2) -> +9.958e-5 (L3) -> {Delta:+.4e} (L4)", flush=True)
    print(f"  R3 27^-L shot: |Delta_4| measured={abs(Delta):.4e}  predicted 2*27^-4={pred27:.4e}  "
          f"ratio meas/pred={abs(Delta)/pred27:.3f}  [ON-RECORD, zero weight until DERIVED; no fit]", flush=True)
    print(f"  R4 gamma-level mass profile (v3 level -> mass frac): {prof}", flush=True)
    print(f"     pre-reg cascade [2/3,2/9,2/27,1/27] = [0.6667,0.2222,0.0741,0.0370]; L3 anchor [0.67,0.22,0.11,0]", flush=True)
    print(f"  R5 SKIPPED (near-degenerate c0 right-eigvec r0 in full op; per F2-4 caveat + spec).", flush=True)
    # dump profile
    with open("g4_profile_L4.json", "w") as f:
        json.dump(out, f, indent=1)

if __name__ == "__main__":
    main()
