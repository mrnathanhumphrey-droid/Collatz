"""
PROBE J2 -- THE L=4 ARM (SpMV projections, NO dense eig). Compute d(1,4), d(2,4) by the same QSD
gauge-Fourier collapse used at L=2,3, but WITHOUT forming a dense operator:
  reduced A[ko,ki] = L_ko . (M . v_ki),  v_ki[state] = R[block] * conj(F[ki, gauge]),
  L_ko[state] = Lp[block] * F[ko, gauge],  F[k,j]=exp(-2pi i k j/D)/sqrt(D).
Lp,R = left/right block-QSD (Perron of the block-marginal B_block, power iteration).
PRE-REGISTERED (direction only): (i) |d(1,4)|/|sigma(th1)| in (1, 1.0201) closer to 1 from above;
(ii) phase ratio arg(d(1,4))/th1 in (0.9434, 1) closer to 1 from below (banked block-level 0.993);
(iii) off-diagonal couplings <= the L=3 scale (diagonal dominance persists). th1=2pi/3^(L-1).
INSTRUMENT LAW: SpMV / within-block power iteration; NO ARPACK/shift-invert; no dense eig; no fit.
"""
import numpy as np, scipy.sparse as sp, time, json, os
from probe_phase2c0 import build_M_tower_and_coords

CACHE = os.path.expanduser("~/j2_L4")


def build_or_load(L=4, lam=0.5):
    if os.path.exists(CACHE + "_Mt.npz"):
        Mt = sp.load_npz(CACHE + "_Mt.npz")
        z = np.load(CACHE + "_coords.npz")
        print(f"  loaded cache: tower {Mt.shape[0]} nnz {Mt.nnz}", flush=True)
        return Mt.tocsr(), z["block_id"], z["gauge_pos"], int(z["D"]), int(z["Nb"])
    t0 = time.time()
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L, lam)
    Mt = Mt.tocsr().astype(np.float64)
    print(f"  built L={L}: tower {Mt.shape[0]} nnz {Mt.nnz}  ({time.time()-t0:.0f}s)", flush=True)
    blocks = {}
    block_id = np.empty(len(twcoords), dtype=np.int64)
    gauge_pos = np.empty(len(twcoords), dtype=np.int64)
    for i, (a, e, g) in enumerate(twcoords):
        key = (e, g)
        if key not in blocks:
            blocks[key] = len(blocks)
        block_id[i] = blocks[key]
        gauge_pos[i] = dl[a]
    Nb = len(blocks)
    sp.save_npz(CACHE + "_Mt.npz", Mt)
    np.savez(CACHE + "_coords.npz", block_id=block_id, gauge_pos=gauge_pos, D=D, Nb=Nb)
    print(f"  cached; Nb={Nb} D={D}", flush=True)
    return Mt, block_id, gauge_pos, D, Nb


def block_marginal(Mt, block_id, D, Nb):
    """B_block[bo,bi] = (1/D) sum_{states in bo,bi} Mt = mean-field block transfer G."""
    coo = Mt.tocoo()
    B = np.zeros((Nb, Nb))
    np.add.at(B, (block_id[coo.row], block_id[coo.col]), coo.data)
    return B / D


def perron(B, tol=1e-13, maxit=5000):
    n = B.shape[0]
    vR = np.abs(np.random.default_rng(1).standard_normal(n)); vR /= np.linalg.norm(vR)
    vL = np.abs(np.random.default_rng(2).standard_normal(n)); vL /= np.linalg.norm(vL)
    rho = 0.0
    for it in range(maxit):
        w = B @ vR; rho = vR @ w; vR = w / np.linalg.norm(w)
        wl = B.T @ vL; vL = wl / np.linalg.norm(wl)
        if it > 5 and np.linalg.norm(B @ vR - rho * vR) < tol:
            break
    Lp = vL / (vL @ vR)                                          # normalize Lp.R = 1
    return rho, Lp, vR, it


def reduced_entry(Mt, Lp, R, block_id, gauge_pos, D, ko, ki):
    F = np.exp(-2j * np.pi * np.arange(D) / D) / np.sqrt(D)      # F[k,j]=exp(-2pi i k j/D)/sqrt(D); use powers
    v = R[block_id] * np.conj(np.exp(-2j * np.pi * ki * gauge_pos / D) / np.sqrt(D))
    w = Mt.dot(v)
    u = Lp[block_id] * (np.exp(-2j * np.pi * ko * gauge_pos / D) / np.sqrt(D))
    return complex(u @ w)


def reduced_matrix(Mt, Lp, R, block_id, gauge_pos, D, rungs):
    n = len(rungs)
    A = np.zeros((n, n), dtype=complex)
    for a, ko in enumerate(rungs):
        for b, ki in enumerate(rungs):
            A[a, b] = reduced_entry(Mt, Lp, R, block_id, gauge_pos, D, ko, ki)
    return A


def main():
    L = 4
    print("# PROBE J2 -- L=4 ARM (SpMV, no dense eig). d(1,4), d(2,4). Pre-committed directions.")
    Mt, block_id, gauge_pos, D, Nb = build_or_load(L)
    print(f"  tower {Mt.shape[0]} = Nb {Nb} x D {D}", flush=True)

    B = block_marginal(Mt, block_id, D, Nb)
    rho, Lp, R, it = perron(B)
    print(f"  block-QSD Perron rho={rho:.10f}  ({it} it)  [g4 partner rho_4=0.3334999 expected]", flush=True)

    rungs1 = [(1 * 3 ** j) % D for j in range(L)]                 # [1,3,9,27]
    rungs2 = [(2 * 3 ** j) % D for j in range(L)]                 # [2,6,18,0]
    A1 = reduced_matrix(Mt, Lp, R, block_id, gauge_pos, D, rungs1)
    A2 = reduced_matrix(Mt, Lp, R, block_id, gauge_pos, D, rungs2)

    ev1 = sorted(np.linalg.eig(A1)[0], key=lambda z: -abs(z))
    ev2 = sorted(np.linalg.eig(A2)[0], key=lambda z: -abs(z))
    # d(1,4) = k=+-1 dominant (doublet-diagonal); d(2,4) = k=+-2 dominant NON-DC (m=2 analog)
    d1 = ev1[0]
    d2 = next((z for z in ev2 if abs(z - rho) > 0.02), ev2[0])    # skip the DC/partner mode
    th1 = 2 * np.pi / 3 ** (L - 1)
    sig1 = (1 / 3) * np.cos(th1 / 2) ** 2

    offdiag1 = max(abs(A1[a, b]) for a in range(L) for b in range(L) if a != b)
    offdiag2 = max(abs(A2[a, b]) for a in range(L) for b in range(L) if a != b)

    print(f"\n  === J2 L=4 ARM ===", flush=True)
    print(f"  rungs k=+-1 {rungs1}   A1 diag: " + ", ".join(f"{A1[i,i].real:+.5f}{A1[i,i].imag:+.5f}j" for i in range(L)), flush=True)
    print(f"  rungs k=+-2 {rungs2}   A2 diag: " + ", ".join(f"{A2[i,i].real:+.5f}{A2[i,i].imag:+.5f}j" for i in range(L)), flush=True)
    print(f"  d(1,4) = {d1.real:+.6f}{d1.imag:+.6f}j  |d|={abs(d1):.6f}  arg={np.angle(d1):.6f}", flush=True)
    print(f"  d(2,4) = {d2.real:+.6f}{d2.imag:+.6f}j  |d|={abs(d2):.6f}  arg={np.angle(d2):.6f}", flush=True)
    print(f"  sigma(th1) = {sig1:.6f}  th1 = {th1:.6f}", flush=True)
    r_i = abs(d1) / sig1
    r_ii = np.angle(d1) / th1
    print(f"\n  (i)  |d(1,4)|/sigma(th1) = {r_i:.5f}   pre-reg (1, 1.0201): "
          f"{'PASS' if 1 < r_i < 1.0201 else 'DEV'}", flush=True)
    print(f"  (ii) phase ratio arg(d1)/th1 = {r_ii:.5f}   pre-reg (0.9434, 1): "
          f"{'PASS' if 0.9434 < r_ii < 1 else 'DEV'}", flush=True)
    print(f"  (iii) off-diag coupling: k=+-1 {offdiag1:.5f}, k=+-2 {offdiag2:.5f}  "
          f"(L=3 scale ~2e-3..1e-3; pre-reg <= L=3 scale)", flush=True)
    out = dict(L=4, tower=int(Mt.shape[0]), Nb=Nb, D=D, block_rho=rho,
               d1=[d1.real, d1.imag], d2=[d2.real, d2.imag], sigma_th1=sig1, th1=th1,
               ratio_mod=r_i, ratio_phase=r_ii, offdiag1=offdiag1, offdiag2=offdiag2,
               A1_diag=[[A1[i,i].real, A1[i,i].imag] for i in range(L)],
               A2_diag=[[A2[i,i].real, A2[i,i].imag] for i in range(L)])
    with open("outputs/judge_L4_arm.json", "w") as f:
        json.dump(out, f, indent=1)
    print("RESULT_JSON " + json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
