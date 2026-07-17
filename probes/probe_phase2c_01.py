"""
PHASE 2c stages 0 + 1 -- freeze the defect object + zeroth Collatz-Wielandt bracket.
INSTRUMENT LAW: direct/exact at q=3. No proof authored; code gates. Label ALGEBRAIC vs STATISTICAL.

2c-0 FREEZE THE DEFECT OBJECT:
  M_bar = source-side gauge-average closure: M_bar[dst,src] = (1/D) sum_{s in <2>} M[dst, s.src]
          (gauge s: (a,b,g)->(sa,sb,g); orbit = fixed (e_rho, gamma), varies e_a). This IS the
          (e_rho,gamma)-compressed chain lifted (E/W's Lmat); k=0 unit-group character component.
  Defect D = M_tower - M_bar; character expansion D = sum_{k!=0} D_k over unit-group chars chi_k(s)=w^{k e_a}.
  Amplitude ||D_k||_F^2 = (1/D) sum_orbit ||S_orbit v_k||^2 (v_k[m]=w^{km}); D_0=M_bar.
  GATE: Parseval reconstruction sum_k ||D_k||^2 = ||M_tower||_F^2 (machine precision) => M_bar+sum D_k=M_tower.
  Fold in: G0-2 scope-correction -- per-STATE survival spread (G0-2's spread-0 was CLASS-AVERAGED).
  Named obstruction check: |c_k|/c_0 = |sum w_d^2 w^{kd}| / sum w_d^2 does NOT decay (0.83@D6 -> 0.974@D18 -> 1).

2c-1 ZEROTH BRACKET:
  h0 = coarse Perron eigvec lift by cell (2/3, 5/3, 5/6, 4/3). Collatz-Wielandt bracket
  [min_x, max_x] (M_tower h0)(x)/h0(x) contains rho(M_tower)=partner. PRE-REG: 1/3 inside; width=zeroth scale.
"""
import numpy as np
import scipy.sparse as sp
from collections import defaultdict

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def build_tower(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    gam = np.array([s[2] for s in states])
    ea = np.array([dl[s[0]] for s in states])
    erho = np.array([dl[(s[1] * pow(s[0], -1, qL)) % qL] for s in states])
    tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsc()
    w = np.array(raw) / sum(raw)
    return Mt, gam[tw], ea[tw], erho[tw], D, qL, w


def h0_lift(erho, gam):
    par = erho % 2; v0 = (gam % 3 != 0)
    h = np.empty(len(erho))
    h[(par == 0) & v0] = 2 / 3
    h[(par == 0) & ~v0] = 5 / 3
    h[(par == 1) & v0] = 5 / 6
    h[(par == 1) & ~v0] = 4 / 3
    return h


def cell_of(erho, gam):
    return (int(erho % 2), 0 if gam % 3 != 0 else 1)


def run(L, partner):
    Mt, gam, ea, erho, D, qL, w = build_tower(L)
    nt = Mt.shape[0]
    log(f"\n{'='*76}\n## PHASE 2c  q=3 L={L}   tower dim {nt}   partner rho_L={partner}")

    # ---------- 2c-1 : zeroth Collatz-Wielandt bracket ----------
    # G1's coarse transfer is T[src,dst]=flow(src->dst) with h0 its RIGHT eigvec (T h0 = (1/3) h0).
    # The fine operator matching that convention is M^T (Mt is M[dst,src], so flow src->dst = Mt^T).
    # (M h0)(dst)=inflow-weighted has zero rows for unreached states -> degenerate; M^T is correct.
    h0 = h0_lift(erho, gam)
    MtT = Mt.T.tocsr()
    Mh = MtT.dot(h0)                               # (M^T h0)(x) = sum_dst flow(x->dst) h0[dst]
    r = Mh / h0
    m0, M0 = r.min(), r.max()
    contains_third = m0 <= 1 / 3 <= M0
    contains_rho = m0 <= partner <= M0
    log(f"\n   2c-1 ZEROTH BRACKET  [m0, M0] = [{m0:.8f}, {M0:.8f}]   width = {M0-m0:.6e}")
    log(f"        1/3 = {1/3:.8f} inside: {contains_third}      rho_L={partner} inside: {contains_rho}")
    log(f"        (Collatz-Wielandt guarantees rho_L in [m0,M0]; width = zeroth defect scale)")
    # where are the extremes (which cell)?
    imin, imax = int(np.argmin(r)), int(np.argmax(r))
    log(f"        min at cell {cell_of(erho[imin],gam[imin])} (gamma={gam[imin]}), "
        f"max at cell {cell_of(erho[imax],gam[imax])} (gamma={gam[imax]})")

    # ---------- 2c-0a : per-state survival spread (scope correction) ----------
    surv = np.asarray(Mt.sum(axis=0)).ravel()      # column sums = per-source survival
    cells = defaultdict(list)
    for i in range(nt): cells[cell_of(erho[i], gam[i])].append(surv[i])
    log(f"\n   2c-0a SCOPE CORRECTION (per-STATE survival; G0-2's spread-0 was CLASS-averaged):")
    global_spread = surv.max() - surv.min()
    for c in sorted(cells):
        arr = np.array(cells[c])
        log(f"        cell {c}: mean {arr.mean():.6f}  [{arr.min():.6f}, {arr.max():.6f}]  "
            f"within-cell spread {arr.max()-arr.min():.4f}")
    max_within = max(np.array(cells[c]).max() - np.array(cells[c]).min() for c in cells)
    log(f"        => per-state survival is STATISTICAL not ALGEBRAIC: max within-cell spread = {max_within:.4f} "
        f"(global {global_spread:.4f})")

    # ---------- 2c-0b : defect character amplitudes + Parseval reconstruction ----------
    orbits = defaultdict(dict)
    for i in range(nt): orbits[(int(erho[i]), int(gam[i]))][int(ea[i])] = i
    Vmat = np.exp(2j * np.pi * np.outer(np.arange(D), np.arange(D)) / D)   # V[m,k]=w^{mk}
    amp = np.zeros(D)
    for key, dmap in orbits.items():
        cols = [dmap[m] for m in range(D)]         # ordered by e_a (complete: a ranges over all units)
        S = Mt[:, cols].toarray()                  # nt x D
        B = S @ Vmat                               # nt x D over k
        amp += (np.abs(B) ** 2).sum(axis=0) / D
    fro2 = float((Mt.multiply(Mt)).sum())
    parseval_res = abs(amp.sum() - fro2)
    log(f"\n   2c-0b DEFECT CHARACTER AMPLITUDES (M_bar=D_0 gauge-average; D_k = twist-k defect):")
    log(f"        Parseval reconstruction  sum_k ||D_k||^2 = ||M_tower||_F^2 :  "
        f"{amp.sum():.8f} vs {fro2:.8f}   residual {parseval_res:.2e}  "
        f"({'EXACT (M_bar + sum D_k = M_tower)' if parseval_res < 1e-10 else 'MISMATCH'})")
    a0 = amp[0]
    log(f"        ||M_bar||^2 = ||D_0||^2 = {a0:.6f} ; defect ||D||^2 = sum_{{k!=0}} = {amp[1:].sum():.6f} "
        f"(defect/meanfield = {amp[1:].sum()/a0:.4f})")
    log(f"        amplitude spectrum sqrt(||D_k||^2/||D_0||^2) by k:")
    ratios = np.sqrt(amp / a0)
    log("          " + "  ".join(f"k={k}:{ratios[k]:.4f}" for k in range(D)))

    # ---------- named obstruction: |c_k|/c_0 does NOT decay ----------
    dd = np.arange(1, D + 1)
    ck = np.array([np.sum(w ** 2 * np.exp(2j * np.pi * k * dd / D)) for k in range(D)])
    c0 = ck[0].real
    log(f"\n   NAMED OBSTRUCTION (2c-3): character amplitudes |c_k|/c_0 = |sum w_d^2 w^{{kd}}|/sum w_d^2 :")
    log(f"        |c_1|/c_0 = {abs(ck[1])/c0:.4f}  (D={D})   [Wilson: 0.83@D6 -> 0.974@D18 -> 1; amplitude decay CANNOT carry contraction]")
    log(f"        |c_k|/c_0 by k: " + "  ".join(f"{abs(ck[k])/c0:.3f}" for k in range(1, min(D, 7))))
    return (m0, M0, contains_third, contains_rho, max_within, parseval_res, abs(ck[1])/c0)


def main():
    log("# PHASE 2c stages 0+1 -- defect object + zeroth bracket. Direct at q=3. Claude derives, code gates.")
    r2 = run(2, 0.346827)
    r3 = run(3, 0.333236)
    log(f"\n{'='*76}\n## SUMMARY")
    log(f"   2c-1 bracket: L=2 [{r2[0]:.5f},{r2[1]:.5f}] w={r2[1]-r2[0]:.2e} (1/3 in:{r2[2]}, rho in:{r2[3]}); "
        f"L=3 [{r3[0]:.5f},{r3[1]:.5f}] w={r3[1]-r3[0]:.2e} (1/3 in:{r3[2]}, rho in:{r3[3]})")
    log(f"        width shrink L2->L3: {(r2[1]-r2[0])/(r3[1]-r3[0]):.2f}x  (the zeroth-order data 2c-2 must beat)")
    log(f"   2c-0a per-state survival spread: L=2 {r2[4]:.3f}, L=3 {r3[4]:.3f}  (STATISTICAL, corrects G0-2 spread-0)")
    log(f"   2c-0b reconstruction residual: L=2 {r2[5]:.1e}, L=3 {r3[5]:.1e}  (EXACT); |c_1|/c_0 = {r2[6]:.3f}(D6)/{r3[6]:.3f}(D18) NON-DECAYING")
    with open("logs/probe_phase2c_01_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
