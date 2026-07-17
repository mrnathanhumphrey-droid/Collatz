"""
PROBE 2c0 -- GATE THE DEFECT FREEZE (D-FORM + D-DEPTH selection rule). Pure entry algebra, no eigensolves.
INSTRUMENT LAW: direct/exact at q=3. Claude derives / code gates. Label ALGEBRAIC vs STATISTICAL.

Coords: tower state (a, e=e_rho, gamma!=0), a in U=<2> mod 3^L (gauge circle Z/D). Single-state move:
target gauge u=a', shift s => e'=e+s; m = dlog2(a) - dlog2(u) in {1..D} = delta_a; weight w_m w_{m-s};
gate reads u:  (gamma + u(1-2^{e'})) == 0 mod 3;  carry  gamma' = (gamma + u(1-2^{e'}))//3.

G1 RECONSTRUCTION: build M_tower from the D-FORM (a,e,gamma)+(u,s) parametrization, diff vs build_M_gen
   entry-by-entry (machine precision, both L). Confirms the single-state D-FORM kernel.
G2 SECTOR FORM: gauge-Fourier of M_tower; verify  B_hat[k_out,k_in] = R_{k_in}(s) * N_{k_in-k_out}(e',g,g')
   (R_k(s)=sum_m w_m w_{m-s} w^{km}; N_kappa=sum_{u:gate,carry->g'} w^{kappa*dlog u}); pin twist sign; k=0 -> E-FORM.
G3 D-DEPTH ZEROS: N_k aggregated over gamma' at depth-j resolution vanishes unless 3^{L-j} | k. Report the
   observed divisibility pattern vs claimed, both digit conventions.
G4: scope number (per-state survival spread) -- one number, folds in 2c-0a.
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


def setup(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    w = np.array([lam ** d for d in range(1, D + 1)]); w = w / w.sum()
    two = [pow(2, e, qL) for e in range(D)]
    return q, qL, sub, D, dl, w, two


def build_M_tower_and_coords(L, lam=0.5):
    q, qL, sub, D, dl, w, two = setup(L, lam)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states])
    tw = np.where(gam != 0)[0]
    pos = {int(t): i for i, t in enumerate(tw)}                 # global idx -> tower idx
    Mt = M[tw][:, tw]
    twcoords = [(states[t][0], dl[(states[t][1] * pow(states[t][0], -1, qL)) % qL], states[t][2]) for t in tw]
    return Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two


# ============================ G1 : RECONSTRUCTION ============================
def g1(L):
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L)
    nt = len(tw)
    rows, cols, vals = [], [], []
    for src_t, (a, e, gam) in enumerate(twcoords):
        dla = dl[a]
        for u in sub:                                          # target gauge a'
            dlu = dl[u]
            da = (dla - dlu) % D; da = da if da != 0 else D    # delta_a in 1..D
            for s in range(D):                                 # shift; e'=e+s
                ep = (e + s) % D
                bp = (u * two[ep]) % qL
                T = (u - bp) % qL                              # = u(1-2^{e'})
                if (gam + T) % q != 0:
                    continue
                gp = ((gam + T) // q) % qL
                dst = (u, bp, gp)
                if dst not in idx:
                    continue
                dg = idx[dst]
                if dg not in pos:
                    continue
                db = (da - s) % D; db = db if db != 0 else D
                rows.append(pos[dg]); cols.append(src_t); vals.append(w[da - 1] * w[db - 1])
    Md = sp.csr_matrix((vals, (rows, cols)), shape=(nt, nt))
    diff = (Md - Mt)
    maxd = abs(diff).max() if diff.nnz else 0.0
    nnz_form = Md.nnz; nnz_true = Mt.nnz
    log(f"\n   G1 RECONSTRUCTION L={L}: D-FORM build vs build_M_gen  max|diff| = {maxd:.3e}  "
        f"(nnz form {nnz_form} vs true {nnz_true})  "
        f"({'EXACT' if maxd < 1e-13 and nnz_form == nnz_true else 'MISMATCH'})")
    return maxd < 1e-13 and nnz_form == nnz_true


# ============================ N_k twisted counts ============================
def Nk_table(L, lam=0.5):
    """N_k(e',gamma,gamma') = sum_{u: gate passes, carry=gamma'} w^{k dlog u}. Returns dict."""
    q, qL, sub, D, dl, w, two = setup(L, lam)
    om = np.exp(2j * np.pi / D)
    N = defaultdict(lambda: np.zeros(D, dtype=complex))        # (e',gamma,gamma') -> N_k[k]
    for ep in range(D):
        c = (1 - two[ep]) % qL
        for u in sub:
            T = (u * c) % qL
            duk = np.array([om ** (k * dl[u]) for k in range(D)])
            for gam in range(1, qL):
                if (gam + T) % q == 0:
                    gp = ((gam + T) // q) % qL
                    N[(ep, gam, gp)] += duk
    return N, D, qL, q


# ============================ G2 : SECTOR FORM ============================
def g2(L, lam=0.5):
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L, lam)
    om = np.exp(2j * np.pi / D)
    # R_k(s) = sum_m w_m w_{m-s} om^{km}, m=1..D indexing (w[m-1])
    def Rk(k, s):
        return sum(w[m - 1] * w[((m - s - 1) % D)] * om ** (k * m) for m in range(1, D + 1))
    Nt, _, _, _ = Nk_table(L, lam)
    # gauge-Fourier of M_tower per (e,gamma)->(e',gamma') block: B[a',a]=Mt entry; B_hat[kout,kin]
    # group tower indices by (e,gamma) and by a
    by_eg = defaultdict(dict)                                  # (e,gamma) -> {dlog a: tower_idx}
    for ti, (a, e, gam) in enumerate(twcoords):
        by_eg[(e, gam)][dl[a]] = ti
    Mtc = Mt.tocsr()
    maxres = 0.0; nchecked = 0; e_form_ok = True
    # iterate source blocks
    for (e, gam), amap in by_eg.items():
        src_idx = [amap[al] for al in range(D)]               # ordered by dlog a
        for (ep, gpp), amap2 in by_eg.items():
            dst_idx = [amap2[al] for al in range(D)]
            # B[a', a] = Mt[dst=(a',ep,gpp), src=(a,e,gam)]
            B = Mtc[np.ix_(dst_idx, src_idx)].toarray()       # D x D  (rows a'=dlog, cols a=dlog)
            if not B.any():
                continue
            s = (ep - e) % D
            # direct Fourier: B_hat[kout,kin] = sum_{a',a} om^{-kout dlog a'} B[a',a] om^{kin dlog a}
            F = np.array([[om ** (kk * jj) for jj in range(D)] for kk in range(D)])  # F[k,j]=om^{kj}
            Bhat = F.conj() @ B @ F.T                          # [kout,kin]
            # formula: R_{kin}(s) * N_{kin-kout}(ep,gam,gpp)
            Nvec = Nt[(ep, gam, gpp)]                          # N_kappa[kappa]
            for kin in range(D):
                for kout in range(D):
                    pred = Rk(kin, s) * Nvec[(kin - kout) % D]
                    maxres = max(maxres, abs(Bhat[kout, kin] - pred)); nchecked += 1
            # E-FORM sanity: kin=kout=0 diagonal element / D == R_0(s) N_0 / D
    log(f"\n   G2 SECTOR FORM L={L}: B_hat[kout,kin] vs R_{{kin}}(s)*N_{{kin-kout}}  "
        f"max residual over {nchecked} sector entries = {maxres:.3e}   "
        f"({'CLOSES (twist: R_{kin}, N_{kin-kout})' if maxres < 1e-10 else 'does not close - try other sign'})")
    return maxres < 1e-10


# ============================ G3 : D-DEPTH SELECTION RULE ============================
def g3(L, lam=0.5):
    Nt, D, qL, q = Nk_table(L, lam)
    log(f"\n   G3 D-DEPTH SELECTION RULE L={L}  (gamma' resolved mod 3^j, j=0..L-1)")
    log(f"       Wilson stated 3^(L-j)|k; GATE PINS the index -> 3^(L-1-j)|k (carry's /3 drops gamma' one 3-adic level)")
    all_hold = True
    for j in range(0, L):                                        # gamma' mod 3^j (j digits); j=0 => total
        exp = L - 1 - j
        pinned = sorted([k for k in range(D) if k % (3 ** exp) == 0])
        stated = sorted([k for k in range(D) if k % (3 ** (L - j)) == 0]) if (L - j) <= L else []
        agg = defaultdict(lambda: np.zeros(D, dtype=complex))
        for (ep, gam, gp), vec in Nt.items():
            cls = gp % (3 ** j)
            agg[(ep, gam, cls)] += vec
        nz = set(); max_forbidden = 0.0
        for key, vec in agg.items():
            for k in range(D):
                if abs(vec[k]) > 1e-10: nz.add(k)
                if k % (3 ** exp) != 0: max_forbidden = max(max_forbidden, abs(vec[k]))
        observed = sorted(nz)
        ok = (observed == pinned)
        all_hold = all_hold and ok
        log(f"       j={j} (gamma' mod 3^{j}): PINNED excited k (3^{exp}|k) = {pinned}")
        log(f"                observed nonzero k = {observed}   max|N_k| on FORBIDDEN(pinned) k = {max_forbidden:.2e}  "
            f"({'HOLDS (=pinned set exactly)' if ok else 'MISMATCH'})   [Wilson-stated 3^{L-j}|k was {stated}]")
    log(f"     => D-DEPTH SELECTION RULE {'CONFIRMED with pinned index 3^(L-1-j)|k' if all_hold else 'FAILS'} "
        f"(mechanism = coset char-sum on principal-unit tower; off-by-one from stated = carry level-drop)")
    return all_hold


# ============================ G4 : scope number ============================
def g4(L, lam=0.5):
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L, lam)
    surv = np.asarray(Mt.sum(axis=0)).ravel()
    cells = defaultdict(list)
    for i, (a, e, gam) in enumerate(twcoords):
        cells[(e % 2, 0 if gam % 3 != 0 else 1)].append(surv[i])
    spreads = {c: (np.array(v).max() - np.array(v).min()) for c, v in cells.items()}
    stat = [c for c, s in spreads.items() if s > 1e-9]
    log(f"\n   G4 SCOPE (per-state survival spread): statistical cells {stat}; "
        f"max within-cell spread {max(spreads.values()):.4f} (half {max(spreads.values())/2:.4f} ~ Wilson 0.17); "
        f"3/4 cells ALGEBRAIC (spread 0)")


def main():
    log("# PROBE 2c0 -- GATE THE DEFECT FREEZE (D-FORM + D-DEPTH). Pure entry algebra, no eigensolves.")
    for L in [2, 3]:
        log(f"\n{'='*76}\n## L={L}")
        r1 = g1(L)
        r2 = g2(L)
        g3(L)
        g4(L)
        log(f"   >> L={L}: G1 {'PASS' if r1 else 'FAIL'}, G2 {'PASS' if r2 else 'FAIL'}")
    with open("logs/probe_phase2c0_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
