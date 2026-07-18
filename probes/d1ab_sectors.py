"""
PROBE D1-A/B -- cap-sector recon (dense/direct at L<=3, sanctioned instrument).
A: dense eig of M_tower (L=2,3); extract partner (real top) + leading complex pair (eigenvectors).
   Gauge-character decomposition: mass Sum|f_k|^2 per sector k (Fourier in log2 a) + (k x carry-level) joint.
   Pre-reg: pair conjugate-symmetric mass(k)=mass(-k); partner k=0-dominated; pair in a single +-k pair (READOUT).
B: for the identified k (+ k=0): N_kappa(e',g,g') tables + R_k(s) tables (raw material for the pen).
"""
import numpy as np, json
from fractions import Fraction
from collections import defaultdict
from probe_phase2c0 import build_M_tower_and_coords, Nk_table, setup

def v3(n):
    if n == 0: return 99
    k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def sector_decomp(vec, twcoords, dl, D):
    blocks = defaultdict(lambda: np.zeros(D, dtype=complex))
    cnt = defaultdict(int)
    for i, (a, erho, gam) in enumerate(twcoords):
        blocks[(erho, gam)][dl[a]] = vec[i]; cnt[(erho, gam)] += 1
    full = all(c == D for c in cnt.values())
    massk = np.zeros(D); massk_lev = defaultdict(lambda: np.zeros(D))
    for (erho, gam), bv in blocks.items():
        fk = np.fft.fft(bv)                    # fk[k] = sum_ea bv[ea] exp(-2pi i k ea/D)
        p = np.abs(fk) ** 2
        massk += p; massk_lev[v3(int(gam))] += p
    tot = massk.sum()
    return massk / tot, {lv: (mk / tot) for lv, mk in massk_lev.items()}, full

def Rk_shell(D, lam=0.5):
    Z = 2 ** D - 1
    w = [Fraction(2 ** (D - m), Z) for m in range(1, D + 1)]
    om = np.exp(2j * np.pi / D)
    def Rk(k, s): return sum(w[m - 1] * w[((m - s - 1) % D)] * om ** (k * m) for m in range(1, D + 1))
    return Rk

def analyze(L):
    out = {"L": L}
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L)
    Md = Mt.toarray()
    ev, VR = np.linalg.eig(Md)
    order = np.argsort(-np.abs(ev))
    ev = ev[order]; VR = VR[:, order]
    # partner = top real; pair = top complex (|imag|>1e-9)
    isreal = np.abs(ev.imag) < 1e-9
    pi_ = int(np.where(isreal)[0][0])                      # first real = partner
    partner = ev[pi_].real
    cplx = np.where(~isreal & (ev.imag > 0))[0]
    pair_i = int(cplx[0]) if len(cplx) else None
    out["partner"] = partner
    out["pair"] = complex(ev[pair_i]) if pair_i is not None else None
    out["top8"] = [complex(e) for e in ev[:8]]
    # buffer = next distinct mode after partner+pair
    used = {pi_} | ({pair_i, pair_i} if pair_i is not None else set())
    buf_i = next((i for i in range(len(ev)) if i not in used and not (pair_i is not None and abs(ev[i]-np.conj(ev[pair_i]))<1e-12)), None)

    def report(name, i):
        mk, mkl, full = sector_decomp(VR[:, i], twcoords, dl, D)
        topk = sorted(range(D), key=lambda k: -mk[k])[:4]
        sym = float(np.max([abs(mk[k] - mk[(D - k) % D]) for k in range(D)]))
        r = {"name": name, "eval": complex(ev[i]), "blocks_full": full,
             "mass_by_k_top": {int(k): float(mk[k]) for k in topk},
             "k0_mass": float(mk[0]), "sym_maxdev": sym,
             "k_by_level": {int(lv): {int(k): round(float(mkl[lv][k]), 5) for k in topk} for lv in sorted(mkl)}}
        return r, topk
    rep = {}
    rep["partner"], tkp = report("partner", pi_)
    if pair_i is not None: rep["pair"], tkpair = report("pair", pair_i)
    else: tkpair = []
    if buf_i is not None: rep["buffer"], _ = report("buffer", buf_i)
    out["sectors"] = rep
    out["pair_dominant_k"] = [int(k) for k in tkpair[:2]] if tkpair else []
    return out, D, (tkpair[:2] if tkpair else [])

def dumpB(L, kids):
    """N_kappa + R_k(s) tables for k in {0} U kids."""
    q, qL, sub, D, dl, w, two = setup(L)
    N, _, _, _ = Nk_table(L)
    Rk = Rk_shell(D)
    ks = sorted(set([0] + list(kids)))
    lines = [f"# D1-B N_kappa + R_k tables, q=3 L={L}, D={D}. sectors k={ks}"]
    # R_k(s) table
    lines.append("## R_k(s)  [k ; s=0..D-1 ; value]")
    for k in ks:
        row = [f"R_{k}(s):"] + [f"{complex(Rk(k,s)).real:+.6e}{complex(Rk(k,s)).imag:+.6e}j" for s in range(D)]
        lines.append("\t".join(row))
    # N_kappa: for a sample of (e',g,g') tower transitions, list N_kappa[kappa] at kappa in ks
    lines.append("## N_kappa[(e',gamma,gamma')]  [kappa in ks]  (first 40 nonzero triples)")
    cnt = 0
    for (ep, g, gp), vecN in N.items():
        if g == 0 or gp == 0: continue
        vals = [f"k{kk}:{vecN[kk].real:+.4f}{vecN[kk].imag:+.4f}j" for kk in ks]
        lines.append(f"e'={ep} g={g} g'={gp}\t" + "\t".join(vals)); cnt += 1
        if cnt >= 40: break
    with open(f"d1b_NR_tables_q3_L{L}.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    return len(N)

def main():
    allout = {}
    for L in [2, 3]:
        o, D, kids = analyze(L)
        nN = dumpB(L, kids)
        o["Ncount"] = nN
        allout[f"L{L}"] = o
        print(f"\n===== L={L} (D={D}) =====", flush=True)
        print(f"  partner = {o['partner']:.10f}", flush=True)
        print(f"  leading pair = {o['pair']}  (|.|={abs(o['pair']):.6f}, arg={np.angle(o['pair']):.6f} rad)"
              if o['pair'] else "  no complex pair", flush=True)
        print(f"  top-8 |eval|: {[round(abs(e),5) for e in o['top8']]}", flush=True)
        s = o["sectors"]
        print(f"  PARTNER sector: k0_mass={s['partner']['k0_mass']:.6f}  top-k={s['partner']['mass_by_k_top']}", flush=True)
        if 'pair' in s:
            print(f"  PAIR sector: top-k={s['pair']['mass_by_k_top']}  k=0 mass={s['pair']['k0_mass']:.6f}  "
                  f"conj-sym maxdev(mass_k vs mass_-k)={s['pair']['sym_maxdev']:.2e}", flush=True)
            print(f"  PAIR (k x carry-level): {s['pair']['k_by_level']}", flush=True)
            print(f"  => pair dominant k = {o['pair_dominant_k']}  (names the effective-model basis)", flush=True)
        if 'buffer' in s:
            print(f"  BUFFER {s['buffer']['eval']}: top-k={s['buffer']['mass_by_k_top']}", flush=True)
    with open("d1a_sectors.json", "w") as f:
        json.dump(allout, f, indent=1, default=lambda z: [z.real, z.imag] if isinstance(z, complex) else str(z))
    print("\nRESULT_JSON_A written d1a_sectors.json + d1b_NR_tables_q3_L{2,3}.txt", flush=True)

if __name__ == "__main__":
    main()
