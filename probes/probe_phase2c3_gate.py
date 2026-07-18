"""
PROBE 2c3-GATE -- judge Wilson's blind rung-2 derivation (T1/T2/T3 -> P1-P4). Direct/exact, no eigensolves.

Branches from a v0 source, per channel e' (dst class): classify each move's destination dst=(a',e',gamma'):
   DEEP: gamma'==0 mod 3 (v'>=1);  UP: gamma'==-gamma mod3 (sigma'=+1);  DOWN: gamma'==gamma mod3 (sigma'=-1).
W-(x) := total pair-weight on DOWN branches (odd channels). ratio r(x) at beta*=3/5.
kappa(pop): (O,0)sig+ -> 3/4 ; (O,0)sig- -> 3 ; (E,0) -> 3/2.
P3 (master): r + kappa*W- == 4/9 exactly, all v0. P2: 3 pops MERGE onto one 3-value ladder {12,17,20}/49.
P1: r constant on W- (tau) level sets within each pop. P4: ladder = 4/9 - kappa*W-, SET predicted.
"""
import numpy as np
from fractions import Fraction
from collections import defaultdict

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def build(q, L, lam=0.5):
    qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    gam = np.array([s[2] for s in states]); a_arr = np.array([s[0] for s in states])
    erho = np.array([dl[(s[1] * pow(s[0], -1, qL)) % qL] for s in states])
    tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsc()
    return Mt, gam[tw], erho[tw], a_arr[tw], D, qL


def hbeta(erho, gam, a_arr, beta):
    nt = len(erho); par = erho % 2; v0 = (gam % 3 != 0)
    h = np.empty(nt)
    h[(par == 0) & v0] = 2/3; h[(par == 0) & ~v0] = 5/3
    h[(par == 1) & v0] = 5/6; h[(par == 1) & ~v0] = 4/3
    ov0 = (par == 1) & v0
    sig = np.zeros(nt); sig[ov0] = np.where((a_arr[ov0] + gam[ov0]) % 3 == 0, 1.0, -1.0)
    return h * (1 + beta * sig)


def run(L):
    Mt, gam, erho, a_arr, D, qL = build(3, L)
    nt = Mt.shape[0]
    beta = 3/5
    hb = hbeta(erho, gam, a_arr, beta)
    r = Mt.T.dot(hb) / hb
    Mc = Mt.tocsc()
    v0 = (gam % 3 != 0)
    # W-(x): DOWN weight over ODD channels; also W-_all over all channels (fallback)
    Wdn_odd = np.zeros(nt); Wdn_all = np.zeros(nt)
    indptr, indices, data = Mc.indptr, Mc.indices, Mc.data
    for x in range(nt):
        if not v0[x]:
            continue
        gx = gam[x] % 3
        for p in range(indptr[x], indptr[x + 1]):
            dst = indices[p]; wgt = data[p]
            gp = gam[dst]
            if gp % 3 == gx:                        # DOWN: gamma' == gamma mod 3 (=> v'=0)
                Wdn_all[x] += wgt
                if erho[dst] % 2 == 1:              # odd channel
                    Wdn_odd[x] += wgt
    # populations + kappa
    par = erho % 2; ov0 = (par == 1) & v0; ev0 = (par == 0) & v0
    sigp = ((a_arr + gam) % 3 == 0)
    pop = np.full(nt, "", dtype=object)
    pop[ov0 & sigp] = "O+"; pop[ov0 & ~sigp] = "O-"; pop[ev0] = "E"
    kappa = {"O+": 3/4, "O-": 3.0, "E": 3/2}
    log(f"\n{'='*72}\n## L={L}")
    for Wname, Wdn in [("W-(odd channels)", Wdn_odd), ("W-(all channels)", Wdn_all)]:
        res = np.zeros(nt)
        for x in range(nt):
            if v0[x]:
                res[x] = r[x] + kappa[pop[x]] * Wdn[x]
        maxdev = max(abs(res[x] - 4/9) for x in range(nt) if v0[x])
        log(f"   P3 [{Wname}]: max|r + kappa*W- - 4/9| over v0 = {maxdev:.3e}  "
            f"({'PASS (common intercept 4/9)' if maxdev < 1e-12 else 'no'})")
    # use the winning W- for the rest
    Wdn = Wdn_odd if max(abs(r[x] + kappa[pop[x]]*Wdn_odd[x] - 4/9) for x in range(nt) if v0[x]) < 1e-12 else Wdn_all
    which = "odd" if Wdn is Wdn_odd else "all"
    log(f"   -> using W-({which} channels) for P1/P2/P4")
    # P1: r constant on W- level sets within each pop; # distinct W- (trit)
    log(f"   P1/P4 per population (r = 4/9 - kappa*W-):")
    allvals = set()
    for pnm in ["O+", "O-", "E"]:
        mask = np.array([pop[i] == pnm for i in range(nt)])
        rv = sorted(set(np.round(r[mask], 9)))
        wv = sorted(set(np.round(Wdn[mask], 9)))
        allvals |= set(rv)
        # check r <-> W- bijective affine
        pairs = sorted(set((round(float(Wdn[i]), 9), round(float(r[i]), 9)) for i in range(nt) if mask[i]))
        log(f"      pop {pnm}: kappa={Fraction(kappa[pnm]).limit_denominator(10)}  #distinct W- = {len(wv)}  "
            f"#distinct r = {len(rv)}  r-values {[str(Fraction(v).limit_denominator(200)) for v in rv]}")
        for w_, r_ in pairs:
            pred = 4/9 - kappa[pnm] * w_
            log(f"          W-={Fraction(w_).limit_denominator(4000)} -> r={Fraction(r_).limit_denominator(200)} "
                f"(pred 4/9-kappa*W- = {Fraction(pred).limit_denominator(200)}, {'ok' if abs(pred-r_)<1e-6 else 'DEV'})")
    # P2: merge -> single 3-value ladder
    ladder = sorted(allvals)
    pred_set = {Fraction(12,49), Fraction(17,49), Fraction(20,49)}
    got_set = {Fraction(v).limit_denominator(200) for v in ladder}
    log(f"   P2 MERGE: union of r-values across pops = {[str(Fraction(v).limit_denominator(200)) for v in ladder]}  "
        f"({'== {12,17,20}/49 PASS' if got_set == pred_set else 'MISMATCH'})")


def main():
    log("# PROBE 2c3-GATE -- judge Wilson's blind rung-2 derivation (P1-P4). Direct/exact.")
    run(2); run(3)
    with open("logs/probe_phase2c3_gate_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
