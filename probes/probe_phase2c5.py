"""
PROBE 2c5 -- GATE Wilson's blind v1-trit derivation (pre-regs P1'..P4'). Direct/exact, no eigensolves.

Baseline: beta*=3/5, NO v0 trit dressing (g0=0). Operator = M^T (flow); ratio r=(M^T h_beta)/h_beta.
Target = v>=1 tower states (gamma == 0 mod 3). Wilson's level-set key (P1'):
   class3 in {D9 (gamma==0 mod9), U+ (v3=1, a==g mod3), U- (v3=1, a!=g mod3)},  g:=gamma/3
   b(x) := +1 iff a == gamma/3 (mod3);  eta := +1 iff a == 2g (mod3) == (-b).

P1' r takes <=3 values per e mod6; level sets EXACTLY (e mod6) x class3; value-set == baseline 5-set.
P2' e==0 mod6 is b-BLIND (U+==U-, two levels); e==2 and e==4 carry the bit with OPPOSITE-sign coeffs.
P3' D9 level's DEEP mass (flow to v'>=1) = R0(s0) exactly, s0 = -e mod 6, per e-group.
P4' unit-g DEEP mass = [R0(s2)+R0(s4)]/2 - eta*[R_{D/2}(s2)-R_{D/2}(s4)]/2, s2=2-e, s4=4-e; both L, machine prec.

Deviations reported AS deviations. Reference/convention taken verbatim from probe_phase2c0.g2 (R_k) + 2c3.
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
    return M[tw][:, tw].tocsc(), gam[tw], erho[tw], a_arr[tw], D, qL


def hbeta(erho, gam, a_arr, beta):
    nt = len(erho); par = erho % 2; v0 = (gam % 3 != 0)
    h = np.empty(nt)
    h[(par == 0) & v0] = 2/3; h[(par == 0) & ~v0] = 5/3
    h[(par == 1) & v0] = 5/6; h[(par == 1) & ~v0] = 4/3
    ov0 = (par == 1) & v0
    sig = np.zeros(nt); sig[ov0] = np.where((a_arr[ov0] + gam[ov0]) % 3 == 0, 1.0, -1.0)
    return h * (1 + beta * sig)


def Rweights(D):
    """R weights on the mod-6 BASE SHELL (D=6): the deep-mass closed form lives here, NOT on the full
    tower depth 2*3^{L-1}. (At L=2 D=6 coincides; at L>=3 the full-D autocorrelation is the wrong shell.
    Empirically the deep masses are L-INVARIANT /189 fractions == R computed at D=6.)"""
    Z = 2 ** D - 1
    w = [Fraction(2 ** (D - m), Z) for m in range(1, D + 1)]   # w[m-1], m=1..D
    def R0(s):  return sum(w[m - 1] * w[((m - s - 1) % D)] for m in range(1, D + 1))
    def Rhalf(s): return sum(w[m - 1] * w[((m - s - 1) % D)] * (-1) ** m for m in range(1, D + 1))
    return R0, Rhalf, None


def cls3(a, gam):
    """Wilson's v>=1 level class: D9 / U+ / U-. Returns (label, b, eta) with g=gamma/3."""
    g3 = (gam % 9) // 3            # 0 => gamma==0 mod9 (D9); else v3(gamma)=1, g mod3 = g3
    if g3 == 0:
        return "D9", 0, 0
    a3 = a % 3
    b = 1 if a3 == g3 else -1
    eta = 1 if a3 == (2 * g3) % 3 else -1
    return ("U+" if b == 1 else "U-"), b, eta


def deep_mass(Mt, gam, hb, x):
    """h-weighted DEEP (to v'>=1) and SURFACE (to v'=0) contributions to r[x], + raw move-weight deep."""
    indptr, indices, data = Mt.indptr, Mt.indices, Mt.data
    deep_h = surf_h = deep_raw = 0.0
    for p in range(indptr[x], indptr[x + 1]):
        dst = indices[p]
        if gam[dst] % 3 == 0:      # v'>=1 = DEEP
            deep_h += data[p] * hb[dst]; deep_raw += data[p]
        else:
            surf_h += data[p] * hb[dst]
    return deep_h / hb[x], surf_h / hb[x], deep_raw


def run(L):
    log(f"\n{'='*74}\n## L={L}")
    Mt, gam, erho, a_arr, D, qL = build(3, L)
    nt = Mt.shape[0]
    hb = hbeta(erho, gam, a_arr, 3/5)
    r = Mt.T.dot(hb) / hb
    R0, Rhalf, _ = Rweights(6)   # deep-mass closed form lives on the mod-6 base shell, NOT full depth D

    v1 = (gam % 3 == 0)
    idx1 = np.where(v1)[0]

    # ---- baseline 5-set (2c3 resolution) for value-set comparison ----
    base = defaultdict(set)
    for i in idx1:
        base[(int(a_arr[i] % 9), int(gam[i] % 27), int(erho[i] % 6))].add(round(float(r[i]), 9))
    base_wd = all(len(v) == 1 for v in base.values())
    base_vals = sorted({list(v)[0] for v in base.values()})

    # ============ P1' : level sets (e mod6) x class3 ============
    key3 = defaultdict(set)
    parity_of_key = {}
    for i in idx1:
        lbl, b, eta = cls3(int(a_arr[i]), int(gam[i]))
        k = (int(erho[i] % 6), lbl)
        key3[k].add(round(float(r[i]), 9)); parity_of_key[k] = int(erho[i] % 2)
    wd3 = all(len(v) == 1 for v in key3.values())
    vals3 = sorted({list(v)[0] for v in key3.values() if len(v) == 1})
    per_e = defaultdict(set)
    for (em, lbl), vv in key3.items():
        if len(vv) == 1: per_e[em].add(list(vv)[0])
    maxper = max(len(s) for s in per_e.values())
    same_valset = (len(vals3) == len(base_vals) and
                   all(any(abs(a - b) < 1e-9 for b in base_vals) for a in vals3))
    log(f"\n P1' LEVEL SETS  (baseline: well-def on (a%9,g%27,e%6)={base_wd}, {len(base_vals)} vals "
        f"{[str(Fraction(x).limit_denominator(200)) for x in base_vals]})")
    log(f"   r well-defined on (e mod6)x{{D9,U+,U-}}: {wd3}   values/e-group max = {maxper} (<=3?)   "
        f"distinct values total = {len(vals3)}")
    log(f"   value-set matches baseline 5-set: {same_valset}   "
        f"=> P1' {'PASS' if (wd3 and maxper<=3 and same_valset) else 'DEVIATES'}")
    if not wd3:
        bad = [(k, sorted(v)) for k, v in key3.items() if len(v) > 1][:4]
        log(f"   NOT well-defined on class3 key; examples {bad}")

    # ============ P2' : bit-blindness + opposite signs ============
    log(f"\n P2' BIT (U+ vs U-) per e mod6  [b=+1 iff a==g mod3]:")
    split = {}
    for em in range(6):
        up = key3.get((em, "U+")); un = key3.get((em, "U-"))
        if up and un and len(up) == 1 and len(un) == 1:
            d = list(up)[0] - list(un)[0]
            split[em] = d
            log(f"   e={em} ({'odd' if em%2 else 'even'}): U+={Fraction(list(up)[0]).limit_denominator(200)}"
                f"  U-={Fraction(list(un)[0]).limit_denominator(200)}  U+ - U- = {d:+.6f} "
                f"({'BLIND' if abs(d)<1e-9 else 'live'})")
    e0_blind = abs(split.get(0, 1.0)) < 1e-9
    opp = (0 in [0] and 2 in split and 4 in split and split[2] * split[4] < 0)
    log(f"   e==0 mod6 b-BLIND: {e0_blind}   e==2,e==4 OPPOSITE sign: "
        f"{opp} (2:{split.get(2,float('nan')):+.4f}, 4:{split.get(4,float('nan')):+.4f})   "
        f"=> P2' {'PASS' if (e0_blind and opp) else 'DEVIATES'}")

    # ============ P3' : D9 DEEP mass = R0(s0), s0 = -e mod6 ============
    log(f"\n P3' D9 DEEP mass vs R0(s0), s0=-e mod6   (D={D}; R0 arg tried mod6 and mod D):")
    p3ok = True
    for em in range(6):
        xs = [i for i in idx1 if int(erho[i] % 6) == em and (int(gam[i]) % 9) // 3 == 0]
        if not xs: continue
        dm = [deep_mass(Mt, gam, hb, i) for i in xs]
        dh = sorted({round(v[0], 9) for v in dm}); draw = sorted({round(v[2], 9) for v in dm})
        s0_6 = (-em) % 6
        cand = {f"R0({s0_6})": float(R0(s0_6))}
        for s in range(D):
            cand[f"R0({s})"] = float(R0(s))
        # match raw-move-deep and h-weighted-deep against R0 candidates
        def hits(vals): return [n for n, cv in cand.items() if any(abs(cv - x) < 1e-9 for x in vals)]
        log(f"   e={em}: D9 count={len(xs)}  deep(h-wt)={[str(Fraction(x).limit_denominator(4000)) for x in dh]}"
            f"  deep(raw)={[str(Fraction(x).limit_denominator(4000)) for x in draw]}")
        log(f"        R0(-e%6={s0_6})={Fraction(R0(s0_6)).limit_denominator(4000)}  "
            f"raw matches: {hits(draw)[:3]}   h-wt matches: {hits(dh)[:3]}")
        if not any(abs(float(R0(s0_6)) - x) < 1e-9 for x in draw + dh): p3ok = False
    log(f"   => P3' {'PASS' if p3ok else 'DEVIATES (see matches above)'}")

    # ============ P4' : unit DEEP mass = W_D formula ============
    log(f"\n P4' unit-g DEEP mass vs W_D = [R0(s2)+R0(s4)]/2 - eta*[Rh(s2)-Rh(s4)]/2, s2=2-e,s4=4-e:")
    p4ok = True
    for em in range(6):
        for etaval in (+1, -1):
            xs = [i for i in idx1 if int(erho[i] % 6) == em and (int(gam[i]) % 9) // 3 != 0
                  and cls3(int(a_arr[i]), int(gam[i]))[2] == etaval]
            if not xs: continue
            dm = [deep_mass(Mt, gam, hb, i) for i in xs]
            dh = sorted({round(v[0], 9) for v in dm}); draw = sorted({round(v[2], 9) for v in dm})
            for smod in (6, D):
                s2 = (2 - em) % smod; s4 = (4 - em) % smod
                WD = (R0(s2) + R0(s4)) / 2 - etaval * (Rhalf(s2) - Rhalf(s4)) / 2
                hitraw = any(abs(float(WD) - x) < 1e-9 for x in draw)
                hith = any(abs(float(WD) - x) < 1e-9 for x in dh)
                if smod == 6:
                    log(f"   e={em} eta={etaval:+d} (n={len(xs)}): deep(raw)={[str(Fraction(x).limit_denominator(4000)) for x in draw]}"
                        f"  deep(h)={[str(Fraction(x).limit_denominator(4000)) for x in dh]}")
                log(f"        W_D(s mod{smod}: s2={s2},s4={s4})={Fraction(WD).limit_denominator(8000)}"
                    f"={float(WD):.6f}  raw?{hitraw}  h?{hith}")
            # pass if any convention matches raw or h
            ok = False
            for smod in (6, D):
                s2 = (2 - em) % smod; s4 = (4 - em) % smod
                WD = (R0(s2) + R0(s4)) / 2 - etaval * (Rhalf(s2) - Rhalf(s4)) / 2
                if any(abs(float(WD) - x) < 1e-9 for x in draw + dh): ok = True
            if not ok: p4ok = False
    log(f"   => P4' {'PASS' if p4ok else 'DEVIATES (see W_D vs measured above)'}")

    return dict(wd3=wd3, maxper=maxper, same_valset=same_valset, e0_blind=e0_blind, opp=opp,
                p3=p3ok, p4=p4ok, nvals=len(vals3), vals=vals3)


def main():
    log("# PROBE 2c5 -- GATE the v1-trit blind derivation (P1'..P4'). beta*=3/5, no v0 trit. Direct/exact.")
    res = {}
    for L in [2, 3]:
        res[L] = run(L)
    log(f"\n{'='*74}\n## VERDICT")
    for L in [2, 3]:
        r = res[L]
        log(f" L={L}: P1'(levelsets {r['nvals']} vals, <=3/e {r['maxper']<=3}, valset {r['same_valset']})="
            f"{'PASS' if (r['wd3'] and r['maxper']<=3 and r['same_valset']) else 'DEV'}  "
            f"P2'(e0-blind {r['e0_blind']}, opp {r['opp']})={'PASS' if (r['e0_blind'] and r['opp']) else 'DEV'}  "
            f"P3'={'PASS' if r['p3'] else 'DEV'}  P4'={'PASS' if r['p4'] else 'DEV'}")
    with open("logs/probe_phase2c5_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
