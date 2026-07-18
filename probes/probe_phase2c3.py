"""
PROBE 2c3 -- rung-2 inputs + judge. Direct/exact, no eigensolves. (C held until Wilson sends tau.)
A) MOD-9 GATE TABLES: v0 sources, per e' mod 6: which u mod 9 pass, carry gamma' mod 3 as f(u mod9,
   gamma mod9); + R(s) mass weights; + the v>=1 analog (residual spread 1/7).
B) THE JUDGE: at beta*=3/5, the ratio r(x)=(M^T h_beta)/h_beta for every v0 tower state, keyed by
   (a mod 9, gamma mod 9, e mod 6). PRE-REG SHAPE: exactly 3 populations, r constant on the level sets
   (=> tau exists at this resolution). Deviations reported as deviations.
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


def build_tower(q, L, lam=0.5):
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
    return M[tw][:, tw].tocsc(), gam[tw], erho[tw], a_arr[tw], D, qL, dl, sub


def hbeta(erho, gam, a_arr, beta):
    nt = len(erho); par = erho % 2; v0 = (gam % 3 != 0)
    h = np.empty(nt)
    h[(par == 0) & v0] = 2/3; h[(par == 0) & ~v0] = 5/3
    h[(par == 1) & v0] = 5/6; h[(par == 1) & ~v0] = 4/3
    ov0 = (par == 1) & v0
    sig = np.zeros(nt)
    sig[ov0] = np.where((a_arr[ov0] + gam[ov0]) % 3 == 0, 1.0, -1.0)
    return h * (1 + beta * sig)


# ============================ B : THE JUDGE ============================
def partB(L):
    Mt, gam, erho, a_arr, D, qL, dl, sub = build_tower(3, L)
    nt = Mt.shape[0]
    hb = hbeta(erho, gam, a_arr, Fraction(3, 5).__float__())
    r = Mt.T.dot(hb) / hb
    v0 = (gam % 3 != 0)
    idx = np.where(v0)[0]
    keymap = defaultdict(set)
    for i in idx:
        key = (int(a_arr[i] % 9), int(gam[i] % 9), int(erho[i] % 6))
        keymap[key].add(round(float(r[i]), 9))
    well_defined = all(len(v) == 1 for v in keymap.values())
    vals = sorted({list(v)[0] for v in keymap.values()})
    log(f"\n   B) JUDGE L={L}: r at beta*=3/5 on v0 states, keyed by (a mod9, gamma mod9, e mod6)")
    log(f"      well-defined on (a mod9,gamma mod9,e mod6): {well_defined}   distinct r-values: {len(vals)} "
        f"= {[str(Fraction(x).limit_denominator(200)) for x in vals]}")
    if not well_defined:
        bad = [(k, sorted(v)) for k, v in keymap.items() if len(v) > 1][:5]
        log(f"      NOT well-defined (needs finer resolution). examples: {bad}")
    # level sets: value -> set of keys, sizes
    lset = defaultdict(list)
    for k, v in keymap.items():
        lset[list(v)[0]].append(k)
    for val in vals:
        keys = lset[val]
        log(f"      population r={Fraction(val).limit_denominator(200)} ({float(val):.6f}): {len(keys)} keys")
    # structural characterization: does r factor through a coarser key? (judge sharpening; tau is Wilson's)
    def factors_through(keyfn):
        km = defaultdict(set)
        for i in idx: km[keyfn(i)].add(round(float(r[i]), 9))
        return all(len(v) == 1 for v in km.values()), len({list(v)[0] for v in km.values() if len(v)==1})
    cands = [("(a+gamma mod9, e mod6)", lambda i: (int((a_arr[i]+gam[i]) % 9), int(erho[i] % 6))),
             ("(a+gamma mod9)",          lambda i: int((a_arr[i]+gam[i]) % 9)),
             ("(a*gamma mod9, e mod6)",  lambda i: (int((a_arr[i]*gam[i]) % 9), int(erho[i] % 6))),
             ("(a mod9, gamma mod9)",    lambda i: (int(a_arr[i] % 9), int(gam[i] % 9)))]
    log(f"      level-set factorization (judge sharpening):")
    for name, fn in cands:
        wd, nv = factors_through(fn)
        log(f"        r factors through {name}: {wd}" + (f" ({nv} values)" if wd else ""))
    # dump the full level-set table
    lines = [f"# 2c3-B judge: v0 ratio at beta*=3/5, q=3 L={L}. key=(a mod9,gamma mod9,e mod6) -> r (exact frac)",
             "# a%9\tgamma%9\te%6\tr"]
    for k in sorted(keymap):
        v = list(keymap[k])[0]
        lines.append(f"{k[0]}\t{k[1]}\t{k[2]}\t{Fraction(v).limit_denominator(2000)}")
    with open(f"outputs/rung2_judge_v0_q3_L{L}.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    # also v>=1 residual shape (bonus)
    v1 = (gam % 3 == 0)
    idx1 = np.where(v1)[0]
    km1 = defaultdict(set)
    for i in idx1:
        key = (int(a_arr[i] % 9), int(gam[i] % 27), int(erho[i] % 6))
        km1[key].add(round(float(r[i]), 9))
    wd1 = all(len(v) == 1 for v in km1.values())
    vals1 = sorted({list(v)[0] for v in km1.values()})
    log(f"      [v>=1 residual] keyed by (a mod9, gamma mod27, e mod6): well-defined {wd1}, "
        f"{len(vals1)} values {[str(Fraction(x).limit_denominator(200)) for x in vals1]}")
    return well_defined, vals


# ============================ A : MOD-9 GATE TABLES ============================
def partA(L):
    qL = 3 ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    # units mod 9 and their mass under R_0(s): the shift s = e'-e; here tabulate the gate/carry at mod 9.
    log(f"\n   A) MOD-9 GATE/CARRY TABLE (v0 sources), q=3 L={L}: per e' mod6, gate (mod3) + carry gamma' mod3")
    log(f"      units mod 9 = {sorted(set(u % 9 for u in subgroup(2 % 9, 9)))}; e' mod6 sets c=1-2^{{e'}} mod9")
    lines = [f"# 2c3-A mod-9 gate/carry, q=3 (v0). rows: (e' mod6, u mod9, gamma mod9) -> gate_pass, gamma' mod3",
             "# e'mod6\tc=1-2^e' mod9\tu%9\tgamma%9\tgate_pass\tgammap%3"]
    for epm in range(6):
        c = (1 - pow(2, epm, 9)) % 9
        passing_u = set()
        for u9 in sorted(set(u % 9 for u in subgroup(2 % 9, 9))):
            for g9 in range(9):
                T = (u9 * c) % 9
                gp = (g9 + T) % 3 == 0
                if gp:
                    passing_u.add(u9)
                    gpm3 = ((g9 + u9 * c) // 3) % 3
                    lines.append(f"{epm}\t{c}\t{u9}\t{g9}\t1\t{gpm3}")
        log(f"      e' mod6={epm}: c={c}; u mod9 that can pass = {sorted(passing_u)} "
            f"(LTE t(e')={ '0(odd)' if epm%2 else ('2' if epm in (0,) else '1') })")
    with open(f"outputs/rung2_gate_mod9_q3.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    # R_0(s) mass weights (the move-weight per shift), exact
    Z = 2 ** D - 1
    wf = [Fraction(2 ** (D - m), Z) for m in range(1, D + 1)]
    R0 = [sum(wf[m - 1] * wf[((m - s - 1) % D)] for m in range(1, D + 1)) for s in range(D)]
    log(f"      R_0(s) mass weights s=0..{min(D,6)-1}: " + "  ".join(f"{s}:{R0[s]}" for s in range(min(D, 6))))
    log(f"      dumped outputs/rung2_gate_mod9_q3.tsv (+ v>=1 one level up in judge dump)")


def main():
    log("# PROBE 2c3 -- rung-2 inputs (A) + judge (B). Direct/exact. C held for Wilson's tau.")
    partA(2); partA(3)
    partB(2); partB(3)
    with open("logs/probe_phase2c3_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
