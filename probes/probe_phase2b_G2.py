"""
MICRO-PROBE G2 -- gate Wilson's two fine-structure lemmas U0 + D2 (per-state / per-edge, exact).
INSTRUMENT LAW: direct/exact at q=3. No proof authored, no rate fit. Combinatorial + machine-precision.

Tower state (a,b,gamma), gamma!=0. Move (da,db): ap=a*2^-da, bp=b*2^-db, T=ap-bp mod 3^L,
gate (gamma+T)%3==0, carry gamma'=((gamma+T)//3)%3^L. e'=(e_rho+da-db) mod D, t(e')=v3(2^e'-1).

LEMMA U0 (v3(gamma)=0 sources are EXACTLY mean-field, PER STATE):
  PRE-REG: per-state v'-split from EVERY v=0 tower state = geometric {2/3,2/9,1/9} (L=3) to machine
  precision. Mechanism claim: surviving new-carry gamma' is uniform on Z/3^{L-1}.

LEMMA D2 (v3(gamma)>=1 fine law; gamma=3g):
  PRE-REG (a): t>=2 targets (e'==0 mod 6) preserve g mod 3 EXACTLY (digit-shift, deterministic).
  PRE-REG (b): t=1 targets (e'==2,4 mod 6) from g==0 mod 3 states land v'=0 with weight 1.
  Bonus: t=1 from g==1,2 mod 3 states split [1/2,1/2] between v'=0 and v'>=1.
"""
import numpy as np
from collections import defaultdict

from probe_phase2a_q2b_q6 import subgroup

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def v3int(g, cap=None):
    if g == 0: return 10 ** 9 if cap is None else cap
    v = 0
    while g % 3 == 0: g //= 3; v += 1
    return v


def run(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    inv2 = pow(2, -1, qL); pw = [pow(inv2, d, qL) for d in range(D + 1)]
    w = np.array([lam ** d for d in range(1, D + 1)]); w = w / w.sum()
    # LTE class of each target phase e'
    def t_of(ep):
        val = 2 ** ep - 1
        if val == 0: return 999                       # e'=0 => T=0, digit-shift (>=2 bucket)
        v = 0
        while val % 3 == 0: val //= 3; v += 1
        return v
    # geometric mean-field target on Z/3^{L-1}: P(v'=j)=2*3^-(j+1), tail at j=L-1
    tgt = {j: 2 * 3 ** (-(j + 1)) for j in range(L - 1)}
    tgt[L - 1] = 3 ** (-(L - 1))

    # U0 accumulators
    u0_maxdev = 0.0; u0_states = 0; u0_worst = None
    u0_gp_uniform_maxspread = 0.0; u0_support = set()
    u0_count_maxdev = 0.0; u0_count_spread = 0       # count-based mean-field + count-uniformity of gamma'
    # D2 accumulators
    d2a_total = 0; d2a_viol = 0
    d2b_total = 0; d2b_viol = 0
    d2_t1_g12_v0 = 0.0; d2_t1_g12_vpos = 0.0     # weighted split for g==1,2 mod3

    for a in sub:
        ainv = pow(a, -1, qL)
        ap = [(a * pw[d]) % qL for d in range(D + 1)]
        for b in sub:
            er = dl[(b * ainv) % qL]
            bp = [(b * pw[d]) % qL for d in range(D + 1)]
            for gam in range(1, qL):                  # tower source carry
                vg = v3int(gam)
                dist = defaultdict(float); surv = 0.0
                distc = defaultdict(int); survc = 0   # unweighted COUNT v'-split
                gpw = defaultdict(float)              # gamma' -> weight (for U0 uniformity)
                gpc = defaultdict(int)                # gamma' -> edge COUNT
                g_src = gam // 3 if gam % 3 == 0 else None
                for da in range(1, D + 1):
                    for db in range(1, D + 1):
                        T = (ap[da] - bp[db]) % qL
                        if (gam + T) % q == 0:
                            gp = ((gam + T) // q) % qL
                            wt = w[da - 1] * w[db - 1]
                            vp = v3int(gp, cap=L - 1)
                            vpj = min(vp, L - 1)
                            surv += wt; dist[vpj] += wt; gpw[gp] += wt
                            survc += 1; distc[vpj] += 1; gpc[gp] += 1
                            if vg >= 1:
                                ep = (er + da - db) % D
                                tt = t_of(ep)
                                if tt >= 2:            # digit-shift channel
                                    d2a_total += 1
                                    if gp % 3 != g_src % 3: d2a_viol += 1
                                elif tt == 1:
                                    if g_src % 3 == 0:
                                        d2b_total += 1
                                        if vpj != 0: d2b_viol += 1
                                    else:
                                        if vpj == 0: d2_t1_g12_v0 += wt
                                        else: d2_t1_g12_vpos += wt
                if vg == 0:
                    u0_states += 1
                    dev = max(abs(dist[j] / surv - tgt[j]) for j in range(L))
                    if dev > u0_maxdev:
                        u0_maxdev = dev
                        u0_worst = (a, b, gam, {j: dist[j] / surv for j in range(L)})
                    # uniformity of gamma': equal WEIGHT per distinct gamma'?
                    vals = np.array(list(gpw.values()))
                    u0_gp_uniform_maxspread = max(u0_gp_uniform_maxspread, vals.max() - vals.min())
                    u0_support.add(len(gpw))
                    # COUNT version: count v'-split vs mean-field + count-uniformity of gamma'
                    cdev = max(abs(distc[j] / survc - tgt[j]) for j in range(L))
                    u0_count_maxdev = max(u0_count_maxdev, cdev)
                    cvals = np.array(list(gpc.values()))
                    u0_count_spread = max(u0_count_spread, int(cvals.max() - cvals.min()))

    log(f"\n{'='*74}\n## G2  q=3 L={L}   (tower dim {D*D*(qL-1)}; mean-field target {{{','.join(f'{tgt[j]:.4f}' for j in range(L))}}})")
    log(f"\n   LEMMA U0 (v=0 sources exactly mean-field, PER STATE):")
    log(f"      [MASS, w-weighted] per-state v'-split vs geometric: MAX dev over {u0_states} v=0 states = {u0_maxdev:.3e}"
        f"   ({'CONFIRMED <=1e-12' if u0_maxdev <= 1e-12 else 'REFUTED (per-state mass NOT mean-field)'})")
    log(f"      [MASS] gamma' support = {sorted(u0_support)} (=3^(L-1)={3**(L-1)}), within-state weight spread = "
        f"{u0_gp_uniform_maxspread:.3e}   ({'uniform' if u0_gp_uniform_maxspread <= 1e-15 else 'NOT uniform (geometric w-weights)'})")
    log(f"      [COUNT, unweighted] per-state v'-split vs geometric: MAX dev = {u0_count_maxdev:.3e}"
        f"   ({'CONFIRMED <=1e-12' if u0_count_maxdev <= 1e-12 else 'dev'})")
    log(f"      [COUNT] gamma' edge-count spread within state = {u0_count_spread}"
        f"   ({'UNIFORM in COUNT: the carry MAP equidistributes on Z/3^(L-1) (Wilson proof correct for the MAP)' if u0_count_spread == 0 else 'NOT count-uniform'})")
    log(f"\n   LEMMA D2 (v>=1 fine law, gamma=3g):")
    log(f"      (a) t>=2 (e'==0 mod6) preserve g mod 3: {d2a_total-d2a_viol}/{d2a_total} edges"
        f"   ({'CONFIRMED (deterministic)' if d2a_viol == 0 and d2a_total > 0 else ('VIOLATIONS='+str(d2a_viol)) if d2a_total>0 else 'no such edges at this L'})")
    if d2b_total > 0:
        log(f"      (b) t=1 from g==0 mod3 states land v'=0: {d2b_total-d2b_viol}/{d2b_total} edges"
            f"   ({'CONFIRMED (weight 1)' if d2b_viol == 0 else 'VIOLATIONS='+str(d2b_viol)})")
    else:
        log(f"      (b) t=1 from g==0 mod3: no g==0 tower states at L={L} (needs v3(gamma)>=2) -- untestable here")
    tot12 = d2_t1_g12_v0 + d2_t1_g12_vpos
    if tot12 > 0:
        log(f"      bonus: t=1 from g==1,2 mod3 split [v'=0, v'>=1] = "
            f"[{d2_t1_g12_v0/tot12:.6f}, {d2_t1_g12_vpos/tot12:.6f}]   "
            f"({'CONFIRMED [1/2,1/2]' if abs(d2_t1_g12_v0/tot12-0.5)<1e-12 else 'DEVIATION from 1/2'})")
    if u0_worst is not None and u0_maxdev > 1e-12:
        log(f"      [U0 worst state (a,b,gamma)={u0_worst[:3]} split={u0_worst[3]}]")
    return u0_maxdev, d2a_viol, d2b_viol, d2b_total


def main():
    log("# MICRO-PROBE G2 -- gate U0 (mean-field per state) + D2 (v>=1 digit-tape law). Exact at q=3.")
    r2 = run(2); r3 = run(3)
    log(f"\n{'='*74}\n## SUMMARY")
    log(f"   U0 per-state mean-field: L=2 maxdev {r2[0]:.1e}, L=3 maxdev {r3[0]:.1e}  "
        f"({'BOTH CONFIRMED' if max(r2[0], r3[0]) <= 1e-12 else 'see deviations'})")
    log(f"   D2(a) g-mod-3 preservation: L=2 viol {r2[1]}, L=3 viol {r3[1]}  "
        f"({'BOTH deterministic' if r2[1]==0 and r3[1]==0 else 'violations'})")
    log(f"   D2(b) t=1/g==0 -> v'=0: L=3 {r3[3]-r3[2]}/{r3[3]} edges  "
        f"({'CONFIRMED' if r3[2]==0 and r3[3]>0 else 'see line'})")
    with open("logs/probe_phase2b_G2_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
