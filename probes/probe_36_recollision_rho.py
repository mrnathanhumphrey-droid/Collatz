"""
PROBE 36 -- the re-collision conditioning factor rho(q), read off the exact structure.

REQUEST (user): among off-diagonal pairs that collide at level 1 (shift d|m), what is the
conditional P(re-collide at level 2), vs the unconditional level-1 rate? Define
   rho(q) := P(collide L2 | collided L1) / P(collide L1)      [ = r_2 / r_1 ]
Also report the carry gamma distribution conditioned on level-1 collision (equidistributed
vs biased is the whole question), and the per-level rates r_1,r_2,r_3.

METHOD (exact, no fit): enumerate address pairs (a=(a_1..a_k)), value = sum_m q^{m-1} 2^{-S_m}
mod q^k (S_m suffix sums), weight w(a)=prod 2^{-a_i}. Group by val mod q^m. Off-diagonal
collision mass at level m: offmass(m) = sum_{res} (sum w)^2 - sum w^2  (distinct-address).
   r_1 = offmass(1)/P_offdiag,  r_2 = offmass(2)/offmass(1),  r_3 = offmass(3)/offmass(2)
   rho = r_2 / r_1.  Compare r_1 to P_coll(d) = (2/3)/(2^d - 1).
gamma dist: within each L1-class c (val=c+q*h mod q^2), gamma_1 = h_i - h_j mod q for a
colliding pair; the gamma distribution = the weighted difference-distribution of h.

PRE-REGISTERED (user's, stated to lose): q=3 (d=2) rho ~ 1 (no penalty, parity self-sim);
d>=3 rho < 1 (carry ~equidistributed -> fresh ~1/(2^d-1) penalty). gamma flat => contraction;
gamma biased => fat rho. Reported raw, no fit.

Secondary: q=1093 (s=2) -- level-2 conditional should look like other primes' level-1 (R35
onset shift). k=2 only (deep shifts underflow; report with caution).

Foreground; exact grouping (O(V^k)); no windows/extrapolation.
"""
import numpy as np
from itertools import product
from collections import defaultdict
from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def analyze(q, k, V, do_gamma=True):
    N = q ** k
    inv2 = pow(2, -1, N)
    qp = [q ** m for m in range(k + 1)]
    inv2pow = [1] * (k * V + 2)
    for i in range(1, len(inv2pow)):
        inv2pow[i] = (inv2pow[i - 1] * inv2) % N
    Sm = [defaultdict(float) for _ in range(k + 1)]
    SSm = [defaultdict(float) for _ in range(k + 1)]
    Wtot = 0.0
    W2tot = 0.0
    Hc = defaultdict(lambda: defaultdict(float))   # Hc[c mod q][h]  (h = (val mod q^2 - c)//q)
    for addr in product(range(1, V + 1), repeat=k):
        s = 0
        val = 0
        for m in range(1, k + 1):
            s += addr[k - m]
            val = (val + qp[m - 1] * inv2pow[s]) % N
        w = 2.0 ** (-sum(addr))
        Wtot += w
        W2tot += w * w
        for m in range(1, k + 1):
            r = val % qp[m]
            Sm[m][r] += w
            SSm[m][r] += w * w
        if do_gamma and k >= 2:
            c = val % q
            h = (val % qp[2] - c) // q
            Hc[c][h] += w

    def offmass(m):
        return sum(sv * sv - SSm[m][rv] for rv, sv in Sm[m].items())

    Poff = Wtot * Wtot - W2tot
    om = {m: offmass(m) for m in range(1, k + 1)}
    r1 = om[1] / Poff if Poff else float('nan')
    r2 = om[2] / om[1] if om[1] > 0 else float('nan')
    r3 = (om[3] / om[2] if k >= 3 and om[2] > 0 else None)
    rho = r2 / r1 if r1 > 0 else float('nan')
    gdist = None
    if do_gamma:
        gd = defaultdict(float)
        diag0 = 0.0
        for c, H in Hc.items():
            items = list(H.items())
            for h, wh in items:
                diag0 += wh * wh
                for h2, wh2 in items:
                    gd[(h - h2) % q] += wh * wh2
        gd[0] -= diag0
        tot = sum(gd.values())
        gdist = {g: gd[g] / tot for g in range(q)} if tot > 0 else None
    return dict(r1=r1, r2=r2, r3=r3, rho=rho, om=om, Poff=Poff, gdist=gdist)


def main():
    log("# PROBE 36 -- re-collision conditioning rho(q) = r_2/r_1, + gamma distribution")
    log("")
    log(f"{'q':>5} {'d':>3} {'r_1(L1)':>10} {'P_coll(d)':>11} {'r_2(L2|L1)':>11} "
        f"{'r_3':>8} {'rho=r2/r1':>10}")
    results = {}
    for q in [3, 5, 7, 11, 13]:
        d = order_of_two(q)
        k = 3
        V = min(4 * d, 40)
        if V ** k > 3_000_000:
            V = int(3_000_000 ** (1 / k))
        res = analyze(q, k, V)
        results[q] = res
        Pcoll = (2 / 3) / (2 ** d - 1)
        r3s = f"{res['r3']:.5f}" if res['r3'] is not None else "-"
        log(f"{q:>5} {d:>3} {res['r1']:>10.6f} {Pcoll:>11.6f} {res['r2']:>11.6f} "
            f"{r3s:>8} {res['rho']:>10.5f}   (V={V},k={k})")
    log("")
    log("## GAMMA distribution | level-1 collision (flat => equidistributed => contraction)")
    for q in [3, 5, 7, 11, 13]:
        gd = results[q]['gdist']
        if gd:
            vals = [gd[g] for g in range(q)]
            unif = 1.0 / q
            mx, mn = max(vals), min(vals)
            log(f"   q={q}: gamma dist = {['%.4f' % v for v in vals]}  (uniform={unif:.4f}; "
                f"max/min={mx/mn if mn>0 else float('inf'):.2f})")
    log("")

    # ---- rho(q) summary + pre-registered check ----
    log("## rho(q) vs pre-registration (q=3 ~1 no-penalty; d>=3 <1):")
    for q in [3, 5, 7, 11, 13]:
        d = order_of_two(q)
        rho = results[q]['rho']
        tag = "~1 (NO penalty)" if abs(rho - 1) < 0.15 else ("<1 (penalty)" if rho < 1 else ">1 (FAT!)")
        log(f"   q={q} (d={d}): rho = {rho:.4f}  -> {tag}")
    log("")

    # ---- secondary: q=1093 (s=2) onset shift ----
    log("## SECONDARY -- q=1093 (s=2): level-2 conditional should mimic level-1 (R35 onset shift)")
    q = 1093
    d = order_of_two(q)
    V = 800
    try:
        res = analyze(q, 2, V, do_gamma=False)
        Pcoll = (2 / 3) / (2 ** d - 1)
        log(f"   q=1093 d={d} (s=2): r_1(L1)={res['r1']:.3e}  P_coll(d)={Pcoll:.3e}  "
            f"r_2(L2|L1)={res['r2']:.3e}  rho={res['rho']:.4f}  (V={V},k=2)")
        log(f"   offmass: L1={res['om'][1]:.3e}  L2={res['om'][2]:.3e}  "
            f"(L2<<L1 or L2=0 => onset shifted past depth 2, per R35 cross(2)=0)")
    except Exception as e:
        log(f"   q=1093: {e}")
    log("")
    log("## READ: rho<1 for d>=3 with ~flat gamma => carry equidistributes => re-collision")
    log("   pays a fresh penalty => contraction. rho~1 at q=3 => no penalty => divergence.")
    flush()


def flush():
    with open("result_36_recollision_rho_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
