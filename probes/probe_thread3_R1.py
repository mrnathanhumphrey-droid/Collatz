"""
PROBE R1 -- THE RENEWAL/LOSS SPINE (thread 3). ONE frozen instrument: build_M_gen (the pair operator).
q=7 (and q=5) control columns woven in. Guards: NO fitting, raw dumps are the deliverable, exact rationals
where feasible, deviations as deviations.

Renewal variable (pinned from R7/R16): Y_k = 3^k * ||pi_k||^2  (the "3" is UNIVERSAL in q, R5/R8), with
  a_k = ||pi_k||^2 = 1^T M^k v0,  v0 = delta(1,1,0)  [frozen M; exact for k<=L, operator dynamics beyond].
  c_k = Y_k - Y_{k-1}  (per-scale loss; R19's increment).  X_k (renewal accumulation) = Y_k = sum_{j<=k} c_j.
PRE-REG (ALGEBRAIC, banked): q=3 -> c_k FLAT at 7/15 (slope of X_k -> 7/15); q=5,7 -> c_k DECAYS, ratio 3/q.

R1-C mass sequence (Wilson's Prediction R shape, sealed): m_n = ||M_tower^n . uniform||_1, n=1..200.
  RAW DUMP ONLY -- no exponent extraction, no regression, no log-log by code.
"""
import numpy as np, scipy.sparse as sp, os
from fractions import Fraction as Fr
from collections import defaultdict
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

np.set_printoptions(linewidth=160, suppress=True)
SEVEN_15 = 7 / 15


def full_M(q, L, lam=0.5):
    qL = q ** L; sub = subgroup(2, qL); ordn = len(sub)
    raw = [lam ** d for d in range(1, ordn + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    return M.tocsr(), idx, n


def tower_block(M, idx):
    states = [None] * len(idx)
    for s, i in idx.items():
        states[i] = s
    tw = np.array([i for i in range(len(states)) if states[i][2] != 0])
    return M[tw][:, tw].tocsr()


def Y_c_sequence(M, idx, K=30):
    """a_k = 1^T M^k v0 (float, frozen); Y_k = 3^k a_k; c_k = Y_k - Y_{k-1}. Y_0 = 1."""
    n = M.shape[0]; v = np.zeros(n); v[idx[(1, 1, 0)]] = 1.0
    Y = [1.0]; c = []
    for k in range(1, K + 1):
        v = M.dot(v)
        Yk = (3.0 ** k) * v.sum()
        c.append(Yk - Y[-1]); Y.append(Yk)
    return Y[1:], c                                              # Y[k-1]=Y_k, c[k-1]=c_k, k=1..K


def exact_a(q, L, K=3, lam=Fr(1, 2)):
    """Exact-rational a_k = ||pi_k||^2 for small k on the frozen instrument (cross-era weld)."""
    qL = q ** L; sub = subgroup(2, qL); ordn = len(sub); inv = pow(2, -1, qL)
    raw = [lam ** d for d in range(1, ordn + 1)]; Z = sum(raw); w = [r / Z for r in raw]
    mult = [(pow(inv, d, qL), w[d - 1]) for d in range(1, ordn + 1)]
    v = defaultdict(Fr); v[(1, 1, 0)] = Fr(1); a = []
    for k in range(K):
        nv = defaultdict(Fr)
        for (a_, b_, g), val in v.items():
            for (ga, wa) in mult:
                ap = (a_ * ga) % qL
                for (gb, wb) in mult:
                    bp = (b_ * gb) % qL
                    T = (ap - bp) % qL
                    if (g + T) % q == 0:
                        gp = ((g + T) // q) % qL
                        nv[(ap, bp, gp)] += val * wa * wb
        v = nv; a.append(sum(v.values()))
    return a                                                    # a[k-1] = a_k exact


def mass_seq(Mt, N=200):
    """m_n = ||Mt^n . uniform||_1 = sum(Mt^n 1) (Mt nonneg). RAW dump."""
    n = Mt.shape[0]; v = np.ones(n); m = []
    for _ in range(N):
        v = Mt.dot(v); m.append(float(v.sum()))
    return m


# ============================ R1-A / R1-B ============================
def renewal_spine():
    print(f"\n{'='*90}\n## R1-A / R1-B  RENEWAL ACCUMULATION X_k=Y_k and PER-SCALE LEDGER c_k  (frozen instrument)")
    cols = {}
    plan = [(3, 3), (5, 2), (7, 2)]                             # (q, L); q=7 L=3 = 7.4M states, infeasible -> L=2
    for q, L in plan:
        M, idx, n = full_M(q, L)
        Y, c = Y_c_sequence(M, idx, K=30)
        cols[q] = (Y, c, L, n)
        print(f"\n   q={q} (L={L}, full-M dim {n}):")
        print(f"     k :   c_k (per-scale)      X_k=Y_k (accum)    X_k-(7/15)k [q3]" if q == 3 else
              f"     k :   c_k (per-scale)      X_k=Y_k (accum)    c_k/c_(k-1) [ratio->3/q={3/q:.4f}]")
        for k in range(1, 31):
            if q == 3:
                print(f"     {k:2d}: {c[k-1]:+.10f}   {Y[k-1]:+.8f}   {Y[k-1]-SEVEN_15*k:+.6f}")
            else:
                rat = (c[k-1]/c[k-2]) if k > 1 and c[k-2] != 0 else float('nan')
                print(f"     {k:2d}: {c[k-1]:+.10f}   {Y[k-1]:+.8f}   {rat:+.6f}")
    # slope readout (q=3): last few c_k vs 7/15 (NO fit -- raw levels)
    Y3, c3, _, _ = cols[3]
    print(f"\n   >> q=3 PRE-REG (7/15={SEVEN_15:.10f}): c_k tail levels "
          f"c_20={c3[19]:.8f}, c_25={c3[24]:.8f}, c_30={c3[29]:.8f} "
          f"(dev from 7/15: {c3[29]-SEVEN_15:+.2e}) [frozen L=3; deviation = L-truncation, not a fit]")
    for q in (5, 7):
        _, cq, _, _ = cols[q]
        print(f"   >> q={q} PRE-REG (ratio 3/q={3/q:.6f}): c_k ratios tail "
              f"c9/c8={cq[8]/cq[7]:.6f}, c12/c11={cq[11]/cq[10]:.6f} [decays -> geometric saturation of X_k]")
    return cols


# ============================ R1-C ============================
def mass_dump():
    print(f"\n{'='*90}\n## R1-C  MASS SEQUENCE m_n = ||M_tower^n . uniform||_1, n=1..200  (RAW DUMP; shape sealed for pen)")
    out = {}
    # q=3 L=3 tower (dense-built, sparse iterate)
    M3, idx3, _ = full_M(3, 3); Mt3 = tower_block(M3, idx3)
    out["q3_L3_tower"] = mass_seq(Mt3, 200)
    # q=3 L=3 FULL operator (kinematic + tower; two-sector interference version)
    out["q3_L3_full"] = mass_seq(M3, 200)
    # q=3 L=4 tower via cached SpMV
    cache = os.path.expanduser("~/j2_L4_Mt.npz")
    if os.path.exists(cache):
        Mt4 = sp.load_npz(cache).tocsr()
        out["q3_L4_tower"] = mass_seq(Mt4, 200)
        print(f"   q3 L4 tower: cached {Mt4.shape[0]} states")
    else:
        print("   q3 L4 tower: cache ~/j2_L4_Mt.npz absent -- SKIP (report)")
    # q=7 control: L=3 is 7.4M states (infeasible CPU) -> L=2 control (gapped band; no slow ringing expected)
    M7, idx7, _ = full_M(7, 2); Mt7 = tower_block(M7, idx7)
    out["q7_L2_tower"] = mass_seq(Mt7, 200)
    print(f"   [!] q=7 control at L=2 (not L=3): L=3 q=7 = 7.4M states, build infeasible on CPU. Deviation noted.")

    # print heads so the ringing PERIOD is visible in the raw dump (no fit, no extraction)
    for key, m in out.items():
        print(f"\n   --- {key}  (m_1..m_20, then m_50,100,150,200) ---")
        print("     " + "  ".join(f"{m[i]:.6e}" for i in range(min(20, len(m)))))
        idxs = [49, 99, 149, 199]
        print("     n=50,100,150,200: " + "  ".join(f"{m[i]:.6e}" for i in idxs if i < len(m)))
    return out


# ============================ R1-D ============================
def bookkeeping(cols):
    print(f"\n{'='*90}\n## R1-D  BOOKKEEPING RIDERS")
    # (i) exact c_k at q=3 vs Era-6 (S1=2/3, S2=10/21)
    a = exact_a(3, 3, K=3)
    Y = [Fr(1)] + [Fr(3) ** k * a[k-1] for k in range(1, 4)]
    cex = [Y[k] - Y[k-1] for k in range(1, 4)]
    print(f"   (i) EXACT c_k (q=3, frozen build): "
          f"c_1={cex[0]} ({float(cex[0]):.8f}), c_2={cex[1]} ({float(cex[1]):.8f}), c_3={cex[2]} ({float(cex[2]):.8f})")
    print(f"       Era-6 R70 sequence: S_1=2/3, S_2=10/21.  "
          f"c_1==2/3: {cex[0]==Fr(2,3)}  |  c_2==10/21: {cex[1]==Fr(10,21)}  "
          f"({'CROSS-ERA WELD CONFIRMED' if cex[0]==Fr(2,3) and cex[1]==Fr(10,21) else 'DEVIATION -- report'})")
    print(f"       (float c_k from the L=3 iteration matches: c_1={cols[3][1][0]:.8f}, c_2={cols[3][1][1]:.8f})")
    # (ii) 7/45 definition transcribed (frozen)
    print(f"   (ii) 7/45 PLANCKerel normalization (frozen def, c_seven_forty_fifth_derivation.py / Paper 5):")
    print(f"        ||d_(k+1)||^2  =  Sum_r' pi_(k+1)(r')^2  -  (1/3) Sum_r pi_k(r)^2   [R74 identity]")
    print(f"        ||d_(k+1)||^2  ~  c * (1/3)^k   with   c = 7/45.")
    print(f"        BRIDGE (zero weight; pen derives): 3c = 3*(7/45) = 7/15 = the flat per-scale level c_k.")
    print(f"        check: 3*(7/45) = {float(3*Fr(7,45)):.10f}  vs  7/15 = {SEVEN_15:.10f}  "
          f"({'EXACT' if 3*Fr(7,45)==Fr(7,15) else 'DEV'})")


def main():
    print("# PROBE R1 -- THE RENEWAL/LOSS SPINE (thread 3). One frozen instrument. No fitting; raw dumps.")
    cols = renewal_spine()
    massdata = mass_dump()
    bookkeeping(cols)
    # dump TSVs
    with open("outputs/thread3_renewal_ledger.tsv", "w", encoding="utf-8") as f:
        f.write("# R1-A/B renewal spine: c_k (per-scale) and X_k=Y_k (accum), frozen instrument\n")
        f.write("k\tq3_c_k\tq3_X_k\tq5_c_k\tq5_X_k\tq7_c_k\tq7_X_k\n")
        for k in range(1, 31):
            row = [str(k)]
            for q in (3, 5, 7):
                Y, c, _, _ = cols[q]
                row += [f"{c[k-1]:.12e}", f"{Y[k-1]:.12e}"]
            f.write("\t".join(row) + "\n")
    with open("outputs/thread3_mass_sequences.tsv", "w", encoding="utf-8") as f:
        keys = list(massdata.keys())
        f.write("# R1-C mass sequences m_n = ||M_tower^n . uniform||_1 (RAW; no fit/extraction)\n")
        f.write("n\t" + "\t".join(keys) + "\n")
        N = max(len(v) for v in massdata.values())
        for i in range(N):
            f.write(f"{i+1}\t" + "\t".join(f"{massdata[k][i]:.12e}" if i < len(massdata[k]) else "" for k in keys) + "\n")
    print("\n   [dump] outputs/thread3_renewal_ledger.tsv + outputs/thread3_mass_sequences.tsv")


if __name__ == "__main__":
    main()
