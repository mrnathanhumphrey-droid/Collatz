"""
Q6-B settled with a RELIABLE collision detector (ESPRIT failed on the q=3 Jordan).

DETECTOR (Jordan-aware, no exponential fitting): a_k = 1^T M^k v0. lambda_1 = Perron (direct eigs,
robust). b_k = a_k / lambda_1^k. Jordan/COLLISION at the top (r=1) <=> a_k ~ (A+Bk)lambda_1^k <=>
b_k grows LINEARLY (slope B != 0). GAPPED (r<1) <=> a_k ~ A lambda_1^k + C mu^k, |mu|<lambda_1 <=>
b_k -> const (B=0). Read B = late-k linear slope of b_k; also Delta b_k tail (->const vs ->0).

VALIDATE on knowns: M(3,2,1/2) collision (B~7/15); M(7,2,1/2) gapped; M(7,-1,1/2) gapped.
Q6-B: M(3,2,lam), lam in {1/4,1/3,2/5,1/2,3/5,0.7}. Is cross-Perron collision lam-BLIND (B>0 all lam)
   or does it unglue at some lam (B->0)? Settles whether the WEIGHT base is inessential at d=2.
"""
import numpy as np
import scipy.sparse.linalg as spla

from probe_phase2a_q2b_q6 import build_M_gen, subgroup
from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    print(m, flush=True); LOG.append(str(m))


def mass(M, idx, K=30):
    n = M.shape[0]; v = np.zeros(n); v[idx[(1, 1, 0)]] = 1.0
    out = []
    for _ in range(K):
        v = M.dot(v); out.append(v.sum())
    return np.array(out)


def detect(M, idx, K=30):
    lam1 = abs(spla.eigs(M, k=1, which='LM', return_eigenvectors=False)[0])
    a = mass(M, idx, K)
    ks = np.arange(1, K + 1)
    b = a / lam1 ** ks                                   # b_k, k=1..K
    db = np.diff(b)
    # late-k linear slope of b_k (last third)
    lo = 2 * K // 3
    B = np.polyfit(ks[lo:], b[lo:], 1)[0]
    # Delta b tail behaviour: ratio of |db| over last 6 (->1 collision, ->0 gapped)
    tail = np.abs(db[-6:])
    tail_ratio = tail[-1] / tail[0] if tail[0] > 1e-300 else float('nan')
    return lam1, B, b, db, tail_ratio


def classify(B, tail_ratio):
    # collision: b_k grows (B not ~0) and db tail stays flat (ratio ~ O(1))
    if abs(B) > 5e-3 and tail_ratio > 0.5:
        return "COLLISION (r=1, Jordan)"
    if abs(B) < 5e-3 and tail_ratio < 0.5:
        return "GAPPED (r<1)"
    return f"AMBIGUOUS (B={B:.4f}, tail_ratio={tail_ratio:.3f})"


def run(q, gen, L, lam, label):
    ordsub = len(subgroup(gen % (q ** L), q ** L))
    raw = [lam ** d for d in range(1, ordsub + 1)]
    M, idx, n = build_M_gen(q, L, gen, raw)
    lam1, B, b, db, tr = detect(M, idx)
    verdict = classify(B, tr)
    log(f"   {label:<22} dim={n:<6} lam1={lam1:.5f}  B(slope)={B:>8.4f}  db_tail_ratio={tr:>6.3f}  -> {verdict}")
    return B, verdict


def main():
    log("# Q6-B RELIABLE COLLISION DETECTOR (Jordan-aware, no ESPRIT). b_k=a_k/lam1^k: linear=collision.")
    log("")
    log("## VALIDATION on knowns (detector must reproduce banked collision/gap):")
    run(3, 2, 2, 0.5, "M(3,2,1/2) [banked r=1]")
    run(7, 2, 2, 0.5, "M(7,2,1/2) [banked gap]")
    run(7, -1, 2, 0.5, "M(7,-1,1/2) [Test A gap]")
    run(3, 2, 3, 0.5, "M(3,2,1/2) L=3 [r=1]")
    log("")
    log("## Q6-B -- is the q=3 cross-Perron collision WEIGHT(lam)-BLIND? M(3,2,lam), gen=2 fixed (LIFTS)")
    log(f"   {'':<22} {'':<10} {'':<11} {'slope B (=k-prefactor; 0=>unglued)':>10}")
    res = {}
    for lam in [0.25, 1 / 3, 0.4, 0.5, 0.6, 0.7]:
        B, v = run(3, 2, 2, lam, f"M(3,2,{lam:.3f})")
        res[lam] = (B, v)
    log("")
    ncoll = sum(1 for lam, (B, v) in res.items() if "COLLISION" in v)
    log(f"## VERDICT Q6-B: {ncoll}/{len(res)} lam values COLLIDE.")
    if ncoll == len(res):
        log("   => cross-Perron collision is WEIGHT-BLIND (lam-inessential). Combined with Q6-A (phase")
        log("      group <2>'s q-adic lift IS essential), the theorem is: d=2 boundary carried by the")
        log("      PHASE structure alone; weights quotient out. Phase 2b proof = pure <2>-lift combinatorics.")
    else:
        log("   => collision UNGLUES at some lam => WEIGHT is load-bearing too; Phase 2b must keep the")
        log("      weight base (2^-v) in the statement, not just the phase group.")
    with open("result_phase2a_q6b_detector_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
