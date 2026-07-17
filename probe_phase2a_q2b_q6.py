"""
PHASE 2a -- Q2b (EP Puiseux signature) + Q6 (Mr. Potato Head decoupling). Recon, no proofs.

Q2b: deform weights p_v -> p_v*(1+eps*f(v)), f(v)=v (KNOB = weight deformation). Measure top-pair
   splitting |l1-l2|(eps), eps=1e-2..1e-8, log-log slope. PRE-REG: q=3 slope 1/2 (sqrt-eps = genuine
   2nd-order EP); q=7 CONTROL slope 1 (honestly distinct). Physics-frame adjudication of Q2.

Q6: generalized M(q, gen, lam): phases <gen> mod q^L, weights ~ lam^v. GATE G6: M(q,2,1/2)==build_M
   (r_3=1 collision, r_7~0.38). ROLE HYPOTHESIS test:
   TEST A (swap phase gen, keep weight): collision/Jordan FOLLOWS ord(gen)=2. M(7,-1,1/2), M(5,-1,1/2)
     critical (r=1) despite q!=3; ord(gen)>=3 gapped. {+-1} case-table is weight-blind + 2-blind.
   TEST B (swap weight lam, keep gen=2 at q=3): cross-Perron degeneracy is lam-BLIND. M(3,2,lam),
     lam=1/3,2/5,3/5 all keep r=1 (cross glued to Perron); within/Perron ratio moves (Sum p_v^2).
   ADJUD: A&B pass -> theorem is d=2 ALONE, proof pure phase-combinatorics (2/weights quotient out).

Cheap by construction (banked matrix sizes). r via ESPRIT on mass a_k=1^T M^k v0 (repeated top = collision).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from probe_25_transfer_operator_Aprime import build_M as build_M_banked
from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    try:
        print(m, flush=True)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def subgroup(gen, qL):
    sub = [1]; x = gen % qL
    while x != 1:
        sub.append(x); x = (x * gen) % qL
        if len(sub) > qL:
            break
    return sub


def build_M_gen(q, L, gen, raw_w):
    """Generalized pair operator: phases <gen> mod q^L, per-coordinate raw weights raw_w (normalized).
       gen=2, raw_w=[lam^d] reproduces build_M for lam=1/2."""
    qL = q ** L
    sub = subgroup(gen, qL)
    ordn = len(sub)
    geninv = pow(gen, -1, qL)
    w = np.array(raw_w[:ordn], dtype=float); w = w / w.sum()
    mult = [(pow(geninv, delta, qL), w[delta - 1]) for delta in range(1, ordn + 1)]
    states = [(a, b, g) for a in sub for b in sub for g in range(qL)]
    idx = {s: i for i, s in enumerate(states)}; n = len(states)
    rows, cols, vals = [], [], []
    for (a, b, g) in states:
        i = idx[(a, b, g)]
        for (ga, wa) in mult:
            ap = (a * ga) % qL
            for (gb, wb) in mult:
                bp = (b * gb) % qL
                T = (ap - bp) % qL
                if (g + T) % q == 0:
                    gp = ((g + T) // q) % qL
                    rows.append(idx[(ap, bp, gp)]); cols.append(i); vals.append(wa * wb)
    M = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return M, idx, n


def mass_seq(M, idx, K=13):
    n = M.shape[0]; v = np.zeros(n); v[idx[(1, 1, 0)]] = 1.0
    a = []
    for _ in range(K):
        v = M.dot(v); a.append(v.sum())
    return a


def esprit(seq, p):
    s = np.array(seq, float); N = len(s); L = max(N // 2, p + 1)
    if N - L < p:
        return None
    H = np.array([s[i:i + L] for i in range(N - L + 1)])
    U, S, Vt = np.linalg.svd(H, full_matrices=False)
    Us = U[:, :min(p, len(S))]
    return np.linalg.eigvals(np.linalg.pinv(Us[:-1, :]) @ Us[1:, :])


def r_from_mass(M, idx):
    a = mass_seq(M, idx)
    lam = esprit(a, 3)
    if lam is None:
        return float('nan'), a
    mags = sorted([abs(x) for x in lam], reverse=True)
    return (mags[1] / mags[0] if len(mags) > 1 and mags[0] > 0 else float('nan')), a


def top_pair_gap(M):
    n = M.shape[0]
    if n <= 1000:
        w = np.linalg.eigvals(M.toarray())
    else:
        w = spla.eigs(M, k=4, which='LM', return_eigenvectors=False)
    w = sorted(w, key=lambda z: -abs(z))
    return abs(w[0] - w[1])


def main():
    log("# PHASE 2a -- Q2b (EP Puiseux) + Q6 (Mr. Potato Head). Recon, no proofs.")
    log("")

    # ============ Q2b -- EP Puiseux signature (weight-deform knob, f(v)=v) ============
    log("## Q2b -- EP splitting |l1-l2|(eps), weight deform p_v->p_v(1+eps*v). PRE-REG q=3 slope~1/2, q=7 slope~1")
    eps_list = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
    for q, L in [(3, 3), (7, 2)]:
        ordL = order_of_two(q ** L)
        base = [0.5 ** d for d in range(1, ordL + 1)]
        gaps = []
        for eps in eps_list:
            raw = [base[d - 1] * (1 + eps * d) for d in range(1, ordL + 1)]
            M, idx, n = build_M_gen(q, L, 2, raw)
            gaps.append(top_pair_gap(M))
        # log-log slope (regression)
        le = np.log(eps_list); lg = np.log(gaps)
        slope = np.polyfit(le, lg, 1)[0]
        log(f"   q={q} L={L}: gaps={[f'{g:.3e}' for g in gaps]}")
        log(f"        log-log slope = {slope:.3f}   ({'~1/2 EP' if abs(slope-0.5)<0.15 else ('~1 distinct' if abs(slope-1)<0.2 else 'other')})")
    log("")

    # ============ Q6 -- Mr. Potato Head ============
    log("## Q6 GATE G6 -- M(q,2,1/2) reproduces banked (r_3=1 collision, r_7~0.38)")
    for q in [3, 7]:
        L = 2
        ordL = order_of_two(q ** L)
        r, a = r_from_mass(*build_M_gen(q, L, 2, [0.5 ** d for d in range(1, ordL + 1)])[:2])
        log(f"   M({q},2,1/2) L={L}: r(ESPRIT) = {r:.4f}  ({'collision r~1' if r>0.9 else f'gapped r={r:.2f}'})")
    log("")

    log("## Q6 TEST A -- swap phase gen, keep weight lam=1/2. PRE-REG: collision FOLLOWS ord(gen)=2")
    #   (q, gen, L, note)
    A_cases = [(7, -1, 2, "ord2 at q=7"), (5, -1, 2, "ord2 at q=5"),
               (13, -1, 2, "ord2 at q=13"), (13, 3, 1, "ord3 at q=13 (control)")]
    for q, gen, L, note in A_cases:
        g = gen % (q ** L)
        ordg = len(subgroup(g, q ** L)) if q ** L < 500000 else -1
        # weights lam^d over the subgroup order
        ordsub = len(subgroup(g, q ** L))
        raw = [0.5 ** d for d in range(1, ordsub + 1)]
        M, idx, n = build_M_gen(q, L, gen, raw)
        r, a = r_from_mass(M, idx)
        ordq = len(subgroup(gen % q, q))
        log(f"   M({q},{gen},1/2) L={L} [{note}]: |<gen> mod q|={ordq}, dim={n}, "
            f"r={r:.4f}  ({'COLLISION r~1' if r>0.9 else f'GAPPED r={r:.3f}'})")
    log("")

    log("## Q6 TEST B -- swap weight lam, keep gen=2 at q=3. PRE-REG: cross-Perron degeneracy lam-BLIND")
    log(f"   {'lam':>6} {'dim':>6} {'r(cross/Perron)':>16} {'lam1(Perron)':>13} {'Sum p_v^2 (dilution)':>20}")
    for lam in [1 / 3, 2 / 5, 3 / 5, 1 / 2]:
        L = 2
        ordL = order_of_two(3 ** L)
        raw = [lam ** d for d in range(1, ordL + 1)]
        M, idx, n = build_M_gen(3, L, 2, raw)
        r, a = r_from_mass(M, idx)
        lam1 = a[-1] / a[-2] if len(a) > 1 and a[-2] != 0 else float('nan')
        # Sum p_v^2 (dilution factor) with these normalized weights
        w = np.array(raw) / np.sum(raw)
        sump2 = float(np.sum(w ** 2))
        log(f"   {lam:>6.3f} {n:>6} {r:>16.4f} {lam1:>13.5f} {sump2:>20.5f}  "
            f"{'CROSS GLUED (r~1)' if r>0.9 else f'unglued r={r:.2f}'}")
    log("")
    log("## ADJUD: A&B pass => d=2 ALONE is the theorem; proof = pure phase combinatorics (weights/2 quotient out).")
    with open("result_phase2a_q2b_q6_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
