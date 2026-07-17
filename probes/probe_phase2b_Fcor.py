"""
Gates on LEMMA D1-MAX (user's proof). Two NEW falsifiable predictions beyond
Request F (which only measured L in {1,2}):

 P-F4  (index law):  nilpotency index of the e=-1 block = 2L, so at L=3 -> 6,
        with rho_- still 0 (q=5,7). Also confirm the surviving T-alphabet is
        exactly {2, q^L-2}.
 P-F5  (q=3 corollary): the e=-1 block at q=3 is NOT nilpotent -- it carries
        the pure-A fixed point gamma=1 at weight s^2, so rho_-(M(3,-1,l)) = s^2
        EXACTLY, and maximality (s^2 < |lam2|) holds iff 2 l^2 < 1 <=> l < 1/sqrt2
        ~ 0.70711. Map the crossover: rho_- < |lam2| below, rho_- = s^2 > |lam2|
        above. lambda=1/2 (Syracuse) is safely below.

Dense exact eig; principal-submatrix extraction (exits dropped, no renorm).
"""
import numpy as np

from probe_phase2a_q2b_q6 import build_M_gen, subgroup
from probe_phase2b_F import e_minus_block


def nil_index(Msub):
    P = Msub.copy()
    for k in range(1, Msub.shape[0] + 2):
        P = P @ Msub
        if np.abs(P).max() < 1e-13:
            return k + 1
    return None


def t_alphabet(q, L, lam):
    """Reconstruct the surviving within-sector T values directly from the block build."""
    qL = q ** L
    gen = (-1) % qL
    ordn = len(subgroup(gen, qL))
    raw = [lam ** d for d in range(1, ordn + 1)]
    M, idx, n = build_M_gen(q, L, gen, raw)
    Ts = set()
    inv = {v: k for k, v in idx.items()}
    Mc = M.tocoo()
    for r, c, v in zip(Mc.row, Mc.col, Mc.data):
        (a, b, g) = inv[c]
        (ap, bp, gp) = inv[r]
        if b == (-a) % qL and bp == (-ap) % qL:       # within e=-1 sector
            Ts.add((ap - bp) % qL)
    return sorted(Ts)


def main():
    print("# GATES on LEMMA D1-MAX\n")

    print("## P-F4 -- nilpotency index at L=3 (predict 2L=6), rho_-=0, T-alphabet={2,q^L-2}")
    for q in [5, 7]:
        for L in [1, 2, 3]:
            Msub, S, idx, states_S = e_minus_block(q, L, 0.5)
            rho = float(np.max(np.abs(np.linalg.eigvals(Msub))))
            ni = nil_index(Msub)
            alpha = t_alphabet(q, L, 0.5)
            pred = {2, (q ** L - 2)}
            ok_alpha = set(alpha) <= pred
            print(f"   q={q} L={L}  dim={Msub.shape[0]:>4}  rho_-={rho:.2e}  nil_index={ni} (2L={2*L})  "
                  f"T-alphabet={alpha} pred{{2,{q**L-2}}} [{'OK' if ok_alpha else 'DEVIATION'}]")
    print()

    print("## P-F5 -- q=3 toy corollary: rho_- = s^2 ? and maximality crossover at l=1/sqrt2")
    print(f"   1/sqrt2 = {1/np.sqrt(2):.6f}")
    for lam in [0.30, 0.50, 0.60, 0.70, 0.7071, 0.72, 0.80]:
        s2 = (lam / (1 + lam)) ** 2
        u2 = (1 / (1 + lam)) ** 2
        lam2 = (1 - lam) / (1 + lam)
        for L in [2, 3]:
            Msub, S, idx, states_S = e_minus_block(3, L, lam)
            rho = float(np.max(np.abs(np.linalg.eigvals(Msub))))
            eq_s2 = abs(rho - s2) < 1e-9
            hold = rho < lam2 - 1e-12
            pred_hold = (2 * lam * lam < 1)
            tag = "OK" if (hold == pred_hold) else "DEVIATION"
            if L == 2:
                print(f"   l={lam:.4f}: rho_-={rho:.6f}  s^2={s2:.6f}[{'=' if eq_s2 else 'NE'}]  u^2={u2:.6f}  "
                      f"|lam2|={lam2:.6f}  maximality {'HOLDS' if hold else 'FAILS'} "
                      f"(pred {'holds' if pred_hold else 'fails'}) [{tag}]", end="")
            else:
                print(f"   | L=3 rho_-={rho:.6f}[{'=s^2' if eq_s2 else 'NE'}]")
    print()

    print("## VERDICT")
    print("   P-F4 pass  <=> L=3 index is 6, rho_-=0, T-alphabet subset {2,q^L-2}.")
    print("   P-F5 pass  <=> rho_-(q=3)=s^2 and the maximality flip lands exactly at l=1/sqrt2.")


if __name__ == "__main__":
    main()
