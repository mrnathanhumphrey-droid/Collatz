"""
PROBE R4 -- THE EDGE DENSITY (thread 3 closing derivation). Dense L=2,3; NO near-EP extraction.
Family-sum: S_shell(m) = 3 * sum_j product_j (3 lambda_j)^{m-1},  product_j = A_j (lambda_j - 1/3),
A_j = <1|r_j><l_j|v0> the exact spectral amplitude (=psi_kin phi_j g_j, verified for the partner in R3-C).
Band modes j = the condensing complex pairs (k=1 doublet, k=2 seat, ...); theta_k = 2 pi k / 3^{L-1}.
R4-A  product ratio k=1 : k=2 = 2 : 1  (the 1/theta law).
R4-B  normalized density (3^{L-1}/2pi) theta_k product  L-invariant across L=2,3.
R4-C  arg(product_{k=1}) -> +-pi/2 (sine transform), both L.
R4-D  convention dump for the pen's prefactor assembly.
INSTRUMENT: dense eig at q=3 (band modes clean, well-separated from the 1/3 EP). No fit.
"""
import numpy as np
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

np.set_printoptions(linewidth=160, suppress=True)


def full_M(q, L, lam=0.5):
    qL = q ** L; sub = subgroup(2, qL); ordn = len(sub)
    raw = [lam ** d for d in range(1, ordn + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    return M.toarray(), idx, n


def spectral_setup(Md, idx, n):
    """Right (VR) and left (VL, from M^T) eigvecs. VR globally singular at the defective 1/3 EP,
       so amplitudes are taken PER MODE via matched left/right + bilinear l.r (band modes non-defective)."""
    ev, VR = np.linalg.eig(Md)
    evL, VL = np.linalg.eig(Md.T)
    v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
    one = np.ones(n)
    return ev, VR, evL, VL, one, v0


def grouped_amps(G, ev, VR, evL, VL, one, v0):
    """Robust amplitudes for a group G of (possibly near-degenerate) modes: invert the small g x g
       biorthogonality S = L^T R (bilinear). Returns A_j (j in G). Well-conditioned when G is separated."""
    R = VR[:, G]
    Ls = [VL[:, int(np.argmin(np.abs(evL - ev[j])))] for j in G]
    Lm = np.column_stack(Ls)
    S = Lm.T @ R                            # g x g bilinear (numpy @ does NOT conjugate)
    rw = one @ R                            # <1|r_j>, g-vector
    lv = Lm.T @ v0                          # <l_j|v0>, g-vector
    Anorm = np.linalg.solve(S, lv)          # normalized left projections
    return rw * Anorm, np.linalg.cond(S)


def seat(k, L):
    th = 2 * np.pi * k / (3 ** (L - 1))
    return (1.0 / 3.0) * np.cos(th / 2) ** 2 * np.exp(1j * th), th


def run(L):
    print(f"\n{'='*86}\n## L={L}")
    Md, idx, n = full_M(3, L)
    ev, VR, evL, VL, one, v0 = spectral_setup(Md, idx, n)
    up = [i for i in range(n) if ev[i].imag > 1e-9]      # upper-half complex representatives
    out = {}
    for k in (1, 2):
        s, th = seat(k, L)
        # k-pair = the 2 nearest upper-half modes to the seat (the ladder's two internal condensing modes)
        G = sorted(up, key=lambda i: abs(ev[i] - s))[:2]
        A, cond = grouped_amps(G, ev, VR, evL, VL, one, v0)
        prods = [A[t] * (ev[G[t]] - 1.0 / 3.0) for t in range(len(G))]
        pk = sum(prods)                                  # group-summed coupling-overlap product
        out[k] = dict(G=G, lams=[ev[j] for j in G], theta=th, product=pk, absprod=abs(pk),
                      arg=np.angle(pk), cond=cond, dens=(3 ** (L - 1) / (2 * np.pi)) * th * abs(pk))
        print(f"   k={k}: seat={s.real:+.4f}{s.imag:+.4f}j th_k={th:.4f} | pair lams "
              + ", ".join(f"{ev[j].real:+.5f}{ev[j].imag:+.5f}j" for j in G) + f"  (cond {cond:.1e})")
        print(f"        group product = sum A_j(lam-1/3) = {pk.real:+.6e}{pk.imag:+.6e}j  |.|={abs(pk):.6e}  arg={np.angle(pk):+.4f}")
    # R4-A ratio
    r = out[1]["absprod"] / out[2]["absprod"]
    print(f"   >> R4-A  |product(k=1)|/|product(k=2)| = {r:.4f}   (PRE-REG 2.00 = the 1/theta law; "
          f"{'MATCH' if abs(r-2) < 0.3 else f'dev {r-2:+.2f}'})")
    # R4-C arg
    print(f"   >> R4-C  arg(product_{{k=1}}) = {out[1]['arg']:+.4f}   (PRE-REG +-pi/2={np.pi/2:.4f}; "
          f"dev from pi/2 = {abs(abs(out[1]['arg'])-np.pi/2):.4f})")
    return out


def main():
    print("# PROBE R4 -- THE EDGE DENSITY. Dense L=2,3 band modes. Scaling SHAPES (no magnitude claim). No fit.")
    res = {}
    for L in (2, 3):
        res[L] = run(L)
    # R4-B L-invariance of normalized density (k=1)
    print(f"\n{'='*86}\n## R4-B  normalized edge density (3^{{L-1}}/2pi) theta_k |product|  (k=1; PRE-REG L-invariant)")
    for L in (2, 3):
        d = res[L][1]["dens"]
        print(f"   L={L}: (3^{L-1}/2pi)*theta_1*|product| = {d:.6e}"
              + ("   [L=2 super-critical caveat]" if L == 2 else ""))
    d2, d3 = res[2][1]["dens"], res[3][1]["dens"]
    print(f"   ratio L3/L2 = {d3/d2:.4f}  (PRE-REG ~1 within the L=2 regime caveat)")
    # R4-D convention dump
    print(f"\n{'='*86}\n## R4-D  CONVENTION DUMP (frozen verbatim for the pen's prefactor assembly)")
    print("   basis: M = build_M_gen(3,L,2,[2^-d]) dense; states (a,b,gamma), a,b in <2> mod 3^L, gamma in Z/3^L.")
    print("   right eigvecs r_j = columns of np.linalg.eig(M)[1]; left l_j = rows of VR^{-1} (=> l_j . r_j = delta_ij).")
    print("   readout <1| = all-ones vector (a_m = 1^T M^m v0 = P(pair agrees to depth m)).")
    print("   init v0 = delta(1,1,0) (index idx[(1,1,0)]); the independent pair both at phase 1, zero carry.")
    print("   spectral amplitude A_j = (1^T VR)_j * (VR^{-1} v0)_j ; product_j = A_j (lambda_j - 1/3).")
    print("   theta_k = 2 pi k / 3^{L-1} (seat/band angle). conjugate-pair: upper-half rep used; pair contrib = 2 Re(product_j).")
    print("   inner product: standard Euclidean (un-conjugated bilinear via VR^{-1}); NO l2-renormalization of eigvecs.")


if __name__ == "__main__":
    main()
