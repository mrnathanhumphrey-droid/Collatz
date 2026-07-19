"""
D-2 ladder judge item (committed before looking): the k=+-4-class pair near modulus (1/3)cos^2(1.40)~0.010,
phase ~2.8, in the dense L=3 spectrum. Also confirm THREE pair families at the top (k=+-1 doublet, k=+-2, k=+-4),
no more. Coprime classes mod 9 = {+-1,+-2,+-4}; top-mode phase of ladder k = 2*pi*k/3^{L-1} = 2*pi*k/9.
Report nearest modes, no fit.
"""
import numpy as np
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

def tower_dense(L, lam=0.5):
    q = 3; qL = q ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    M, idx, n = build_M_gen(q, L, 2, [lam ** d for d in range(1, D + 1)])
    gam = np.array([s[2] for s in sorted(idx, key=lambda s: idx[s])])
    tw = np.where(gam != 0)[0]
    return M[tw][:, tw].toarray(), D

def seat(k, L):
    th = k * 2 * np.pi / (3 ** (L - 1))
    return (1 / 3) * np.cos(th / 2) ** 2, th

def main():
    L = 3
    print(f"# D-2 k=+-4 hunt + three-family confirm, dense L={L} spectrum. Committed pre-reg.")
    Md, D = tower_dense(L)
    ev = np.linalg.eig(Md)[0]
    up = [z for z in ev if z.imag > 1e-9]        # upper-half conjugate representatives
    print(f"  {len(ev)} eigenvalues, {len(up)} upper-half complex pairs")
    # committed seats for the three coprime families
    for k in (1, 2, 4):
        mod, ph = seat(k, L)
        tgt = mod * np.exp(1j * ph)
        near = sorted(up, key=lambda z: abs(z - tgt))[:3]
        # also nearest purely by phase (the family fingerprint), among modulus < 0.5
        byph = sorted(up, key=lambda z: abs(((np.angle(z) - ph + np.pi) % (2*np.pi)) - np.pi))[:3]
        tag = "  <== COMMITTED k=+-4 HUNT" if k == 4 else ""
        print(f"\n  k=+-{k}: seat modulus={mod:.5f} phase={ph:.4f}{tag}")
        print(f"     nearest to seat point: " +
              ", ".join(f"{z.real:+.5f}{z.imag:+.5f}j (|.|={abs(z):.5f}, arg={np.angle(z):.4f})" for z in near))
        print(f"     nearest by PHASE:      " +
              ", ".join(f"{z.real:+.5f}{z.imag:+.5f}j (|.|={abs(z):.5f}, arg={np.angle(z):.4f})" for z in byph))
    # top-of-spectrum family census: the 3 largest-modulus complex pairs and their phases (should be k=1 doublet then k=2)
    print("\n  TOP complex pairs by modulus (family fingerprint via phase):")
    for z in sorted(up, key=lambda z: -abs(z))[:6]:
        print(f"     {z.real:+.6f}{z.imag:+.6f}j  |.|={abs(z):.5f}  arg={np.angle(z):.4f}")

if __name__ == "__main__":
    main()
