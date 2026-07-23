"""
probe_warp_W13.py -- the subordination "costume check" (Wilson's warp).

Claim: T_ours = E_v[ scale by 2^{-v} ] = (1/2 B)(I - 1/2 B)^{-1}, B = single fair-coin (v=1) step.
The certified R28 operator is:  nu_hat_r(xi) = e(xi/3^{r+1}) * E_v[ nu_hat_{r-1}(xi 2^{-v}) ].
Split into LINEAR CORE (the E_v resolvent, drop phase) and AFFINE DRESSING A_* (mult by e(xi/3^{r+1})).

W1 (core): build the unshifted core operator C on functions on Z/3^r, C f(xi)=Sum_v 2^{-v} f(xi 2^{-v}),
           read its subdominant |eigenvalue|.  Wilson pre-reg: 1/3 (resolvent of cos(pi/3)=1/2).
W3 (Jacobian): does the affine phase A_* = diag(e(xi/3^{r+1})) MOVE the spectrum (deformation, |lam2| real)
           or PRESERVE it (conjugation/similarity, |lam2|=core value, and 1/2 would be an artifact)?
           Compare spec(C) vs spec(T=A_* C).

CAVEAT stated plainly in output: this is the FIRST-MOMENT operator on nu_hat (fixed-level model). The
S_r rate |lam2|~1/2 is a SECOND-moment quantity. So this locates where 1/2 is NOT, as much as where it is.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def build_core(r, vmax_tol=1e-18):
    """C f(xi) = Sum_{v>=1} 2^{-v} f(xi * 2^{-v} mod 3^r).  Returns dense N x N (N=3^r)."""
    N = 3 ** r
    inv2 = pow(2, -1, N)
    C = np.zeros((N, N))
    for xi in range(N):
        u = inv2; v = 1; w = 0.5
        while w > vmax_tol:
            tgt = (xi * u) % N
            C[xi, tgt] += w
            u = (u * inv2) % N; v += 1; w *= 0.5
    return C


def phase_diag(r):
    """A_* = diag(e(2 pi i xi / 3^{r+1}))  (the '+1' affine shift, R28's e(xi/3^{r+1}))."""
    N = 3 ** r
    return np.array([cmath.exp(2j * math.pi * xi / (3 * N)) for xi in range(N)])


def topmods(M, k=6):
    ev = np.linalg.eigvals(M)
    ev = sorted(ev, key=lambda z: -abs(z))
    return ev[:k]


def fmt(ev):
    return "  ".join(f"{abs(z):.5f}" + (f"@{math.degrees(cmath.phase(z)):+.0f}" if abs(z.imag) > 1e-9 else "") for z in ev)


def main():
    print("# PROBE WARP W1/W3 -- subordination costume check (FIRST-MOMENT fixed-level model).\n")
    print("  1/3 = 0.33333   1/2 = 0.50000   1/sqrt7 = 0.37796   cos(pi/3)=0.5\n")
    for r in (3, 4, 5, 6):
        N = 3 ** r
        C = build_core(r)
        D = phase_diag(r)
        T = D[:, None] * C                      # A_* C  (left-multiply rows by phase)
        evC = topmods(C); evT = topmods(T)
        # subdominant = 2nd-largest modulus
        subC = abs(evC[1]); subT = abs(evT[1])
        print(f"r={r} (N={N}):")
        print(f"   CORE   (unshifted, W1)  lead={abs(evC[0]):.5f}  |lam2|={subC:.5f}   top: {fmt(evC)}")
        print(f"   SHIFTED(A_* C, W3)      lead={abs(evT[0]):.5f}  |lam2|={subT:.5f}   top: {fmt(evT)}")
        # W3 verdict: does the phase move |lam2|?
        moved = abs(subT - subC) / max(subC, 1e-12)
        print(f"   -> W3: |lam2| core={subC:.5f} vs shifted={subT:.5f}  (rel move {moved:.2%}) "
              f"[deformation if it moves; conjugation if not]")
        # restrict to UNITS only (3 not | xi): the multiplicative-character eigenbasis of the core
        units = [xi for xi in range(N) if xi % 3]
        Cu = C[np.ix_(units, units)]
        Du = D[units]
        Tu = Du[:, None] * Cu
        evCu = topmods(Cu); evTu = topmods(Tu)
        print(f"   (units-only) CORE |lam2|={abs(evCu[1]):.5f}   SHIFTED |lam2|={abs(evTu[1]):.5f}   "
              f"coreTop: {fmt(evCu[:4])}")
        print()
    print("  [analytic core eigenvalues on units = 1/(2*chi(2)-1) over characters chi; smallest-order-3 char"
          " chi(2)=e^{2pi i/3} -> 1/(2e^{2pi i/3}-1), |.|=1/sqrt7=0.378, NOT 1/3. Reported to expose the frame.]")


if __name__ == "__main__":
    main()
