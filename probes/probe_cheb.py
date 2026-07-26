"""
PROBE CHEB -- gate Wilson's two pen results (2026-07-25).

(1) CHEBYSHEV/COVARIANCE FORM: with a_j=rho_r(t+j m'), b_j=rho_r(t+k+j m'), m'=3^{r-1}:
      sum_t Cov_j(a,b) = (1/3)p_r(k) - (1/9)p_{r-1}(k) = A_r(k)/3^{r+1}.
    So q_r(k)>=1/3 <=> sum_t Cov_j >= 0. At k=0 Cov=Var (Cauchy-Schwarz => m=0 free, strict unless
    fiber-constant). At k=1 genuine covariance, sign unforced. GATE: identity to machine precision; k=0 all-Var>=0.

(2) x4 = NEAREST-NEIGHBOUR DIGIT OPERATION: 4X = X + 3X => digit_k(4X) = [digit_k(X)+digit_{k-1}(X)+carry] mod 3.
    GATE: verify digitwise on random integers. And 4^k=(1+3)^k => lag k couples k+1 digits w/ binomial weights.
    Print the measured channel ladder |Re dhat(k)| = |A_16(k)|/S_16 vs k (coupling width k+1) -- the
    "more digits coupled => faster equilibration => smaller excess" prediction, read honestly (incl. 3|k folds).
"""
import os, sys, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def main():
    t0 = time.time()
    print("# PROBE CHEB -- gate the Chebyshev/covariance form + the digit-coupling claim\n")

    # ---------- (1) covariance identity ----------
    print("## (1) sum_t Cov_j(a,b) == A_r(k)/3^{r+1}   (and k=0 => all variances >= 0)")
    for r in (12, 14, 16):
        N = 3 ** r; mp = N // 3
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho = rho / rho.sum()
        rp = np.load(os.path.join(SCRATCH, f"rho_r{r-1}.npy")) if r - 1 >= 12 else None
        A = rho.reshape(3, mp)                       # A[j, t] = rho(t + j m')
        for k in (0, 1, 2):
            b = np.roll(rho, -k).reshape(3, mp)
            cov_t = (A * b).sum(0) / 3.0 - A.sum(0) * b.sum(0) / 9.0
            lhs = float(cov_t.sum())
            p_r = float(np.dot(rho, np.roll(rho, -k)))
            fold = A.sum(0)                          # rho_{r-1} via tower
            p_rm = float(np.dot(fold, np.roll(fold, -k)))
            rhs = (3 * p_r - p_rm) / 9.0             # = A_r(k)/3^{r+1} * 3^{r-1} ... check: (1/3)p_r-(1/9)p_rm
            rel = abs(lhs - rhs) / abs(rhs) if rhs != 0 else abs(lhs)
            extra = ""
            if k == 0:
                nneg = int((cov_t < -1e-20).sum())
                extra = f"  | k=0: #negative variances = {nneg} (must be 0)"
            print(f"   r={r} k={k}: sum_t Cov = {lhs:+.6e}  vs (1/3)p_r-(1/9)p_(r-1) = {rhs:+.6e}  rel {rel:.1e}{extra}")
        # sign summary at k=1: fraction of t with positive covariance (the average must be +, individuals need not)
        b1 = np.roll(rho, -1).reshape(3, mp)
        cov1 = (A * b1).sum(0) / 3.0 - A.sum(0) * b1.sum(0) / 9.0
        print(f"   r={r} k=1: frac(t) with Cov_t>0 = {(cov1 > 0).mean():.4f}  (avg + but pointwise mixed = genuine covariance)")
        del rho, A
    print()

    # ---------- (2) digit identity ----------
    print("## (2) digit_k(4X) = [digit_k(X) + digit_(k-1)(X) + carry] mod 3  (4X = X + 3X)")
    rng = random.Random(7)
    bad = 0
    for _ in range(2000):
        X = rng.randrange(1, 3 ** 30)
        Y = 4 * X
        carry = 0
        for k in range(32):
            dk = (X // 3 ** k) % 3
            dkm = (X // 3 ** (k - 1)) % 3 if k >= 1 else 0
            s = dk + dkm + carry
            if s % 3 != (Y // 3 ** k) % 3:
                bad += 1; break
            carry = s // 3
    print(f"   2000 random X < 3^30: digit-recursion failures = {bad}  [{'PASS' if bad == 0 else 'FAIL'}]")
    print("   4^k binomial check: 4^2=16=(1,2,1)_3? ", [(16 // 3 ** i) % 3 for i in range(3)],
          "  4^3=64 digits", [(64 // 3 ** i) % 3 for i in range(4)], "(binomial (1,3,3,1) after carries)")
    print()

    # ---------- channel ladder vs coupling width ----------
    print("## channel ladder at r=16: |Re dhat(k)| = |A_16(k)|/S_16 vs lag k (coupling width k+1)")
    r = 16; N = 3 ** r
    rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho = rho / rho.sum()
    rp = np.load(os.path.join(SCRATCH, f"rho_r{r-1}.npy")); rp = rp / rp.sum()
    S16 = 3 ** r * float(np.dot(rho, rho)) - 3 ** (r - 1) * float(np.dot(rp, rp))
    print(f"   S_16 = {S16:.6f}")
    print(f"   {'k':>3} {'A_16(k)':>12} {'|A|/S':>11}  3|k?")
    for k in range(1, 14):
        Ak = 3 ** r * float(np.dot(rho, np.roll(rho, -k))) - 3 ** (r - 1) * float(np.dot(rp, np.roll(rp, -k)))
        print(f"   {k:>3} {Ak:>+12.4e} {abs(Ak)/S16:>11.4e}  {'fold' if k % 3 == 0 else ''}")
    print("   [prediction: |excess| decreasing in coupling width for 3-nmid k; 3|k rows fold to lower level (different animal).]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
