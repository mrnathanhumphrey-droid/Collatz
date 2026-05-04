"""
result_77_2_nisoli_application.py
=================================
R77.2 Stage 2: Nisoli Theorem 2.15 / Lemma 2.9 application with explicit
epsilon_N from Tao Prop 1.17.

NOTE on execution: this script was NOT run during the R77.2 session because
the harness denied python execution.  All numerical values below are
analytical predictions / explicit lower bounds where Tao's constants would
plug in.  See result_77_2_nisoli_certification.md for full discussion.

================================================================================

Plan
----
Apply Nisoli Lemma 2.9 to the order-3 companion T_3 (built in
result_77_2_T_N_construction.py, flavor B).  T_3 has spec {1/2, 1/4, 1/8}.
We want a contour gamma enclosing only the lambda = 1/2 eigenvalue.

Choose gamma = circle of radius 1/8 centered at 1/2:
- distance from 1/2 to 1/4 is 1/4 > 1/8  (1/4 stays outside)
- distance from 1/2 to 1/8 is 3/8 > 1/8  (1/8 stays outside)
- spectral isolation OK.

Compute M = sup_{z in gamma} ||R(z, T_3)||.

Then  ||P - P_K|| <= eps_K * M^2 * ell(gamma) / (2 (1 - eta))
  where eta = eps_K * M.

And  |lambda - lambda_K| <= [eps_K (1 + alpha) + 2 ||L_K|| alpha] / (1 - alpha)
  where alpha = ||P - P_K||.

================================================================================

Step 1: certified resolvent norm of T_3 along gamma.
=====================================================
Let T_3 = companion of (7/8, -7/32, 1/64).  Eigendecomposition:
  T_3 = V D V^{-1},  D = diag(1/2, 1/4, 1/8),
  V = Vandermonde [[1/4, 1/16, 1/64], [1/2, 1/4, 1/8], [1, 1, 1]]
    Wait, companion eigenvectors are (lambda^{r-1}, ..., lambda, 1), so
    V = [[1/4, 1/16, 1/64], [1/2, 1/4, 1/8], [1, 1, 1]]  (r=3)

Resolvent:  (z I - T_3)^{-1} = V (z I - D)^{-1} V^{-1}.
||(z I - T_3)^{-1}||_op (operator 2-norm)  =  ||V|| * ||V^{-1}|| * max_i 1/|z - lambda_i|.
                                          (upper bound; actual value <= this)

For z on circle of radius 1/8 around 1/2:
  |z - 1/2| = 1/8
  |z - 1/4| in [1/4 - 1/8, 1/4 + 1/8] = [1/8, 3/8]  ==> 1/|z-1/4| <= 8
  |z - 1/8| in [3/8 - 1/8, 3/8 + 1/8] = [1/4, 1/2]  ==> 1/|z-1/8| <= 4
  max_i 1/|z - lambda_i| <= 8.

cond_2(V) numerical (computed by hand):
  V[0] = [1/4, 1/16, 1/64], V[1] = [1/2, 1/4, 1/8], V[2] = [1, 1, 1]
  This is a Vandermonde at (1/2, 1/4, 1/8) scaled by powers.
  Numerically the condition number is moderate (a few dozen).
  We bound: cond_2(V) <= 100 conservatively (verifiable from explicit Vandermonde).

So  M_3 := sup_{z in gamma} ||R(z, T_3)|| <= 100 * 8 = 800.

This is the Nisoli quantity M for T_3 with gamma = circle radius 1/8 around 1/2.
ell(gamma) = 2*pi*1/8 = pi/4.

================================================================================

Step 2: epsilon_N from Tao Prop 1.17 -- where the trail goes cold
===================================================================

Tao's Proposition 1.17 (Tao 2022, Eq 1.25):
    |E e^{-2*pi*i*xi*Syrac(Z/3^n)/3^n}|  <=_A  n^{-A}      for any A > 0.

  Statement of constants (Tao 2022, page 13):
    "the implied constant in (1.25) is uniform in n and xi, though as
     indicated we permit it to depend on A".

Translation to our setting.
  The bilinear pair-form moment M_n(eta) = sum_xi mu_hat(xi) bar(mu_hat(xi*eta))
  has  |M_n(eta)|  <=  (2*3^{n-1}) * sup_xi |mu_hat(xi)|^2  <=  (2*3^{n-1}) * (C_A n^{-A})^2.

  This gives a bound  ||T - T_N||  in terms of the highest-frequency tail.  But:

   * Tao's C_A is NOT GIVEN EXPLICITLY in the paper.  It depends on:
        - the renewal-process white-point density (Section 7.2),
        - constants in Lemma 7.4 (triangle structure),
        - constants in Section 7.3 (probabilistic part).
     None are isolated as numerical bounds.

   * To produce an effective C_A (and hence epsilon_N = C_A * N^{-A}) one would
     need to redo Tao's combinatorial / probabilistic argument with effective
     constants tracked through:
        - the white-point lemma (7.2),
        - the triangle decomposition lemma (7.4),
        - the Pascal random-walk renewal estimate (7.5+).
     This is a non-trivial standalone project.

  STATUS: trail cold at the C_A extraction.  This blocks Stage 2 closure
  in the strict sense.  See result_77_2_nisoli_certification.md
  Outcome (delta).

================================================================================

Step 3: what we CAN do
======================
Even without C_A, we can do two things:

(a) Use empirical envelope from the project itself.  R75 measures
    |epsilon_n| * 2^n in [0.032, 0.041] for n = 2..6.  This gives a
    POSTERIORI rate-1/2 with constant <= 0.041, but it's not a Tao Prop 1.17
    application; it's an empirical fit through k = 6.

(b) Construct the conditional Theorem:
    "If C_A in Tao Prop 1.17 satisfies C_A <= F(A) for some explicit F, then
     ||T - T_N|| <= 2 * 3^N * F(A)^2 * N^{-2A} =: epsilon_N(A, F),
     and Nisoli Lemma 2.9 gives ||P - P_K|| <= eps_K(A,F) * 800^2 * pi/4 / (2(1-eta)),
     with eta < 1 once eps_K M < 1, i.e. once F(A)^2 N^{-2A} <= 1/(2*3^N * 800)."

    Conditional on Tao Prop 1.17 with explicit C_A, Nisoli machinery closes.

================================================================================
"""

import sys
import os
import math
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def companion(coeffs):
    r = len(coeffs)
    C = [[Fraction(0)] * r for _ in range(r)]
    C[0] = list(coeffs)
    for i in range(1, r):
        C[i][i - 1] = Fraction(1)
    return C


def vandermonde(roots):
    """Companion eigenvectors at given roots: V[i][j] = roots[j]^(r-1-i)."""
    r = len(roots)
    V = []
    for i in range(r):
        row = [Fraction(roots[j]) ** (r - 1 - i) for j in range(r)]
        V.append(row)
    return V


def matmul(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[Fraction(0)] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = Fraction(0)
            for t in range(k):
                s += A[i][t] * B[t][j]
            C[i][j] = s
    return C


def matinv_3x3(M):
    """Inverse of 3x3 over Q."""
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if det == 0:
        raise ValueError("singular")
    inv = [
        [(e * i - f * h) / det, -(b * i - c * h) / det, (b * f - c * e) / det],
        [-(d * i - f * g) / det, (a * i - c * g) / det, -(a * f - c * d) / det],
        [(d * h - e * g) / det, -(a * h - b * g) / det, (a * e - b * d) / det],
    ]
    return inv


def fro_norm(M):
    s = Fraction(0)
    for row in M:
        for x in row:
            s += x * x
    return math.sqrt(float(s))


def main():
    print("=" * 78)
    print("R77.2 Stage 2: Nisoli application (analytical, code not run).")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------ #
    # T_3 companion + eigendecomposition                                  #
    # ------------------------------------------------------------------ #
    coeffs = [Fraction(7, 8), Fraction(-7, 32), Fraction(1, 64)]
    T3 = companion(coeffs)
    print("Companion T_3 (over Q):")
    for row in T3:
        print(f"  {[str(x) for x in row]}")
    print(f"  spec(T_3) = {{1/2, 1/4, 1/8}}  (analytic)")
    print()

    roots = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8)]
    V = vandermonde(roots)
    print("V (companion eigenvectors columns):")
    for row in V:
        print(f"  {[str(x) for x in row]}")
    Vinv = matinv_3x3(V)
    print("V^{-1}:")
    for row in Vinv:
        print(f"  {[str(x) for x in row]}")
    print()

    # ------------------------------------------------------------------ #
    # Resolvent on circle gamma = circle radius 1/8 around 1/2            #
    # ------------------------------------------------------------------ #
    print("Contour gamma: circle of radius 1/8 around lambda_2 = 1/2.")
    print("  dist(gamma, 1/4) >= 1/4 - 1/8 = 1/8")
    print("  dist(gamma, 1/8) >= 3/8 - 1/8 = 1/4")
    print("  -> max_i sup_{z in gamma} 1/|z - lambda_i| = max(1/(1/8), 1/(1/8), 1/(1/4))")
    print("                                              = max(8, 8, 4) = 8")
    print(f"  ell(gamma) = 2 pi * (1/8) = pi/4 ~= {math.pi/4:.6f}")
    print()

    # Frobenius-norm bounds on V and V^{-1} (cheap upper bound on op-norm).
    V_fro = fro_norm(V)
    Vinv_fro = fro_norm(Vinv)
    cond_V_fro = V_fro * Vinv_fro
    print(f"||V||_F = {V_fro:.6f}")
    print(f"||V^-1||_F = {Vinv_fro:.6f}")
    print(f"||V||_F * ||V^-1||_F = {cond_V_fro:.6f}    (upper bound on cond_2 is this)")
    print()

    # M = max_z ||R(z, T_3)|| <= cond * 8
    M_bound = cond_V_fro * 8
    print(f"M_3 := sup_{{z in gamma}} ||R(z, T_3)|| <= {M_bound:.4f}")
    print(f"      (so for the Nisoli condition, eta = eps_N * M_3 < 1 requires "
          f"eps_N < 1/{M_bound:.0f}.)")
    print()

    # ------------------------------------------------------------------ #
    # epsilon_N from Tao Prop 1.17 -- DOCUMENT THE GAP                    #
    # ------------------------------------------------------------------ #
    print("Tao Prop 1.17 effective constant C_A:")
    print("  STATUS: NOT EXPLICITLY GIVEN in Tao 2022.")
    print("  Tao states the constant is 'uniform in n and xi, depending on A',")
    print("  i.e. the bound is qualitative.  Extracting an explicit C_A requires")
    print("  redoing Sections 7.2 and 7.3 of Tao 2022 with effective bookkeeping.")
    print("  -> Stage 2 closure: outcome (delta), see writeup.")
    print()

    # ------------------------------------------------------------------ #
    # CONDITIONAL bound: if C_A is known                                  #
    # ------------------------------------------------------------------ #
    print("Conditional bound (if C_A in Prop 1.17 is bounded by F(A)):")
    print("  ||T - T_N||_op <= 2 * 3^{N-1} * F(A)^2 * N^{-2A}")
    print("                 =: eps_N(A, F)")
    print(f"  Nisoli ||P - P_N|| <= eps_N * M^2 * ell(gamma) / (2 (1 - eta))")
    print(f"                     = eps_N * {M_bound**2:.0f} * (pi/4) / (2 (1 - eta))")
    print()

    # Save a single-row CSV summary
    out_csv = os.path.join(OUTDIR, "result_77_2_nisoli_bounds.csv")
    with open(out_csv, 'w') as f:
        f.write("quantity,value,note\n")
        f.write(f"contour_gamma_radius,1/8,circle around lambda_2=1/2\n")
        f.write(f"contour_gamma_length,{math.pi/4:.10f},pi/4\n")
        f.write(f"max_inv_dist_to_spec,8,1/(1/8) for nearest neighbor\n")
        f.write(f"V_fro,{V_fro:.6f},Frobenius norm\n")
        f.write(f"Vinv_fro,{Vinv_fro:.6f},Frobenius norm\n")
        f.write(f"cond_V_fro,{cond_V_fro:.6f},upper bound on cond_2\n")
        f.write(f"M_3_upper,{M_bound:.4f},sup_z ||R(z, T_3)|| upper bound\n")
        f.write("eps_N_constant,UNKNOWN,Tao Prop 1.17 C_A not effective in Tao 2022\n")
    print(f"[save] {out_csv}")


if __name__ == "__main__":
    main()
