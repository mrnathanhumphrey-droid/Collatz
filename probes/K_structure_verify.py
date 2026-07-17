"""
K_structure_verify.py — verify the structural claim:
  K_k has only one nonzero eigenvalue (= 1, the stationary), and K_k maps W_{k-1} -> 0.

The previous K_W_restricted_spectrum.py showed numerically K_W ~ 0 (~1e-17, machine epsilon).
This confirms via:
  1. Exact-rational K_k construction up to k=3 (V_MAX=20 cap, exact Q).
  2. Trace(K_k^m) for m=1, 2, 3 (sum of eigenvalues to m-th power).
     If true spec = {1, 0, ..., 0}, then trace(K_k^m) = 1 for all m >= 1.
  3. Direct exact-rank check (Gaussian elimination over Q).

Cross-check: row equality K_k(r, ·) = K_k(r + 3^{k-1}, ·) for r in (Z/3^k)*.
"""
import sys, os, json
from fractions import Fraction
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_K_exact(k, V_MAX=16):
    """Exact-rational K_k, normalized by Z = 1 - 2^{-V_MAX}."""
    N = 3 ** k
    inv2 = pow(2, -1, N)
    pow_inv2 = [pow(inv2, v, N) for v in range(1, V_MAX + 1)]
    coprime = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime)}
    n = len(coprime)
    Z = Fraction(2 ** V_MAX - 1, 2 ** V_MAX)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for r in coprime:
        i = state_idx[r]
        for v in range(1, V_MAX + 1):
            w = Fraction(1, 2 ** v) / Z
            target = ((3 * r + 1) * pow_inv2[v - 1]) % N
            j = state_idx.get(target)
            if j is not None:
                K[i][j] += w
    return K, coprime, state_idx


def matmul(A, B):
    n = len(A)
    p = len(B[0])
    m = len(B)
    C = [[Fraction(0) for _ in range(p)] for _ in range(n)]
    for i in range(n):
        for j in range(p):
            s = Fraction(0)
            for kk in range(m):
                s += A[i][kk] * B[kk][j]
            C[i][j] = s
    return C


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def rank_exact(A):
    """Exact rank over Q via Gaussian elimination."""
    n = len(A)
    m = len(A[0])
    M = [row[:] for row in A]
    r = 0
    col = 0
    while r < n and col < m:
        piv = -1
        for i in range(r, n):
            if M[i][col] != 0:
                piv = i
                break
        if piv == -1:
            col += 1
            continue
        if piv != r:
            M[r], M[piv] = M[piv], M[r]
        pv = M[r][col]
        for i in range(n):
            if i != r and M[i][col] != 0:
                f = M[i][col] / pv
                for j in range(col, m):
                    M[i][j] -= f * M[r][j]
        r += 1
        col += 1
    return r


def check_row_equality(K, coprime, k):
    """Check whether K_k(r, ·) = K_k(r + 3^{k-1}, ·) for all r."""
    state_idx = {r: i for i, r in enumerate(coprime)}
    N_km = 3 ** (k - 1)
    bad = 0
    for r in coprime:
        for shift in (N_km, 2 * N_km):
            rp = (r + shift) % (3 ** k)
            if rp % 3 == 0 or rp not in state_idx:
                continue
            i, j = state_idx[r], state_idx[rp]
            diff = max(abs(K[i][col] - K[j][col]) for col in range(len(K)))
            if diff != 0:
                bad += 1
                if bad <= 3:
                    print(f"  row r={r} != row r'={rp}, max diff = {diff}")
    if bad == 0:
        print(f"  k={k}: all 3-fiber-row equalities hold EXACTLY (K(r,.) = K(r+3^{k-1},.))")
    else:
        print(f"  k={k}: {bad} row-equality violations")
    return bad


def main():
    out = {}
    for k in (2, 3):
        print(f"\n=== k={k} ===")
        K, coprime, _ = build_K_exact(k, V_MAX=16)
        n = len(coprime)
        print(f"  N_k = {n}")
        # Row equality
        check_row_equality(K, coprime, k)
        # Exact rank
        r = rank_exact(K)
        print(f"  exact rank(K_{k}) over Q = {r}  (expected {n})")
        # Trace of K^m
        tr1 = trace(K)
        K2 = matmul(K, K)
        tr2 = trace(K2)
        K3 = matmul(K2, K)
        tr3 = trace(K3)
        print(f"  trace(K)   = {tr1} (= {float(tr1):.10f})")
        print(f"  trace(K^2) = {tr2} (= {float(tr2):.10f})")
        print(f"  trace(K^3) = {tr3} (= {float(tr3):.10f})")
        # rank(K^2)
        r2 = rank_exact(K2)
        print(f"  exact rank(K^2) = {r2}")
        if k == 3:
            r3 = rank_exact(K3)
            print(f"  exact rank(K^3) = {r3}")
        out[f"k={k}"] = {
            "N_k": n,
            "exact_rank_K": r,
            "trace_K": [tr1.numerator, tr1.denominator],
            "trace_K2": [tr2.numerator, tr2.denominator],
            "trace_K3": [tr3.numerator, tr3.denominator],
            "exact_rank_K2": r2,
        }
    with open(os.path.join(OUTDIR, "K_structure_verify.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "K_structure_verify.json"))
    print("""
INTERPRETATION
==============
If trace(K^m) = 1 for all m >= 1 (sum of eigenvalues^m), then the only
nonzero eigenvalue is 1 (with algebraic multiplicity 1). All other
eigenvalues are 0, possibly with Jordan blocks.

If rank(K_k) > rank(K_{k+1}^∞) = 1, then K_k is non-diagonalizable
(Jordan structure at 0). The rank tells us the size of the largest
Jordan block at 0.
""")


if __name__ == "__main__":
    main()
