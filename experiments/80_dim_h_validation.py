"""
Bowen-Sullivan-Pollicott dimension validation for Result 23's M_closed.

Steps:
  1. Rebuild M_closed from natural-density rules at k in {5, 6, 8, 10}.
  2. Verify spectral properties: lambda_max real positive, spectral gap,
     uniqueness of leading eigvec.
  3. Compute dim_H candidate formulas:
       (a) log(lambda_max) / log(2)         -- Furstenberg-Hutchinson natural
       (b) 2 log(lambda_max) / log(2)       -- our claimed formula (factor 2)
       (c) log(lambda_max) / log(4/3)       -- using K_h scale
       (d) log(lambda_max) / log(3)         -- using 3 scale
  4. Compare to Chang's rigorous dim_H(C) = 0.6942 (Theorem-level, P(s)=0).
  5. Construct Chang's-style parametrized T(s) on inverse-tree residues:
       M_closed(s)[r, child] = weight * 2^(-v_step*s)
     where v_step = 1 for doubling (always) and v_step = 0 for inverse-3 (no
     2-adic valuation cost going backward via doubling = "+1 bit", inverse-3 = "0 bit").
     Actually for backward Collatz: doubling INCREASES n by factor 2, so "cost"
     = log_2(2) = 1; inverse-3 step n -> (n-1)/3 keeps n same magnitude order
     but has 2-adic no-op (since 1/3 is unit in Z_2), so "cost" = 0.
     Solve P(s) = log_2 lambda(M_closed(s)) = 0 for s.
  6. Sensitivity: vary k, check stability of all four estimates.

Output: dim_h_calculations.csv + dim_h_validation_log.txt
"""
import sys
import io
import numpy as np
from pathlib import Path

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path("C:/Collatz/experiments_output")


def build_M_closed(k):
    """Result 23's M_closed at modulus 2^k. Rows = residues, cols = residues.
    M[r, 2r mod M]  += 1                  doubling, ALWAYS
    M[r, child_r]   += 1/3                inverse-3, only EVEN r
    """
    Mod = 1 << k
    inv3 = pow(3, -1, Mod)
    M_mat = np.zeros((Mod, Mod), dtype=np.float64)
    for r in range(Mod):
        r_left = (2 * r) % Mod
        M_mat[r, r_left] += 1.0
        if r % 2 == 0:
            for j in range(3):
                lift = r + j * Mod
                if lift % 6 == 4:
                    r_right = ((lift - 1) * inv3) % Mod
                    M_mat[r, r_right] += 1.0 / 3.0
                    break
    return M_mat


def build_M_closed_parametrized(k, s, c_double=1.0, c_inv3=0.0):
    """
    Bowen-style parametrized: M(s)[r, child] = (multiplicity) * 2^(-c_step * s)
    c_double = "cost" of the doubling backward step
    c_inv3   = "cost" of the inverse-3 backward step
    Solve lambda_max(M(s)) = 1 for s = dim_H.
    """
    Mod = 1 << k
    inv3 = pow(3, -1, Mod)
    M_mat = np.zeros((Mod, Mod), dtype=np.float64)
    w_double = 2.0 ** (-c_double * s)
    w_inv3 = 2.0 ** (-c_inv3 * s)
    for r in range(Mod):
        r_left = (2 * r) % Mod
        M_mat[r, r_left] += w_double
        if r % 2 == 0:
            for j in range(3):
                lift = r + j * Mod
                if lift % 6 == 4:
                    r_right = ((lift - 1) * inv3) % Mod
                    M_mat[r, r_right] += (1.0 / 3.0) * w_inv3
                    break
    return M_mat


def leading_eig(M_mat):
    eigvals, eigvecs = np.linalg.eig(M_mat.T)
    order = np.argsort(-np.abs(eigvals))
    lam_max = eigvals[order[0]]
    v = eigvecs[:, order[0]]
    # Sub-dominant
    lam_sub = eigvals[order[1]] if len(order) > 1 else 0.0
    return lam_max, lam_sub, v, eigvals, order


def solve_pressure(k, c_double=1.0, c_inv3=0.0, lo=0.0, hi=2.0, tol=1e-9):
    """Bisect for s such that lambda_max(M_closed(s)) = 1."""
    def f(s):
        M = build_M_closed_parametrized(k, s, c_double, c_inv3)
        lam, _, _, _, _ = leading_eig(M)
        return float(lam.real) - 1.0
    f_lo = f(lo); f_hi = f(hi)
    if f_lo * f_hi > 0:
        return None  # no sign change in interval
    while hi - lo > tol:
        mid = (lo + hi) / 2
        f_mid = f(mid)
        if f_lo * f_mid <= 0:
            hi = mid; f_hi = f_mid
        else:
            lo = mid; f_lo = f_mid
    return (lo + hi) / 2


def main():
    print(f"=== Bowen-Sullivan-Pollicott Validation of dim_H ~ 0.6755 claim ===\n", flush=True)
    print(f"M_closed structure (Result 23):", flush=True)
    print(f"  M[r, 2r mod 2^k]   += 1       (doubling, always; n -> 2n in inverse tree)", flush=True)
    print(f"  M[r, child_r]      += 1/3     (inverse-3, only EVEN r)", flush=True)
    print(f"  Leading eigenvalue lambda_max governs per-layer GROWTH RATE of inverse tree.\n", flush=True)
    print(f"Chang 2603.11066v6 rigorous dim_H(C) = 0.6942 (verified 10 digits, line 7737).", flush=True)
    print(f"Chang's frequently-cited heuristic: ~0.68 (rounded).\n", flush=True)

    rows = []

    # ============ Step 1+2+3: Spectral verification + dim candidates ============
    print(f"=== Spectral verification + dimension candidates ===\n", flush=True)
    print(f"  {'k':>3}  {'lambda_max':>12}  {'lambda_sub':>12}  {'gap':>8}  {'real_pos?':>10}", flush=True)

    lam_max_values = {}
    for k in [5, 6, 8, 10]:
        Mod = 1 << k
        M = build_M_closed(k)
        lam, sub, v, all_eigs, order = leading_eig(M)
        lam_real = float(lam.real)
        lam_imag = float(lam.imag)
        sub_abs = float(np.abs(sub))
        gap = sub_abs / lam_real
        is_real_pos = abs(lam_imag) < 1e-10 and lam_real > 0
        print(f"  {k:>3}  {lam_real:>12.6f}  {sub_abs:>12.6f}  {gap:>8.4f}  {str(is_real_pos):>10}", flush=True)
        lam_max_values[k] = lam_real

    # Use k=10 reference (largest, most accurate)
    lam_max = lam_max_values[10]
    log_lam = np.log(lam_max)
    print(f"\n  Reference lambda_max (k=10): {lam_max:.9f}", flush=True)
    print(f"  log(lambda_max):              {log_lam:.9f}", flush=True)

    # All candidate formulas
    print(f"\n=== Candidate dimension formulas ===\n", flush=True)
    candidates = [
        ("log(lambda)/log(2)        [Furstenberg-Hutchinson natural]", log_lam / np.log(2)),
        ("2*log(lambda)/log(2)      [Our claimed formula, Result 57]", 2 * log_lam / np.log(2)),
        ("log(lambda)/log(4/3)      [K_h scale = log(2)/log(4/3)]",     log_lam / np.log(4/3)),
        ("log(lambda)/log(3)        [3x+1 scale]",                      log_lam / np.log(3)),
        ("log(lambda)/log(3/2)      [contraction rate]",                log_lam / np.log(3/2)),
        ("log(lambda)/log(2)+1      [related to density codim]",        log_lam / np.log(2) + 1),
        ("1 - log(lambda)/log(2)    [codimension]",                     1 - log_lam / np.log(2)),
    ]
    chang_dim = 0.6942
    chang_round = 0.68
    print(f"  {'Formula':<55}  {'Value':>10}  {'vs 0.6942':>10}  {'vs 0.68':>9}", flush=True)
    for name, val in candidates:
        d_chang = val - chang_dim
        d_round = val - chang_round
        print(f"  {name:<55}  {val:>10.6f}  {d_chang:>+10.6f}  {d_round:>+9.4f}", flush=True)
        rows.append({'formula': name, 'value': val, 'gap_vs_chang_0.6942': d_chang, 'gap_vs_round_0.68': d_round})

    # ============ Step 4: Bowen pressure equation, parametrized M_closed(s) ============
    print(f"\n=== Bowen pressure equation: M_closed(s)[r, c] = w * 2^(-cost * s) ===", flush=True)
    print(f"Solve lambda_max(M_closed(s)) = 1 for s = dim_H\n", flush=True)

    # Configuration 1: Standard self-similar IFS interp
    # doubling: cost = 1 (n -> 2n adds one bit, scale factor 2)
    # inverse-3: cost = 0 (no 2-adic scale change)
    print(f"  Config 1: doubling cost=1, inverse-3 cost=0 (natural 2-adic interp)", flush=True)
    s_solved = solve_pressure(8, c_double=1.0, c_inv3=0.0)
    if s_solved is not None:
        print(f"    s = {s_solved:.6f}  (vs Chang 0.6942)", flush=True)
        rows.append({'formula': 'Bowen P(s)=0, cost=(1, 0)', 'value': s_solved, 'gap_vs_chang_0.6942': s_solved - chang_dim, 'gap_vs_round_0.68': s_solved - chang_round})

    # Configuration 2: doubling cost=1, inverse-3 cost=log_2(3) (3-adic-like)
    print(f"  Config 2: doubling cost=1, inverse-3 cost=log_2(3)=1.585", flush=True)
    s_solved = solve_pressure(8, c_double=1.0, c_inv3=np.log2(3))
    if s_solved is not None:
        print(f"    s = {s_solved:.6f}  (vs Chang 0.6942)", flush=True)
        rows.append({'formula': 'Bowen P(s)=0, cost=(1, log_2(3))', 'value': s_solved, 'gap_vs_chang_0.6942': s_solved - chang_dim, 'gap_vs_round_0.68': s_solved - chang_round})

    # Configuration 3: both cost=1 (doubling and inverse-3 both contribute s)
    print(f"  Config 3: doubling cost=1, inverse-3 cost=1 (both contribute equally)", flush=True)
    s_solved = solve_pressure(8, c_double=1.0, c_inv3=1.0)
    if s_solved is not None:
        print(f"    s = {s_solved:.6f}  (vs Chang 0.6942)", flush=True)
        rows.append({'formula': 'Bowen P(s)=0, cost=(1, 1)', 'value': s_solved, 'gap_vs_chang_0.6942': s_solved - chang_dim, 'gap_vs_round_0.68': s_solved - chang_round})

    # Configuration 4: Match what we'd need to land at 0.6942
    # If formula is log(lambda)/log(scale), what scale gives 0.6942?
    needed_scale = np.exp(log_lam / chang_dim)
    print(f"\n  Reverse-engineer: log(lambda)/log(scale) = 0.6942 ⟹ scale = exp(log_lam/0.6942) = {needed_scale:.4f}", flush=True)
    print(f"    (compare scale=2 = 2.0000, scale=4/3 = 1.3333, scale=3 = 3.0000, scale=sqrt(2) = 1.4142)", flush=True)
    print(f"    Nearest: 4/3 = 1.3333 -- gap from {needed_scale:.4f} is {needed_scale - 4/3:+.4f}", flush=True)

    # ============ Step 5: Sensitivity analysis ============
    print(f"\n=== Sensitivity: dim estimates across k, candidate formulas ===\n", flush=True)
    print(f"  Best-natural-formula = log(lambda)/log(2):", flush=True)
    print(f"  {'k':>3}  {'lambda':>10}  {'log/log2':>10}  {'2*log/log2':>11}", flush=True)
    for k in [5, 6, 8, 10]:
        l = lam_max_values[k]
        ll = np.log(l)
        d_natural = ll / np.log(2)
        d_doubled = 2 * ll / np.log(2)
        print(f"  {k:>3}  {l:>10.6f}  {d_natural:>10.6f}  {d_doubled:>11.6f}", flush=True)

    # ============ Step 6: What does Chang's transfer matrix give for k=10? ============
    print(f"\n=== Chang's actual T(s) at safe classes {{1, 3, 7}} mod 8 ===", flush=True)
    print(f"T(s) = [[2^-2s, 2^-s,  0   ],", flush=True)
    print(f"        [2^-2s, 0,     2^-s],", flush=True)
    print(f"        [2^-2s, 0,     2^-s]]", flush=True)
    print(f"Indexed by safe forward Syracuse classes; weight 2^(-v*s) for valuation cost v.", flush=True)

    def chang_T(s):
        a = 2**(-2*s); b = 2**(-s)
        return np.array([[a, b, 0],
                         [a, 0, b],
                         [a, 0, b]])

    def chang_pressure(s):
        T = chang_T(s)
        eigvals = np.linalg.eigvals(T)
        return np.max(np.abs(eigvals)) - 1.0

    # bisect Chang's
    lo, hi = 0.5, 0.9
    while hi - lo > 1e-10:
        mid = (lo + hi) / 2
        if chang_pressure(mid) > 0: lo = mid
        else: hi = mid
    s_chang = (lo + hi) / 2
    print(f"\n  Chang's solved dim_H = {s_chang:.10f}  (paper claims 0.6942)", flush=True)
    rows.append({'formula': 'Chang T(s) Bowen pressure (verification)', 'value': s_chang, 'gap_vs_chang_0.6942': s_chang - chang_dim, 'gap_vs_round_0.68': s_chang - chang_round})

    # ============ Verdict ============
    print(f"\n=== VERDICT ===\n", flush=True)
    print(f"Our 2*log(lam)/log(2) = {2*log_lam/np.log(2):.6f}", flush=True)
    print(f"Chang's rigorous dim_H = {s_chang:.6f} (= 0.6942 in paper, verified)", flush=True)
    print(f"Heuristic 0.68 (rounded) = 0.6800", flush=True)
    print(f"", flush=True)
    print(f"Gap: |Our - Chang_rigorous| = {abs(2*log_lam/np.log(2) - s_chang):.4f}", flush=True)
    print(f"Gap: |Our - 0.68|           = {abs(2*log_lam/np.log(2) - 0.68):.4f}", flush=True)

    # save
    import polars as pl
    pl.DataFrame(rows).write_csv(OUT / "dim_h_calculations.csv")
    print(f"\n[save] dim_h_calculations.csv", flush=True)


if __name__ == "__main__":
    main()
