"""
R3a2 — REFORMULATION: the 2^{-e} edge-weighted inverse-tree Plancherel mass.

R3a killed the uniform-count inverse measure: its D_n(k) is E_MAX-cutoff-dependent
by orders of magnitude, so "decays at 1/9" was a small-window + cutoff artifact.
The RIGHT object is the natural Syracuse harmonic/pushforward measure: each inverse
edge y -> g_-(y;e)=(2^e y-1)/3 carries weight 2^{-e} (the forward map's halving
probability). Cutoff-free because sum_e 2^{-e} < infinity.

Weighted measure at depth n, level k:
  W_n(r) = sum of path-weights 2^{-(sum e_i)} over depth-n vertices == r mod 3^k,
  mu_{n,k}(r) = W_n(r) / W_n   (r coprime to 3),  0 for r==0 mod 3.
  D_n(k) = 3^k sum_r mu(r)^2 - 3^{k-1} sum_s Q(s)^2   (same Plancherel identity).

Same exact precision-tower as R3a (residue transfer mod 3^{m} -> 3^{m-1}), but the
state carries WEIGHT not COUNT: child_weight += parent_weight * 2^{-e}.

PRE-REGISTERED HYPOTHESES (stated before running):
  H_WEIGHTED_CLEAN   : D_n(k) decays at a clean, E_MAX-STABLE geometric rate.
                       (If that rate is 1/9, the sibling 3^2 result is rescued
                        under the correct object; if another rate, report it.)
  H_STABLE_NONCLEAN  : E_MAX-stable but no clean geometric rate.
  H_STILL_CUTOFF     : still E_MAX-sensitive (would indicate a deeper problem;
                       unexpected given sum_e 2^{-e} converges).

VALIDATION: with edge weight set to 1 (not 2^{-e}), the tower must reproduce the
committed uniform-count D_n(k), n<=6 (result_inverse_tree_residue.csv). Gate on that.
"""
import csv, sys
from fractions import Fraction
from collections import defaultdict
sys.stdout.reconfigure(encoding="utf-8")

CSV_REF = r"C:\Collatz\result_inverse_tree_residue.csv"

def D_tower(n, k, E_MAX, weighted):
    """D_n(k). weighted=False -> exact unit-weight (Fraction, for validation);
    True -> 2^{-e} edge weight in FLOAT (rate/stability don't need exactness;
    exact Fractions blow up denominators ~2^{E_MAX*n})."""
    m = k + n
    if weighted:
        state = {1 % (3**m): 1.0}; zero = 0.0
        ew = [2.0**(-e) for e in range(E_MAX+1)]
    else:
        state = {1 % (3**m): Fraction(1)}; zero = Fraction(0)
    for d in range(n):
        mod_p = 3**m; mod_c = 3**(m-1)
        pow2 = [pow(2, e, mod_p) for e in range(E_MAX+1)]
        new = defaultdict(float) if weighted else defaultdict(Fraction)
        for rho, w in state.items():
            r3 = rho % 3
            if r3 == 0:
                continue
            e_start = 2 if r3 == 1 else 1
            for e in range(e_start, E_MAX+1, 2):
                if d == 0 and rho == 1 and e == 2:
                    continue
                child = ((pow2[e]*rho - 1) % mod_p // 3) % mod_c
                new[child] += w * ew[e] if weighted else w
        state = dict(new); m -= 1
    N = 3**k; Nm1 = 3**(k-1)
    total = sum(state.values())
    if total == 0:
        return zero
    sum_w2 = zero; grouped = defaultdict(float) if weighted else defaultdict(Fraction)
    for r, w in state.items():
        if r % 3 == 0:
            continue
        sum_w2 += w*w; grouped[r % Nm1] += w
    sum_Q2 = sum(w*w for w in grouped.values())
    return (N*sum_w2 - Nm1*sum_Q2) / (total*total)

if __name__ == "__main__":
    # ---- validation: unit weight reproduces committed uniform values ----
    ref = {}
    with open(CSV_REF) as f:
        for row in csv.DictReader(f):
            ref[(int(row["n"]), int(row["k"]))] = Fraction(int(row["D_n_k_num"]), int(row["D_n_k_den"]))
    ok = all(D_tower(n, k, 30, False) == v for (n, k), v in ref.items())
    print(f"VALIDATION (unit weight == committed uniform, {len(ref)} cells): {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("validation failed")

    # ---- weighted object: E_MAX stability + decay rate ----
    K_VALUES = [2, 3, 4, 5]
    NMAX = 14
    print("\n=== 2^{-e}-WEIGHTED D_n(k): E_MAX stability (k=4) ===")
    print(f"{'n':>2} " + "".join(f"{'E'+str(E):>14}" for E in (20, 30, 40)))
    for n in range(1, NMAX+1):
        vals = [float(D_tower(n, 4, E, True)) for E in (20, 30, 40)]
        print(f"{n:>2} " + "".join(f"{v:>14.6e}" for v in vals))

    print("\n=== weighted decay ratios D_{n+1}/D_n at E_MAX=40 ===")
    for k in K_VALUES:
        seq = [float(D_tower(n, k, 40, True)) for n in range(NMAX+1)]
        rats = [seq[i+1]/seq[i] if seq[i] > 0 else float('nan') for i in range(len(seq)-1)]
        print(f"k={k}: ratios(n>=3) = {[f'{r:.5f}' for r in rats[3:]]}")
    print("\nreference rates: 1/9=0.11111  1/3=0.33333  1/6=0.16667  1/12=0.08333")
