"""
R3a — the 1/9 inverse-tree decay (PRE_REG_3A_INVERSE_TREE_NINTH).

Confirms the 3x+1 single-basin inverse-tree Plancherel mass D_n(k) decays at
EXACTLY 1/9, extending past the n<=6 explicit-tree wall via an exact precision-tower.

Method (extends result_inverse_tree_residue.py definitions, does not rebuild them):
  D_n(k) needs only vertex counts mod 3^k. Child residue mod 3^{m-1} is a function
  of parent residue mod 3^m and e (child=(2^e*y-1)/3). So track a count-vector over
  Z/3^m with m running k+n (depth 0, single vertex {1}) down to k (depth n),
  reducing one power of 3 per forward step. State size <= min(|V_d|, 3^m): small.
  Same inverse map, same E_MAX=30, same D-formula as the committed script.

Validated against result_inverse_tree_residue.csv (n=0..6, k=1..5) before trusting n>6.
"""
import csv, sys
from fractions import Fraction
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")
E_MAX = 30
CSV_REF = r"C:\Collatz\result_inverse_tree_residue.csv"

def D_tower(n, k):
    """Exact D_n(k) via the precision tower. Returns Fraction."""
    m = k + n                      # starting modulus exponent
    state = {1 % (3**m): 1}        # depth 0: single vertex value 1
    for d in range(n):             # step depth d -> d+1, mod 3^m -> mod 3^{m-1}
        mod_p = 3**m
        mod_c = 3**(m-1)
        pow2 = [pow(2, e, mod_p) for e in range(E_MAX+1)]   # 2^e mod 3^m
        new = defaultdict(int)
        for rho, cnt in state.items():
            r3 = rho % 3
            if r3 == 0:            # leaf: no predecessors
                continue
            e_start = 2 if r3 == 1 else 1
            for e in range(e_start, E_MAX+1, 2):
                if d == 0 and rho == 1 and e == 2:          # skip trivial self-edge
                    continue
                t = (pow2[e]*rho - 1) % mod_p               # == 0 mod 3 by validity
                child = (t // 3) % mod_c
                new[child] += cnt
        state = dict(new)
        m -= 1
    # state now = counts mod 3^k at depth n
    N = 3**k; Nm1 = 3**(k-1)
    total = sum(state.values())
    if total == 0:
        return Fraction(0)
    # coprime residues only for mu; sum_r mu^2 and Q(s)=sum over fibre mod 3^{k-1}
    sum_c2 = 0
    grouped = defaultdict(int)
    for r, c in state.items():
        if r % 3 == 0:
            continue
        sum_c2 += c*c
        grouped[r % Nm1] += c
    sum_g2 = sum(c*c for c in grouped.values())
    return Fraction(N*sum_c2 - Nm1*sum_g2, total*total)

def load_ref():
    ref = {}
    try:
        with open(CSV_REF) as f:
            for row in csv.DictReader(f):
                ref[(int(row["n"]), int(row["k"]))] = Fraction(int(row["D_n_k_num"]), int(row["D_n_k_den"]))
    except FileNotFoundError:
        pass
    return ref

if __name__ == "__main__":
    K_VALUES = [1, 2, 3, 4, 5]
    NMAX = 13
    print("R3a: inverse-tree 1/9 decay via exact precision tower")
    print("="*70)

    # ---- validate against committed n<=6 exact values ----
    ref = load_ref()
    print("Validation vs result_inverse_tree_residue.csv (n=0..6):")
    all_ok = True
    for (n, k), val in sorted(ref.items()):
        got = D_tower(n, k)
        ok = (got == val)
        all_ok = all_ok and ok
        if not ok:
            print(f"  MISMATCH n={n} k={k}: tower={got} ref={val}")
    print(f"  -> {'ALL MATCH' if all_ok else 'FAILURES ABOVE'} "
          f"({len(ref)} cells checked)")
    if not all_ok:
        sys.exit("Tower does not reproduce committed values; aborting.")

    # ---- extend ----
    print(f"\nExact D_n(k), n=0..{NMAX}:")
    D = {}
    for k in K_VALUES:
        for n in range(NMAX+1):
            D[(n, k)] = D_tower(n, k)

    rows = []
    print(f"\n{'k':>2} | ratios rho_n(k) = D_{{n+1}}/D_n  (target 1/9 = 0.11111)")
    for k in K_VALUES:
        print(f"k={k}:")
        for n in range(NMAX):
            a, b = D[(n, k)], D[(n+1, k)]
            if a == 0:
                r = None; rf = float('nan')
            else:
                r = b/a; rf = float(r)
            eq19 = (r == Fraction(1, 9)) if r is not None else False
            rows.append((n, k, a, b, rf, eq19))
            if n >= 1:  # skip the n=0 boot
                print(f"    n={n:>2}->{n+1:<2}: {rf:.8f}   {'= 1/9 EXACT' if eq19 else ''}")

    # ---- verdict inputs ----
    print("\nVERDICT INPUTS (late-n ratios, per k):")
    for k in K_VALUES:
        rs = [rf for (n, kk, a, b, rf, eq) in rows if kk == k and n >= NMAX-4]
        exacts = [eq for (n, kk, a, b, rf, eq) in rows if kk == k and n >= 2]
        print(f"  k={k}: last-4 ratios = {[f'{x:.6f}' for x in rs]}   "
              f"any exact 1/9 (n>=2)? {any(exacts)}   all exact 1/9 (n>=2)? {all(exacts)}")

    # ---- save ----
    with open(r"C:\Collatz\result_3a_dn_tables.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["n","k","D_num","D_den","D_float"])
        for k in K_VALUES:
            for n in range(NMAX+1):
                fr = D[(n,k)]
                w.writerow([n, k, fr.numerator, fr.denominator, float(fr)])
    print("\n[saved] result_3a_dn_tables.csv")
