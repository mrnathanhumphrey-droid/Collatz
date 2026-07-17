"""R3a E_MAX sensitivity: is the wild n>6 behavior intrinsic or an E_MAX=30 artifact?"""
import sys
from fractions import Fraction
from collections import defaultdict
sys.stdout.reconfigure(encoding="utf-8")

def D_tower(n, k, E_MAX):
    m = k + n
    state = {1 % (3**m): 1}
    for d in range(n):
        mod_p = 3**m; mod_c = 3**(m-1)
        pow2 = [pow(2, e, mod_p) for e in range(E_MAX+1)]
        new = defaultdict(int)
        for rho, cnt in state.items():
            r3 = rho % 3
            if r3 == 0: continue
            e_start = 2 if r3 == 1 else 1
            for e in range(e_start, E_MAX+1, 2):
                if d == 0 and rho == 1 and e == 2: continue
                t = (pow2[e]*rho - 1) % mod_p
                new[(t // 3) % mod_c] += cnt
        state = dict(new); m -= 1
    N = 3**k; Nm1 = 3**(k-1); total = sum(state.values())
    if total == 0: return Fraction(0)
    sum_c2 = 0; grouped = defaultdict(int)
    for r, c in state.items():
        if r % 3 == 0: continue
        sum_c2 += c*c; grouped[r % Nm1] += c
    return Fraction(N*sum_c2 - Nm1*sum(c*c for c in grouped.values()), total*total)

k = 4
print(f"D_n(k={k}) and ratio D_n/D_(n-1) at varying E_MAX  (target 1/9=0.11111)")
for E in (30, 40, 50, 60):
    print(f"\nE_MAX={E}:")
    prev = None
    for n in range(4, 13):
        D = float(D_tower(n, k, E))
        r = D/prev if prev else float('nan')
        print(f"  n={n:>2}: D={D:.4e}   ratio={r:.6f}")
        prev = D
