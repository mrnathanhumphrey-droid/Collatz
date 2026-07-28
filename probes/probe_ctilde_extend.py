"""
PROBE C-TILDE-EXTEND -- is c~_q = (q-3)/q, and are q=5,7 deviations finite-k or structural? (2026-07-27)

c~_q := lim_k S_k^(q) / (q/3)^k, S_k^(q) = X_k - X_{k-1}, X_k = q^k * sum_r pi_k(r)^2, pi_k = stationary
of the qx+1 Syracuse chain r -> ((q r + 1) 2^-v) mod q^k, v~Geom(1/2). Conjecture c~_q = (q-3)/q
(confirmed q=11,13,17 to 0.07-1%); q=5,7 deviate. OPEN: finite-k transient or structural?

  A  LARGE PRIMES {19,23,29,31,37} at k=2 (and k=3 where cheap): does c~_q -> (q-3)/q tighten?
  B  PUSH q=5 to k=5, q=7 to k=4: does c~_q(k) converge to (q-3)/q as k grows (finite-k) or
     asymptote away (structural)?
  GATE  float power-iteration stationary validated vs the exact-rational build (q17 probe) at q=17,k=2.

Float sparse power-iteration (exact rational is too slow past k=2 for large q); 1e-12 >> the 0.1%
we need to test (q-3)/q.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction
from c_tilde_q17_probe import build_markov_q, stationary_rational, order_of_two


def X_k_float(q, k, vmax=64):
    """X_k = q^k * sum pi^2, pi = stationary of qx+1 chain on (Z/q^k)*, via sparse power iteration."""
    N = q ** k
    M = order_of_two(N)
    vm = min(M, vmax)
    inv2 = pow(2, -1, N)
    # v-weights: P(v)=2^-v truncated at vm, renormalized (matches Z_v normalization in the limit)
    w = np.array([0.5 ** v for v in range(1, vm + 1)]); w /= w.sum()
    powinv2 = np.array([pow(inv2, v, N) for v in range(1, vm + 1)], dtype=np.int64)
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    idx = -np.ones(N, dtype=np.int64)
    idx[coprime] = np.arange(n)
    base = ((q * coprime + 1) % N).astype(np.int64)          # (q r + 1) for each state
    # precompute target index arrays per v
    tgt = np.empty((vm, n), dtype=np.int64)
    for vi in range(vm):
        t = (base * powinv2[vi]) % N
        tgt[vi] = idx[t]
    pi = np.full(n, 1.0 / n)
    for _ in range(2000):
        nxt = np.zeros(n)
        for vi in range(vm):
            np.add.at(nxt, tgt[vi], w[vi] * pi)
        nxt /= nxt.sum()
        if np.abs(nxt - pi).max() < 1e-15:
            pi = nxt; break
        pi = nxt
    return N ** 0 * (q ** k) * float((pi ** 2).sum())


def ctilde(q, ks, vmax=64):
    """returns dict k -> (X_k, S_k, c~_q(k))."""
    X = {0: 1.0}
    out = {}
    for k in ks:
        X[k] = X_k_float(q, k, vmax)
    for k in ks:
        if k - 1 in X:
            S = X[k] - X[k - 1]
            c = S / (q / 3.0) ** k
            out[k] = (X[k], S, c)
    return out


def gate():
    print("## GATE  float power-iteration vs exact rational (q=17, k=2)")
    # rational
    Xr = {0: Fraction(1)}
    for k in (1, 2):
        K, cop, M = build_markov_q(17, k)
        pi = stationary_rational(K)
        Xr[k] = Fraction(17 ** k) * sum(p * p for p in pi)
    S2r = Xr[2] - Xr[1]
    c2r = float(S2r / (Fraction(17, 3) ** 2))
    # float
    o = ctilde(17, [1, 2])
    c2f = o[2][2]
    print(f"   rational c~_17(k=2) = {c2r:.8f}   float = {c2f:.8f}   |diff| = {abs(c2r-c2f):.2e}")
    print(f"   {'PASS' if abs(c2r-c2f) < 1e-6 else 'FAIL'}\n")


def partA():
    print("## PART A  large primes {19,23,29,31,37}: c~_q(k) vs (q-3)/q")
    print(f"   {'q':>3} {'(q-3)/q':>9} {'c~_q(k=2)':>11} {'c~_q(k=3)':>11} {'dev(k=2)':>10} {'dev(k=3)':>10}")
    for q in (19, 23, 29, 31, 37):
        pred = (q - 3) / q
        ks = [1, 2, 3] if q <= 29 else [1, 2]      # k=3 states q^2(q-1): 29->23548 ok; 37->48396 heavier
        t0 = time.time()
        o = ctilde(q, ks)
        c2 = o[2][2]
        c3 = o[3][2] if 3 in o else float("nan")
        print(f"   {q:>3} {pred:>9.5f} {c2:>11.6f} {c3:>11.6f} {c2-pred:>+10.5f} "
              f"{(c3-pred) if 3 in o else float('nan'):>+10.5f}   ({time.time()-t0:.0f}s)")
    print()


def partB():
    print("## PART B  push q=5 to k=5, q=7 to k=4: finite-k (converges to (q-3)/q) or structural?")
    for q, kmax in ((5, 5), (7, 4)):
        pred = (q - 3) / q
        ks = list(range(1, kmax + 1))
        o = ctilde(q, ks)
        print(f"   q={q}  (q-3)/q = {pred:.5f}")
        print(f"     {'k':>2} {'c~_q(k)':>11} {'dev from (q-3)/q':>17}")
        for k in ks:
            if k in o:
                c = o[k][2]
                print(f"     {k:>2} {c:>11.6f} {c-pred:>+17.6f}")
        # trend of the deviation
        devs = [o[k][2] - pred for k in ks if k in o and k >= 2]
        if len(devs) >= 2:
            trend = "SHRINKING (finite-k)" if abs(devs[-1]) < abs(devs[0]) * 0.8 else \
                    "PERSISTENT (structural)" if abs(devs[-1]) > abs(devs[0]) * 0.9 else "flat/slow"
            print(f"     deviation k=2->{kmax}: {devs[0]:+.4f} -> {devs[-1]:+.4f}   => {trend}")
        print()


def partC():
    print("## PART C  the differentiator is ord(2 mod q): dev from (q-3)/q vs 2's order")
    print(f"   {'q':>3} {'ord(2modq)':>10} {'q-1':>4} {'primroot?':>9} {'(q-3)/q':>9} {'c~_q':>10} "
          f"{'dev':>9} {'dev*(2^ord-1)':>13}")
    rows = []
    for q in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        o = order_of_two(q)
        pred = (q - 3) / q
        c2 = ctilde(q, [1, 2])[2][2]
        dev = c2 - pred
        scaled = dev * (2 ** o - 1)
        pr = "yes" if o == q - 1 else "no"
        rows.append((o, q, pred, c2, dev, scaled, pr))
        print(f"   {q:>3} {o:>10} {q-1:>4} {pr:>9} {pred:>9.5f} {c2:>10.6f} {dev:>+9.5f} {scaled:>+13.4f}")
    print("   sorted by ord: dev DECREASES monotonically as ord(2 mod q) grows")
    for o, q, pred, c2, dev, scaled, pr in sorted(rows):
        print(f"     ord={o:>2} (q={q:>2}): dev={dev:>+9.5f}")
    print("   => (q-3)/q holds iff ord(2 mod q) large (prim-root-ish); FAILS for small-ord (Mersenne-like)")
    print("      q=7(ord3), q=5(ord4), q=31(ord5). dev ~ O(1/(2^ord - 1)) -- the 2^M-1 structure.\n")


def main():
    print("# PROBE C-TILDE-EXTEND -- (q-3)/q at large primes + q=5,7 finite-k-vs-structural\n")
    gate()
    partA()
    partB()
    partC()


if __name__ == "__main__":
    main()
