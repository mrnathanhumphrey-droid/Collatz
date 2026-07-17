"""
probe_doob_sigma_chain_2026_05_30.py

DOOB FRAMEWORK for c_inf. Derivation:
  D_{m+1} = (2^-a - 2^-b)(1 + qY_m) + q 2^-a D_m
where D_m = X_m - Y_m. Given collision D_m = q^m sigma_m, the next-level collision
condition is v_q(2^-a - 2^-b) >= m+1. For m >= 1 this dominantly requires a = b
(prob 1/3, matching empirical lambda).

In the dominant a=b regime: sigma_{m+1} = 2^-a * sigma_m, a multiplicative random walk
on (Z/q)*. So sigma_m mod q stays in <2>-coset of sigma_1. The limit dist on coset:
sigma_m mod q = sigma_1 * 2^{-(a_2+...+a_m)} mod q.  S_m mod ord_2(q) is sum of iid
Geom(2) draws, converging to uniform on Z/ord_2(q) by CLT/ergodicity.  Hence
sigma_inf mod q = uniform on sigma_1's <2>-coset.

For q=17: -1 = 2^4 in <2>, so both <2>-cosets are negation-closed. chi(Legendre)
is constant on each coset: +1 on <2> (= QR), -1 on the other.
=> c_inf = 2 * P(sigma_inf in <2>) - 1.
   P(sigma_inf in <2>) = P(sigma_1 in <2>) (coset is preserved under the multiplicative
   walk).

To test/use: compute P(sigma_1 in <2>) and see if it matches (1 + c_inf)/2 = 0.5765.

sigma_1 is set at depth 1. By the recursion analysis, sigma_1 depends on (a_1, b_1)
and the initial (X_0, Y_0) -- which for our deep-collision-conditioned analysis trace
back to the iid Syracuse pair (X, Y). Under the (a_1 = b_1) sub-case, sigma_1 mod q
= 2^-a_1 * (X-Y) mod q. So coset(sigma_1) = coset(X - Y) under this sub-case.

We need to be careful: conditioning on deep collision means we're sampling (a_i, b_i)
sequences and (X, Y) tail such that the chain survives. This re-weights things.

Strategy: simulate iid Syracuse pairs (X, Y), look for deep collisions, measure
coset(sigma_inf mod q) empirically. Check against c_inf prediction.
"""
from __future__ import annotations
import numpy as np
import sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8
TARGET_DEPTH = 6  # target deep-collision depth
N_TRIALS = 5_000_000

def legendre(x, q):
    x = int(x) % q
    return 0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)

def sample_syracuse_qadic(K, rng):
    """Sample one Syracuse-distributed q-adic integer to precision q^K.
    X = sum_{j=1..K} q^{j-1} * 2^{-S_j} mod q^K, S_j = a_1+...+a_j, a_i ~ Geom(2)."""
    N = q ** K
    inv2 = pow(2, -1, N)
    a = rng.geometric(0.5, size=K).astype(np.int64)  # Geom(p=1/2), values >= 1
    S = np.cumsum(a)
    X = 0
    # X mod q^K. 2^{-S_j} mod q^K: compute via pow.
    pow_inv2 = 1
    val = 1  # 2^{-0} = 1
    # iterate: 2^{-S_j} for j = 1..K. Need 2^{-a_1} for j=1, 2^{-(a_1+a_2)} for j=2, ...
    pow_q = 1  # q^{j-1}
    factor = pow(inv2, int(a[0]), N)  # 2^{-a_1}
    cur_inv2 = factor
    X = (X + pow_q * cur_inv2) % N
    for j in range(1, K):
        pow_q = (pow_q * q) % N
        cur_inv2 = (cur_inv2 * pow(inv2, int(a[j]), N)) % N
        X = (X + pow_q * cur_inv2) % N
    return int(X), S

def vq(d, q, max_v=20):
    """q-adic valuation of d (integer); returns 0..max_v."""
    if d == 0: return max_v
    v = 0
    while d % q == 0 and v < max_v:
        d //= q; v += 1
    return v

# build <2> subgroup mod q
H = set()
x = 1
while x not in H:
    H.add(x); x = (x * 2) % q
print(f"<2> mod {q} = {sorted(H)} (size {len(H)})")
# Check if -1 in <2>
neg1 = (q - 1) % q
print(f"-1 ({neg1}) in <2>? {neg1 in H}")

# Generate iid Syracuse pairs at precision q^K with K = TARGET_DEPTH + 2
K = TARGET_DEPTH + 2
rng = np.random.default_rng(42)
print(f"\nSampling {N_TRIALS} iid Syracuse pairs (X, Y) mod q^{K} = mod {q**K}...")

depth_counts = {m: 0 for m in range(TARGET_DEPTH + 2)}
deep_data = {m: {"unit_part": [], "chi_bar": []} for m in range(TARGET_DEPTH + 2)}

# batched sampling for efficiency
import time
t0 = time.time()
BATCH = 100000
n_batches = N_TRIALS // BATCH
for batch in range(n_batches):
    Xs = np.empty(BATCH, dtype=np.int64)
    Ys = np.empty(BATCH, dtype=np.int64)
    for i in range(BATCH):
        Xs[i], _ = sample_syracuse_qadic(K, rng)
        Ys[i], _ = sample_syracuse_qadic(K, rng)
    D = (Xs - Ys) % (q**K)
    # signed difference for v_q
    for i in range(BATCH):
        d = int(D[i])
        v = vq(d, q, K)
        if v >= TARGET_DEPTH + 2: v = TARGET_DEPTH + 1
        depth_counts[v] += 1
        if v >= 1 and v <= TARGET_DEPTH:
            unit_part = (d // (q ** v)) % q
            if unit_part != 0:
                deep_data[v]["unit_part"].append(unit_part)
                deep_data[v]["chi_bar"].append(legendre(unit_part, q))
    if batch % 5 == 4 or batch == n_batches - 1:
        elapsed = time.time() - t0
        print(f"  batch {batch+1}/{n_batches}: elapsed {elapsed:.1f}s")

print(f"\nDepth distribution P(v_q(X-Y) = m):")
total = sum(depth_counts.values())
for m in sorted(depth_counts):
    print(f"  m={m}: count={depth_counts[m]}, frac={depth_counts[m]/total:.6f}")

print(f"\nchi-bar moment at depth m (empirical, target c_inf):")
for m in range(0, TARGET_DEPTH + 1):
    if m in deep_data and len(deep_data[m]["chi_bar"]) > 100:
        mean = np.mean(deep_data[m]["chi_bar"])
        se = np.std(deep_data[m]["chi_bar"]) / np.sqrt(len(deep_data[m]["chi_bar"]))
        # P(unit_part in <2>)
        unit_parts = np.array(deep_data[m]["unit_part"])
        in_H = np.array([int(int(u) in H) for u in unit_parts])
        p_in_H = in_H.mean()
        # Predict chi_bar = 2*p_in_H - 1 (since <2> = QR mod 17, chi=+1 on <2>, -1 on other)
        pred_chi = 2 * p_in_H - 1
        print(f"  m={m}: n={len(deep_data[m]['chi_bar'])}, chi_bar={mean:+.6f} ± {se:.6f}, "
              f"P(unit in <2>)={p_in_H:.6f}, prediction 2P-1={pred_chi:+.6f}, match={abs(mean - pred_chi) < 1e-12}")

# Compare to known values:
print(f"\nKnown:")
print(f"  c(0) (exact)     = 19/127 = {19/127:.8f}")
print(f"  c_inf (Shanks)   = 0.152988994")
print(f"  prediction:  P(sigma_inf in <2>) such that 2P-1 = c_inf  ->  P = {(1 + 0.152988994)/2:.8f}")
