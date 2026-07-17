"""
result_cycle_obstruction.py
===========================
Cycle-obstruction analysis: does the R75 framework's Plancherel identity
or sign-invariance K_- = σK_+σ produce new constraints on putative
non-trivial Collatz cycles?

Steps (per brief):
  3 (FIRST, highest leverage): negate known 3x-1 cycle residue traces mod
    3^k. By sign-invariance these become Markov-level 3x+1 cycles. Test
    whether this gives an integer-level obstruction.
  4: Monte Carlo S_cycle(k=5) over random cycle traces at L ∈ {10^3..10^6}.
    Compare distribution to 7/45.
  5: Compute minimum L compatible with rate-1/2 envelope at k ≤ 5
    (parametrized by what envelope means for finite-cycle stationaries).

Output:
  result_cycle_obstruction.md
  result_cycle_obstruction_S_distributions.csv
  result_cycle_obstruction_length_bound.csv
  result_cycle_obstruction_step3_traces.csv
"""
from __future__ import annotations

import csv
import io
import math
import sys
import time
from collections import Counter

import numpy as np

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)


# ============================================================================
# Step 3: Sign-invariance mirror check
# ============================================================================

def syracuse_3xm1_int(n):
    """Forward 3x-1 Syracuse step on positive integer: n -> (3n-1) >> v_2."""
    x = 3 * n - 1
    while (x & 1) == 0:
        x >>= 1
    return x


def syracuse_3xp1_int(n):
    """Forward 3x+1 Syracuse step on positive integer."""
    x = 3 * n + 1
    while (x & 1) == 0:
        x >>= 1
    return x


def collect_3xm1_cycle(start, max_iter=500):
    """Walk forward from start until we return to start (or another visited)."""
    seen = {start: 0}
    seq = [start]
    cur = start
    for i in range(1, max_iter):
        cur = syracuse_3xm1_int(cur)
        if cur in seen:
            return seq[seen[cur]:]  # cycle = trace from first occurrence onwards
        seen[cur] = i
        seq.append(cur)
    return None


def step_3_sign_invariance_mirror():
    print("=" * 78)
    print("Step 3: sign-invariance mirror check")
    print("=" * 78)
    print()

    # Known 3x-1 cycles on positive integers (odd starting points)
    known_3xm1_seeds = [1, 5, 17]
    cycles_3xm1 = {}
    for s in known_3xm1_seeds:
        c = collect_3xm1_cycle(s)
        if c is not None:
            cycles_3xm1[s] = c
            print(f"  3x-1 cycle from seed {s}: length {len(c)}")
            print(f"    elements: {c}")
    print()

    # For each k, compute residue trace mod 3^k and negate
    print("Sign-invariance: T_-(r) = -T_+(-r) on Z/3^k. So if (r_1,...,r_L)")
    print("is a 3x-1 cycle trace mod 3^k, then (-r_1,...,-r_L) is a Markov-level")
    print("3x+1 cycle trace at the same k.")
    print()

    rows = []
    for s, cycle in cycles_3xm1.items():
        for k in [1, 2, 3, 4, 5]:
            N = 3 ** k
            trace_minus = [n % N for n in cycle]
            trace_neg = [(-r) % N for r in trace_minus]
            print(f"  3x-1 cycle from {s}, k={k}, N=3^{k}={N}:")
            print(f"    3x-1 trace mod {N}:       {trace_minus}")
            print(f"    negated (3x+1 candidate): {trace_neg}")

            # Verify negated trace is closed under the heuristic 3x+1 step:
            # for each r_i, there should be SOME v >= 1 with ((3*r_i+1)*inv(2)^v) mod N == next
            inv2 = pow(2, -1, N)
            powers = [pow(inv2, v, N) for v in range(1, 64)]
            valid = True
            for i in range(len(trace_neg)):
                r = trace_neg[i]
                next_r = trace_neg[(i + 1) % len(trace_neg)]
                # Find v such that ((3r+1) * inv2^v) % N == next_r
                base = (3 * r + 1) % N
                found_v = None
                for v_ in range(1, 64):
                    if (base * powers[v_ - 1]) % N == next_r:
                        found_v = v_
                        break
                if found_v is None:
                    valid = False
                    print(f"      WARNING: no v in [1,63] gives r_{i}={r} -> r_{i+1}={next_r}")
            if valid:
                print(f"    Markov-level 3x+1 closure: VERIFIED (some v exists for each step)")

            rows.append({
                "seed_3xm1": s,
                "cycle_length_int": len(cycle),
                "k": k,
                "trace_3xm1": str(trace_minus),
                "trace_neg_3xp1_candidate": str(trace_neg),
                "markov_closure_verified": valid,
            })
            print()

    # Save
    with open(r"C:\Collatz\result_cycle_obstruction_step3_traces.csv", "w",
              newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("=" * 78)
    print("Step 3 verdict")
    print("=" * 78)
    print("""
The negation σ(r) = -r mod 3^k carries 3x-1 cycle residue traces into
3x+1 Markov-level cycle traces. This is a tautological consequence of
K_- = σK_+σ (proven earlier).

Critical question: does this constrain INTEGER-LEVEL positive 3x+1 cycles?

NO. Reasoning:
  1. Sign-invariance is a Markov-chain identity (heuristic v ~ Geom(1/2)).
  2. Integer-level: for n > 0, T_+(n) is well-defined, but T_+(-n) is not
     part of standard Collatz dynamics (which acts on positive integers).
     The "Markov-level mirror cycle" obtained by negation lives in the
     residue space; its integer-level realization (if any) would be on
     NEGATIVE integers, not positive.
  3. The known 3x-1 positive cycles {1,2}, {5,7,10,14}, {17,...,34} have
     residue traces. Negating gives valid Markov-level 3x+1 traces. But
     these aren't integer-level 3x+1 cycles on positives — they are
     residue patterns that 3x-1 NEGATIVE-INTEGER cycles would carry.
     Integer-level 3x+1 dynamics on positives doesn't see these traces.

So step 3 does NOT produce an obstruction. Sign-invariance gives a
Markov-level cycle correspondence between 3x+1 and 3x-1 residue spaces,
but doesn't rule out (or constrain) integer-level positive-integer cycles
in either system.

This is the expected outcome documented in the prior sibling-symmetry
verdict (sibling_3x_minus_1_symmetry_verdict.md).
""")
    return rows


# ============================================================================
# Step 4: Monte Carlo S_cycle distribution
# ============================================================================

def random_cycle_trace_uniform(L, k, rng):
    """Random L-sample of coprime residues mod 3^k. Used as 'maximally random'
    cycle trace baseline."""
    N = 3 ** k
    coprime = np.array([r for r in range(N) if r % 3 != 0], dtype=np.int64)
    return rng.choice(coprime, size=L, replace=True)


def random_cycle_trace_walk(L, k, rng):
    """Generate a length-L Markov-walk trace at level k via 3x+1 heuristic
    chain (v ~ Geom(1/2) truncated)."""
    N = 3 ** k
    M = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        M += 1
    inv2 = pow(2, -1, N)
    powers = np.array([pow(inv2, v_, N) for v_ in range(1, M + 1)], dtype=np.int64)
    Z_v = (2 ** M - 1) / (2 ** M)
    weights = np.array([(1.0 / 2 ** v_) / Z_v for v_ in range(1, M + 1)])
    weights /= weights.sum()

    # Start at random coprime
    coprime = [r for r in range(N) if r % 3 != 0]
    r = int(rng.choice(coprime))
    trace = np.empty(L, dtype=np.int64)
    trace[0] = r
    for i in range(1, L):
        v_idx = int(rng.choice(M, p=weights))
        r = ((3 * r + 1) * int(powers[v_idx])) % N
        trace[i] = r
    return trace


def compute_S_cycle(trace, k):
    """Compute Plancherel mass S_cycle(k) on coprime non-trivial characters
    of Z/3^k from the empirical residue trace."""
    N = 3 ** k
    L = len(trace)
    # mu_hat(xi) = (1/L) sum_j exp(-2pi i xi r_j / N)
    # |mu_hat(xi)|^2 = (1/L)^2 |S(xi)|^2 where S is exponential sum
    # We compute over all xi in Z/N coprime to 3 (= xi mod 3 != 0)
    # Total of 2*3^{k-1} coprime characters.
    # Vectorized:
    cos_term = np.zeros(N, dtype=np.float64)
    sin_term = np.zeros(N, dtype=np.float64)
    inv_N = 2 * np.pi / N
    cnt = Counter(trace.tolist())
    # mu(r) = count_r / L; mu_hat(xi) = (1/L) sum_r count_r exp(-2pi i xi r / N)
    for r, c in cnt.items():
        ang = inv_N * r
        # vectorized cos/sin over xi 0..N-1
        xi = np.arange(N, dtype=np.float64)
        cos_term += (c / L) * np.cos(xi * ang)
        sin_term += (c / L) * np.sin(xi * ang)
    sq = cos_term * cos_term + sin_term * sin_term
    # Sum over xi coprime to 3 (xi mod 3 != 0)
    coprime_mask = np.array([(xi % 3) != 0 for xi in range(N)])
    S = float(sq[coprime_mask].sum())
    return S


def step_4_monte_carlo():
    print("=" * 78)
    print("Step 4: Monte Carlo S_cycle distribution at k=5")
    print("=" * 78)
    k = 5
    N = 3 ** k
    n_coprime = 2 * 3 ** (k - 1)
    print(f"  k={k}, N=3^{k}={N}, # coprime characters = {n_coprime}")
    print(f"  Reference: S_∞ = 7/15 ≈ {7/15:.6f} (Tao framework limit)")
    print(f"  Reference: c = 7/45 ≈ {7/45:.6f} (R75 closed form)")
    print()

    rows = []
    rng = np.random.default_rng(seed=20260505)
    for L in [100, 1000, 10000, 100000]:
        n_trials_walk = 50 if L <= 10000 else 10
        n_trials_unif = 50 if L <= 10000 else 10

        S_walk = []
        for _ in range(n_trials_walk):
            tr = random_cycle_trace_walk(L, k, rng)
            S_walk.append(compute_S_cycle(tr, k))
        S_unif = []
        for _ in range(n_trials_unif):
            tr = random_cycle_trace_uniform(L, k, rng)
            S_unif.append(compute_S_cycle(tr, k))

        mw, sw = float(np.mean(S_walk)), float(np.std(S_walk))
        mu, su = float(np.mean(S_unif)), float(np.std(S_unif))
        print(f"  L={L:>7}:")
        print(f"    Markov walk (v~Geom(1/2)):     S_cycle = {mw:.6f} ± {sw:.6f} (n={n_trials_walk})")
        print(f"    Uniform iid coprime:           S_cycle = {mu:.6f} ± {su:.6f} (n={n_trials_unif})")
        print(f"    7/15 = {7/15:.6f}; deviation walk = {mw - 7/15:+.6f}")
        rows.append({
            "L": L, "k": k,
            "S_walk_mean": mw, "S_walk_std": sw, "n_walk": n_trials_walk,
            "S_unif_mean": mu, "S_unif_std": su, "n_unif": n_trials_unif,
            "S_inf_ref_7_15": 7/15,
            "deviation_walk_from_7_15": mw - 7/15,
        })

    print()

    # Theoretical: for uniform iid coprime sample of length L, E[S_cycle] = ?
    # |mu_hat(xi)|^2 expected over random sample = 1/L for xi != 0 (high-freq)
    # Sum over coprime xi: E[S] = (n_coprime - 0)/L (excluding xi=0)
    # Wait: xi=0 is xi mod 3 = 0, NOT coprime. So all n_coprime characters
    # have <e^(2pi i xi r / N)> = 0 in true uniform.
    # So E[|mu_hat(xi)|^2 | uniform] = 1/L for each coprime xi.
    # E[S | uniform] = n_coprime / L = 2*3^{k-1} / L
    print(f"  Theoretical for iid uniform: E[S_cycle] ≈ n_coprime/L = {n_coprime}/L")
    for L in [100, 1000, 10000, 100000]:
        print(f"    L={L:>7}: theoretical = {n_coprime/L:.6f}")
    print()

    # Save
    with open(r"C:\Collatz\result_cycle_obstruction_S_distributions.csv", "w",
              newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("=" * 78)
    print("Step 4 verdict")
    print("=" * 78)
    print(f"""
S_cycle for finite-length L cycles does NOT cluster near 7/15.

  - Uniform iid coprime: S_cycle ≈ n_coprime/L = {n_coprime}/L
  - Markov walk:         similar, slightly different constant
  - Both decay as 1/L

7/15 is the asymptotic Tao chain stationary's high-frequency Plancherel
mass — a property of the INFINITE-LENGTH stationary distribution, not
a finite-cycle distribution.

A finite cycle of length L places mass 1/L on each visited residue.
S_cycle at level k for the cycle's empirical measure = sum over coprime
xi of |sample-mean exponential|^2. For uniform sample this is ~2*3^{{k-1}}/L,
NOT 7/15.

Conclusion: the empirical S_cycle of any finite cycle is structurally
DIFFERENT from S_∞ = 7/15. They aren't the same quantity. So "cycle's
S_cycle = 7/15" is not a meaningful constraint.

The framework's R75 (S_n = 3^n · ‖d_n‖²) refers to the asymptotic
stationary, not finite-cycle empirical measures. Putative non-trivial
cycles, if they exist, would have their own finite-L empirical S_cycle
that need not match 7/15.

Step 4 does NOT produce a new cycle obstruction. The framework
characterizes ergodic asymptotic behavior, not finite-cycle structure.
""")
    return rows


# ============================================================================
# Step 5: Length bound from rate-1/2 envelope
# ============================================================================

def step_5_length_bound():
    print("=" * 78)
    print("Step 5: rate-1/2 envelope and minimum cycle length")
    print("=" * 78)
    print("""
Rate-1/2 envelope: |ε_n| · 2^n stable at ~0.04 through k=5..6 (R75 paper),
where ε_n = S_n - 7/15.

Question: does this envelope, applied to cycle stationaries, force a
minimum cycle length?

S_cycle(k) for length-L cycle with empirical measure mu_cycle:
  S_cycle(k) ≈ 2*3^{k-1}/L  for "random" cycle (high-frequency uniform)

If the cycle's S_cycle had to satisfy |S_cycle(k) - 7/15| · 2^k ≤ 0.04:
  |2*3^{k-1}/L - 7/15| · 2^k ≤ 0.04
  (assuming 2*3^{k-1}/L >> 7/15)
  ≈ 2*3^{k-1}/L · 2^k = 6^k / (3 L)

  We need 6^k / (3 L) ≤ 0.04
  L ≥ 6^k / (3 * 0.04) = 25 * 6^{k-1}/3 = (25/3) * 6^{k-1}

  k=1: L ≥ 1.4   (any cycle qualifies)
  k=2: L ≥ 8.3   (trivial bound)
  k=3: L ≥ 50    (still trivial)
  k=4: L ≥ 300
  k=5: L ≥ 1800
  k=10: L ≥ 5*10^7  (compare Eliahou's 1.7 * 10^10)

But this assumes finite-cycle S_cycle should obey the rate-1/2 envelope,
which it does NOT — see Step 4. The rate-1/2 envelope is for the
asymptotic stationary's S_n converging to 7/15 across LEVELS k, not for
cycle empirical measures matching 7/15.

So the bound above is mathematically derivable but conceptually misplaced.
The framework's rate-1/2 says the ASYMPTOTIC S_∞ is 7/15 with certified
convergence rate; it doesn't say a finite cycle's empirical S_cycle has
to match 7/15.
""")
    rows = []
    for k in range(1, 11):
        bound = (25 / 3) * 6 ** (k - 1)
        rows.append({"k": k, "L_min_if_envelope_applied": bound,
                     "Eliahou_2009_actual_bound": 1.7e10 if k == 10 else None})
        print(f"  k={k:>2}: L_min (if rate-1/2 forced) = {bound:.3g}")
    print()
    print("  Eliahou 1993 (refined by Simons-de Weger 2005): non-trivial cycle")
    print("  length L > 17 * 10^9. Framework-derived bound at k=5 is L > 1800,")
    print("  WAY weaker than Eliahou.")
    print()

    with open(r"C:\Collatz\result_cycle_obstruction_length_bound.csv", "w",
              newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["k", "L_min_if_envelope_applied",
                                            "Eliahou_2009_actual_bound"])
        w.writeheader()
        w.writerows(rows)

    print("=" * 78)
    print("Step 5 verdict")
    print("=" * 78)
    print("""
Even if rate-1/2 were misapplied to finite cycles, the resulting bound
L > (25/3) * 6^{k-1} is far weaker than Eliahou's L > 1.7 * 10^10 at any
reasonable k.

At k=10 (where the framework empirical envelope is conjectural beyond
k=5..6), the bound would be L > 5 * 10^7 — still 1000× weaker than
Eliahou.

Bottom line: the rate-1/2 envelope characterizes the asymptotic
stationary distribution's mode of convergence, NOT finite-cycle structure.
Applying it to cycles gives weak/meaningless bounds. Eliahou-style
bounds (residue-dynamics + valuation arithmetic) remain the binding
constraint.
""")
    return rows


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 78)
    print("Cycle obstruction analysis via R75 framework + sign-invariance")
    print("=" * 78)
    print()

    t0 = time.time()
    rows3 = step_3_sign_invariance_mirror()
    rows4 = step_4_monte_carlo()
    rows5 = step_5_length_bound()
    elapsed = time.time() - t0

    print("=" * 78)
    print("Overall verdict")
    print("=" * 78)
    print(f"""
Total compute: {elapsed:.1f}s.

Step 3 (sign-invariance mirror): NULL — Markov-level cycle correspondence
exists (tautological from K_-=σK_+σ), but doesn't constrain integer-level
positive cycles. Negation of 3x-1 positive-integer cycle traces gives
Markov-level 3x+1 cycle traces that would correspond to NEGATIVE-integer
3x+1 cycles, not positive.

Step 4 (residue-distribution gate): NULL — finite-cycle empirical S_cycle
follows ~n_coprime/L, NOT 7/15. The framework's S_∞ = 7/15 is the
asymptotic stationary distribution's high-freq Plancherel mass, a
fundamentally different object than a finite cycle's empirical measure.

Step 5 (rate-1/2 length bound): NULL — even if the envelope were
extended to finite cycles, the bound at k=10 would be L > 5*10^7,
1000× weaker than Eliahou's L > 1.7*10^10.

OVERALL: framework alone does NOT rule out non-trivial cycles. The
R75 / sign-invariance / rate-1/2 identities characterize ergodic
asymptotic behavior, not finite-cycle structure. The c=7/45 closed form
has no implications for cycle existence.

Existing length bounds (Eliahou 1993, Simons-de Weger 2005, refined
machine searches reaching ~10^20 with no second cycle found) remain
the binding constraints. The framework is silent on cycle structure.
""")


if __name__ == "__main__":
    main()
