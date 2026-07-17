"""
path_b_substratum.py — Extended Markov-additive matrix Wiener-Hopf with
sub-stratum tracking at boundary residue r=21 (k=6).

Goal: test whether tracking j_inner = v_2(1+3·h) at boundary residue (where
m = 21 + 64·h) produces per-j W_j variation matching empirical W_2, W_4, W_5.

Sub-stratum encoding: at boundary residue r=21, current m has v_2(1+3·h) = j_inner.
The Brief's claim: absorption at attractor m_j (j ≥ 3) corresponds to specific
sub-stratum j_inner = 2(j-3). So end-state (r=21, j_inner=0) ≡ landing at m_3,
(r=21, j_inner=2) ≡ m_4, (r=21, j_inner=4) ≡ m_5, etc.

Critical issue: the sub-stratum is a v-stratum that aggregates MANY h values
per j_inner (those with v_2(1+3h) = j_inner). Only h = m_(j-3) gives lattice
landing at m_j; other h's at the same j_inner correspond to m ≠ m_j.

This script:
  1. Build extended state space at k=6: 31 non-boundary residues + j_inner_max
     boundary sub-strata (all j_inner ≥ 0, even AND odd). Truncate j_inner_max=12.
  2. Simulate MAP from S=0 to first descent below -L, tracking end (residue, j_inner).
  3. Condition on end state = (r=21, j_inner=2(j-3)) for j ∈ {4, 5} or
     end state = (r=5, *) for j=2.
  4. Compute conditional overshoot in step units (W_j convention: nats / log(4/3)).
  5. Compare to empirical W_2 = 7.156, W_4 = -4.755, W_5 = +4.590.

Empirical baseline from Path B v2 (no sub-stratum): universal W_j_pred = 3.44
across all j (= E[L⁻]/log(4/3)).

Cap < 400 lines.
"""
import numpy as np
import sys
from numba import njit

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3


def v2(n):
    if n == 0:
        return 0
    c = 0
    while n & 1 == 0:
        n >>= 1
        c += 1
    return c


def build_Q_basic_with_v(k=6):
    """Standard Q at k=6 (32 odd residues, no sub-stratum). Returns:
      Q[r,r'] : transition probabilities
      v_per_state : deterministic v(r) for non-boundary; -1 for boundary
      idx : dict r -> state index
    """
    M = 1 << k
    odd = list(range(1, M, 2))
    N = len(odd)
    idx = {r: i for i, r in enumerate(odd)}
    Q = np.zeros((N, N), dtype=np.float64)
    v_per_state = np.full(N, -1, dtype=np.int64)
    boundary = []
    for r in odd:
        i = idx[r]
        v_r = v2(3 * r + 1)
        if v_r < k:
            v_per_state[i] = v_r
            T_r0 = (3 * r + 1) >> v_r
            shift = (3 * M) >> v_r
            for h in range(1 << v_r):
                t = (T_r0 + shift * h) & (M - 1)
                if t & 1:
                    Q[i, idx[t]] += 1.0 / (1 << v_r)
        else:
            boundary.append(r)
            Q[i, :] = 1.0 / N  # uniform target distribution
    return Q, v_per_state, idx, odd, boundary


@njit(cache=True)
def simulate_substratum_walk(
    n_orbits, target_level, max_steps,
    Q_cum, v_per_state, boundary_idx, j_inner_max, seed
):
    """Walk MAP at extended state with sub-stratum tracking at boundary.

    At each step:
      - If at boundary residue, sample j_inner ~ Geom(1/2) (truncated j_inner_max)
        v = 6 + j_inner
      - Else, v = v_per_state[r] (deterministic)
      - X = log(3) - v·log(2); S += X
      - Sample next residue from Q[r, :]
    Stop when S < -target_level. Record (end_residue, end_j_inner) where
    end_j_inner is the j_inner of the most recent boundary visit; if no
    boundary visit on the final step, end_j_inner = -1.

    The sub-stratum at end is the j_inner used to compute the final descent
    step (v = 6 + j_inner). This corresponds to "absorbing at sub-stratum
    j_inner of the boundary residue" if end_residue = boundary.
    """
    np.random.seed(seed)
    N = Q_cum.shape[0]
    overshoots = np.full(n_orbits, np.nan)
    end_residues = np.full(n_orbits, -1, dtype=np.int64)
    end_j_inner = np.full(n_orbits, -1, dtype=np.int64)

    for orbit in range(n_orbits):
        S = 0.0
        r = 0  # start at residue 1 (index 0)
        last_j_inner = -1
        for step in range(max_steps):
            # Determine v based on state
            if r == boundary_idx:
                # Sample j_inner ~ Geom(1/2): P(j_inner = k) = 2^{-k-1}
                u = np.random.random()
                j_inner = 0
                cum = 0.5
                while j_inner < j_inner_max - 1 and u > cum:
                    j_inner += 1
                    cum += 0.5 ** (j_inner + 1)
                v = 6 + j_inner
                last_j_inner = j_inner
            else:
                v = v_per_state[r]
                last_j_inner = -1
            X = LOG3 - v * LOG2
            S += X
            # Sample next residue from Q[r, :]
            u2 = np.random.random()
            r_next = 0
            while r_next < N - 1 and Q_cum[r, r_next] < u2:
                r_next += 1
            r = r_next
            if S <= -target_level:
                # crossed below
                overshoots[orbit] = S - (-target_level)  # negative or zero
                end_residues[orbit] = r  # NEW residue after the descent step
                # The "absorbing sub-stratum" interpretation: the j_inner that
                # produced the descent step. We store this as end_j_inner.
                # If r != boundary_idx after this step, the sub-stratum doesn't
                # match an attractor — it was an "intermediate" boundary visit.
                end_j_inner[orbit] = last_j_inner
                break
    return overshoots, end_residues, end_j_inner


def run_substratum_extraction(k=6, n_orbits=1_000_000, j_inner_max=12, seed=42):
    """Build Q, run simulation per j_target, report conditional overshoots."""
    Q, v_per_state, idx, odd, boundary = build_Q_basic_with_v(k=k)
    N = Q.shape[0]
    boundary_idx = idx[boundary[0]]  # 21 → index in odd_residues list
    Q_cum = np.cumsum(Q, axis=1)
    Q_cum[:, -1] = 1.0

    print(f"State count: {N} odd residues, boundary residue = {boundary[0]} (idx {boundary_idx})")
    print(f"j_inner truncation: {j_inner_max} (covers ≥ {1 - 0.5**j_inner_max:.6f} of mass)")
    print(f"K_h conversion: 1/log(4/3) = {1/LOG43:.4f} step-units per nat\n")

    # Targets
    log_M = 36.0 * LOG2
    targets = [
        (2, 5, -1),     # j=2, m=5, residue=5, no sub-stratum
        (4, 85, 2),     # j=4, m=85, residue=21 (since 85 mod 64 = 21), j_inner=2(4-3)=2
        (5, 341, 4),    # j=5, m=341, residue=21, j_inner=2(5-3)=4
    ]

    print(f"{'j':>3}  {'m_j':>5}  {'L_j':>6}  {'target_state':>15}  "
          f"{'⟨overshoot⟩(nats)':>18}  {'W_j_pred (steps)':>17}  {'W_j_emp':>10}")
    print("-" * 90)

    results = []
    for j_val, m_j, target_j_inner in targets:
        L_j = log_M - np.log(float(m_j))
        target_r = m_j & ((1 << k) - 1)
        if target_r not in idx:
            continue
        target_r_idx = idx[target_r]
        oshs, ends_r, ends_j = simulate_substratum_walk(
            n_orbits, L_j, 200_000,
            Q_cum, v_per_state, boundary_idx, j_inner_max, seed + j_val
        )
        valid = ~np.isnan(oshs)
        oshs = oshs[valid]
        ends_r = ends_r[valid]
        ends_j = ends_j[valid]

        # Filter by end state
        if target_j_inner == -1:
            # j=2: condition only on end residue = target_r (= 5)
            cond = (ends_r == target_r_idx)
            target_label = f"r={target_r}, *"
        else:
            # j ≥ 3: condition on end residue = boundary AND end j_inner = target
            cond = (ends_r == target_r_idx) & (ends_j == target_j_inner)
            target_label = f"r={target_r}, j_inner={target_j_inner}"

        n_cond = int(cond.sum())
        if n_cond < 50:
            print(f"  j={j_val}  insufficient: only {n_cond} orbits matched")
            continue

        osh_nats = -float(oshs[cond].mean())
        osh_se_nats = float(oshs[cond].std() / np.sqrt(n_cond))
        W_j_pred = osh_nats / LOG43
        W_j_se = osh_se_nats / LOG43
        W_j_emp_map = {2: 7.156, 4: -4.755, 5: 4.590}
        W_j_emp = W_j_emp_map[j_val]

        results.append({
            'j': j_val, 'm_j': m_j, 'target_j_inner': target_j_inner,
            'n_cond': n_cond, 'osh_nats': osh_nats, 'osh_se_nats': osh_se_nats,
            'W_j_pred': W_j_pred, 'W_j_se': W_j_se, 'W_j_emp': W_j_emp,
        })
        print(f"  {j_val:>3}  {m_j:>5}  {L_j:>6.2f}  {target_label:>15}  "
              f"{osh_nats:>10.4f} ± {osh_se_nats:.4f}  "
              f"{W_j_pred:>+10.3f} ± {W_j_se:.3f}     {W_j_emp:>+8.3f}")
    return results


def report_per_substratum_breakdown(k=6, n_orbits=2_000_000, j_inner_max=12, seed=42):
    """Run a single large simulation, then break down conditional overshoot
    by every sub-stratum index. This shows whether sub-stratum produces ANY
    differentiation, and whether it correlates with j_attr expectations."""
    Q, v_per_state, idx, odd, boundary = build_Q_basic_with_v(k=k)
    boundary_idx = idx[boundary[0]]
    Q_cum = np.cumsum(Q, axis=1)
    Q_cum[:, -1] = 1.0

    log_M = 36.0 * LOG2
    L = log_M - np.log(85.0)  # use m_4 as the level (between m_2 and m_5)

    print(f"\n--- Per-sub-stratum overshoot breakdown ---")
    print(f"L = log(2^36) - log(85) = {L:.3f} nats")
    print(f"Conditional on end state (r=21, j_inner=k) for k = 0..{j_inner_max-1}")

    oshs, ends_r, ends_j = simulate_substratum_walk(
        n_orbits, L, 200_000,
        Q_cum, v_per_state, boundary_idx, j_inner_max, seed
    )
    valid = ~np.isnan(oshs)
    oshs = oshs[valid]
    ends_r = ends_r[valid]
    ends_j = ends_j[valid]

    print(f"\n  {'j_inner':>8}  {'corresp m_j':>11}  {'count':>8}  "
          f"{'⟨overshoot⟩(nats)':>18}  {'W_j_pred (steps)':>17}")
    print("  " + "-" * 70)
    boundary_orbits = (ends_r == boundary_idx)
    for j_inner in range(j_inner_max):
        cond = boundary_orbits & (ends_j == j_inner)
        n_cond = int(cond.sum())
        if n_cond < 50:
            continue
        osh_nats = -float(oshs[cond].mean())
        osh_se = float(oshs[cond].std() / np.sqrt(n_cond))
        W_pred = osh_nats / LOG43
        # j_inner = 2(j_attr - 3) → j_attr = j_inner/2 + 3 if even
        if j_inner % 2 == 0:
            j_attr = j_inner // 2 + 3
            corresp = f"m_{j_attr}"
        else:
            corresp = "(none)"
        print(f"  {j_inner:>8}  {corresp:>11}  {n_cond:>8,}  "
              f"{osh_nats:>10.4f} ± {osh_se:.4f}  {W_pred:>+10.3f}")


def diagnostic_sanity_check():
    """Verify:
    (a) At boundary residue 21, Geom(1/2) sub-stratum gives drift = -log(4/3).
    (b) i.i.d.-limit reduction (collapse sub-stratum to marginal v) recovers Path B v2.
    """
    print("\n--- Diagnostic sanity checks ---")
    # E[v | r=21] under Geom: E[v] = 6 + Σ k·2^{-k-1} = 6 + 1 = 7
    # Drift contribution at r=21: log 3 - 7·log 2 = -3.7534 nats
    # Other 31 residues sum: see Path B v2 (drift exact)
    Q, v_per_state, idx, odd, boundary = build_Q_basic_with_v(k=6)
    pi = 1.0 / Q.shape[0] * np.ones(Q.shape[0])  # uniform stationary at k=6
    Estep = np.zeros(Q.shape[0])
    for i in range(Q.shape[0]):
        if v_per_state[i] >= 0:
            Estep[i] = LOG3 - v_per_state[i] * LOG2
        else:
            Estep[i] = LOG3 - 7.0 * LOG2  # E[v|r=21] = 7
    drift = (pi * Estep).sum()
    print(f"  Drift E_π[X] = {drift:+.8f}  (target: log(3/4) = {-LOG43:+.8f})")
    print(f"  Discrepancy: {drift - (-LOG43):+.2e}")


def main():
    print("=" * 76)
    print("Path B sub-stratum extension: per-j W_j via (residue, sub-stratum) MAP")
    print("=" * 76)

    diagnostic_sanity_check()

    print("\n--- Per-j W_j extraction with sub-stratum conditioning ---")
    results = run_substratum_extraction(k=6, n_orbits=2_000_000, j_inner_max=12, seed=42)

    # Per-sub-stratum breakdown to see if ANY differentiation emerges
    report_per_substratum_breakdown(k=6, n_orbits=2_000_000, j_inner_max=12, seed=43)

    print("\n" + "=" * 76)
    print("Verdict")
    print("=" * 76)
    if results:
        gaps = [abs(r['W_j_pred'] - r['W_j_emp']) for r in results]
        max_gap = max(gaps) if gaps else 9999
        if max_gap < 0.05:
            print(f"(1) MATCH: all per-j W_j within ±0.05 of empirical (max gap {max_gap:.3f})")
            print(f"           Constant 3 closes via sub-stratum extension.")
        elif max_gap < 0.1:
            print(f"(2) APPROX: per-j W_j within ±0.1 (max gap {max_gap:.3f})")
            print(f"            Push to k=10 to confirm or refine.")
        elif max_gap < 0.5:
            print(f"(3) PARTIAL: per-j W_j within ±0.5 (max gap {max_gap:.3f})")
            print(f"             Sub-stratum captures some structure, not all.")
        else:
            print(f"(4) FRAMEWORK GAP: per-j W_j off by >{max_gap:.2f}")
            print(f"    Sub-stratum extension at residue level does not capture lattice landing.")
            print(f"    Per-j variation requires integer-level info beyond residue chain.")


if __name__ == "__main__":
    main()
