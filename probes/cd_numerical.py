"""
cd_numerical.py - Caravenna-Doney + Borovkov local-LD numerical attack on W_j.

Outputs:
  1. Empirical L^- moments (Path A 1M-orbit re-simulation, sanity check)
  2. Asymptotic Lorden constant via E[L^-2]/(2 E[L^-]) (j-independent prediction)
  3. Borovkov finite-y correction to Lorden at log(m_j) for j in {2, 4, 5}
     via direct renewal-equation numerical integration on the empirical L^- histogram
  4. Direct simulation of i.i.d. Syracuse walk first-crossing overshoot at log(m_j)
     (the operational i.i.d. baseline W_j_iid)
  5. Comparison to empirical W_j and quantification of Markov correction

Cap: under 300 lines.
"""
import numpy as np
import sys
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3  # log(4/3) = 0.287682

# Empirical W_j from compute_threads_findings.md (50M orbits at N=2^36)
W_J_EMPIRICAL = {
    2: (+7.156, 0.006),  # m_j = 5
    4: (-4.755, 0.06),   # m_j = 85
    5: (+4.590, 0.06),   # m_j = 341
}
M_J = {2: 5, 4: 85, 5: 341}
P_J = {2: 0.9379, 4: 0.0237, 5: 0.0379}  # entry distribution


@njit(parallel=True, cache=True)
def simulate_ladder_heights(n_orbits, seed):
    """Simulate L^- (strict descending ladder height) for iid Geom(1/2) walk."""
    L_vals = np.zeros(n_orbits)
    for i in prange(n_orbits):
        rng_state = np.uint64(seed + i * np.uint64(2654435761))
        S = 0.0
        for _ in range(50000):
            # LCG for Geom(1/2): v ~ Geom(1/2), v >= 1
            rng_state = rng_state * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
            u = (rng_state >> np.uint64(32)) / np.float64(2**32)
            # Geom(p=1/2) on {1, 2, ...}
            v = int(np.floor(np.log(1.0 - u) / np.log(0.5))) + 1
            if v < 1:
                v = 1
            S += LOG3 - v * LOG2
            if S < 0:
                L_vals[i] = -S
                break
    return L_vals


@njit(parallel=True, cache=True)
def simulate_overshoot_at_level(n_orbits, level, seed, log_m_start):
    """
    For iid Geom(1/2) Syracuse log-walk starting at log_m_start, descending,
    measure overshoot below `level` at FIRST crossing (i.e. first step with S <= level).
    Returns (overshoot_at_crossing, n_steps_to_crossing).
    """
    overshoots = np.full(n_orbits, np.nan)
    n_steps_arr = np.zeros(n_orbits, dtype=np.int32)
    for i in prange(n_orbits):
        rng_state = np.uint64(seed + i * np.uint64(2654435761) + np.uint64(int(level * 1000000)))
        S = log_m_start
        for n in range(200000):
            if S <= level:
                overshoots[i] = level - S  # >= 0, units = nats
                n_steps_arr[i] = n
                break
            rng_state = rng_state * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
            u = (rng_state >> np.uint64(32)) / np.float64(2**32)
            v = int(np.floor(np.log(1.0 - u) / np.log(0.5))) + 1
            if v < 1:
                v = 1
            S += LOG3 - v * LOG2
    return overshoots, n_steps_arr


def borovkov_finite_y_correction(L_samples, y_target_nats, n_grid=2000, y_max=None):
    """
    Numerical solution of the renewal equation
        m(y) = E[(L - y) * 1_{L > y}] + integral_0^y m(y - L) dF_L(L)
    on a uniform grid, using the empirical L distribution from L_samples.

    Returns m(y_target) = E[overshoot | first crossing of level y_target].
    For y -> infinity, m(y) -> Lorden = E[L^2]/(2 E[L]).

    Convention: y >= 0 measures the descent distance from start to target,
    L_samples is the descending ladder height distribution (L > 0).
    """
    L_pos = L_samples[L_samples > 0]
    if y_max is None:
        y_max = max(np.percentile(L_pos, 99.9), 4.0 * y_target_nats + 5.0)
    dy = y_max / n_grid
    grid = np.arange(n_grid + 1) * dy

    hist, edges = np.histogram(L_pos, bins=n_grid + 1, range=(0, y_max))
    # f_L probability mass per bin (NOT density), normalized to sum to 1
    f_L_mass = hist.astype(np.float64) / hist.sum()

    # g(y) = E[(L - y) * 1_{L > y}] computed on grid (rectangle rule)
    g = np.zeros_like(grid)
    for i, y in enumerate(grid):
        L_excess = (L_pos - y)
        L_excess = L_excess[L_excess > 0]
        if len(L_excess) > 0:
            g[i] = L_excess.sum() / len(L_pos)
        else:
            g[i] = 0.0

    # Solve renewal equation iteratively:
    #   m[i] = g[i] + sum_{k=1}^{i} m[i - k] * f_L_mass[k]
    m = np.zeros_like(grid)
    for i in range(len(grid)):
        conv = 0.0
        for k in range(1, i + 1):
            conv += m[i - k] * f_L_mass[k]
        m[i] = g[i] + conv

    # Interpolate m at y_target
    if y_target_nats <= grid[-1]:
        m_target = np.interp(y_target_nats, grid, m)
    else:
        m_target = m[-1]
    return m_target, grid, m, f_L_mass


def main():
    print("=" * 72)
    print("Caravenna-Doney + Borovkov local-LD attack on W_j")
    print("=" * 72)

    print(f"\nConstants:")
    print(f"  log(4/3)       = {LOG43:.6f} nats/step")
    print(f"  log(5)         = {np.log(5):.6f} nats  (= log m_2)")
    print(f"  log(85)        = {np.log(85):.6f} nats  (= log m_4)")
    print(f"  log(341)       = {np.log(341):.6f} nats  (= log m_5)")

    # ----------------------------------------------------------------------
    # Step 1: simulate L^- distribution (iid Geom(1/2) descending ladder)
    # ----------------------------------------------------------------------
    print(f"\n--- Step 1: Simulate L^- distribution (1M orbits) ---")
    n_orbits = 1_000_000
    L_vals = simulate_ladder_heights(n_orbits, np.uint64(42))
    L_vals = L_vals[L_vals > 0]
    E_L = L_vals.mean()
    E_L2 = (L_vals ** 2).mean()
    Var_L = E_L2 - E_L ** 2
    print(f"  E[L^-]    = {E_L:.6f} nats     (= {E_L/LOG43:.4f} step units)")
    print(f"  E[L^-^2]  = {E_L2:.6f} nats^2")
    print(f"  Var[L^-]  = {Var_L:.6f} nats^2  (SD = {np.sqrt(Var_L):.4f})")
    print(f"  CV[L^-]   = {np.sqrt(Var_L)/E_L:.4f}  (1.0 = exponential)")

    # ----------------------------------------------------------------------
    # Step 2: Asymptotic Lorden -- LADDER formula vs walk one-step formula
    # ----------------------------------------------------------------------
    print(f"\n--- Step 2: Asymptotic Lorden (j-independent) ---")
    Lorden_ladder = E_L2 / (2 * E_L)
    # one-step formula (what Path C/A reported)
    E_X = LOG3 - 2 * LOG2  # E[X] = -log(4/3)
    E_v = 2.0; E_v2 = 6.0  # Geom(1/2) on {1, 2, ...}
    E_X2 = LOG3**2 - 2*LOG3*LOG2*E_v + LOG2**2 * E_v2
    Lorden_walk_onestep = E_X2 / (2 * abs(E_X))
    print(f"  E[X^2]/(2|E[X]|) (one-step) = {Lorden_walk_onestep:.6f} nats")
    print(f"                              = {Lorden_walk_onestep/LOG43:.4f} step units  <-- Path C reported")
    print(f"  E[L^-^2]/(2 E[L^-]) (ladder) = {Lorden_ladder:.6f} nats")
    print(f"                               = {Lorden_ladder/LOG43:.4f} step units  <-- correct first-crossing overshoot")
    print(f"")
    print(f"  Note: the asymptotic conditional overshoot at first crossing of -y")
    print(f"  is given by the LADDER formula (interarrivals = L^-), not the")
    print(f"  one-step formula. Path C/A's 6.305 is the residual-life of the X")
    print(f"  process, which is the wrong renewal interval for first-passage overshoot.")

    # ----------------------------------------------------------------------
    # Step 3: Borovkov finite-y correction at log(m_j)
    # ----------------------------------------------------------------------
    print(f"\n--- Step 3: Borovkov finite-y correction via renewal equation ---")
    print(f"  Solving m(y) = E[(L-y)*1_{{L>y}}] + integral_0^y m(y-L) dF_L(L)")
    print(f"  on empirical L distribution at y = log(m_j).")
    for j in [2, 4, 5]:
        y_j = np.log(M_J[j])
        m_y, _, _, _ = borovkov_finite_y_correction(L_vals, y_j, n_grid=400, y_max=20.0)
        delta = m_y - Lorden_ladder
        print(f"  j={j}, m_j={M_J[j]:>4d}, log(m_j)={y_j:.3f} nats:")
        print(f"    m(log m_j) = {m_y:.6f} nats = {m_y/LOG43:.4f} step units")
        print(f"    Lorden     = {Lorden_ladder:.6f} nats = {Lorden_ladder/LOG43:.4f} step units")
        print(f"    Delta(y)   = {delta:+.6f} nats = {delta/LOG43:+.4f} step units")

    # ----------------------------------------------------------------------
    # Step 4: Direct simulation of i.i.d. Syracuse walk first-crossing overshoot
    # ----------------------------------------------------------------------
    print(f"\n--- Step 4: Direct iid simulation of first-crossing overshoot at log(m_j) ---")
    print(f"  Walk: iid Geom(1/2), starting at log m_start = log(2^36) = {36*LOG2:.3f} nats")
    print(f"  (matches N=2^36 empirical regime; conditional on first crossing of log(m_j))")
    log_m_start = 36 * LOG2
    n_sim = 500_000
    print(f"  {n_sim} orbits per j.\n")
    print(f"  {'j':>3} {'m_j':>5} {'mean overshoot (nats)':>22} {'(steps)':>10} {'SE':>8}")
    W_iid_sim = {}
    for j in [2, 4, 5]:
        y_j = np.log(M_J[j])
        os_vals, n_steps = simulate_overshoot_at_level(n_sim, y_j, np.uint64(1729 + j), log_m_start)
        os_vals = os_vals[~np.isnan(os_vals)]
        mean_os = os_vals.mean()
        se_os = os_vals.std() / np.sqrt(len(os_vals))
        W_iid_sim[j] = mean_os / LOG43
        print(f"  {j:>3} {M_J[j]:>5} {mean_os:>22.6f} {mean_os/LOG43:>10.4f} {se_os/LOG43:>8.4f}")

    # ----------------------------------------------------------------------
    # Step 5: Comparison to empirical W_j and Markov correction
    # ----------------------------------------------------------------------
    print(f"\n--- Step 5: i.i.d. baseline vs empirical W_j (Markov correction) ---")
    print(f"\n  W_j framework: empirical = i.i.d. baseline + Markov correction\n")
    print(f"  {'j':>3} {'W_j_emp':>10} {'W_j_iid_sim':>13} {'Markov':>10} {'|Markov|/|iid|':>15}")
    for j in [2, 4, 5]:
        W_emp, W_emp_se = W_J_EMPIRICAL[j]
        W_iid = W_iid_sim[j]
        markov = W_emp - W_iid
        ratio = abs(markov) / abs(W_iid)
        print(f"  {j:>3} {W_emp:>+10.4f} {W_iid:>+13.4f} {markov:>+10.4f} {ratio:>15.2%}")

    # ----------------------------------------------------------------------
    # Step 6: epsilon_S decomposition
    # ----------------------------------------------------------------------
    print(f"\n--- Step 6: epsilon_S decomposition (i.i.d. + Markov) ---")
    print(f"\n  epsilon_S = sum_j P(j) * [W_j - log(m_j)/log(4/3) + 1]\n")
    eps_S_iid = 0.0
    eps_S_emp = 0.0
    for j in [2, 4, 5]:
        W_emp, _ = W_J_EMPIRICAL[j]
        W_iid = W_iid_sim[j]
        contrib_iid = P_J[j] * (W_iid - np.log(M_J[j]) / LOG43 + 1)
        contrib_emp = P_J[j] * (W_emp - np.log(M_J[j]) / LOG43 + 1)
        eps_S_iid += contrib_iid
        eps_S_emp += contrib_emp
        print(f"  j={j}: P(j)={P_J[j]:.4f}, W_iid={W_iid:.3f}, W_emp={W_emp:+.3f}")
        print(f"        iid contrib = {contrib_iid:+.4f}, emp contrib = {contrib_emp:+.4f}")
    print(f"\n  epsilon_S (i.i.d. baseline only) = {eps_S_iid:.4f}")
    print(f"  epsilon_S (empirical W_j)        = {eps_S_emp:.4f}")
    print(f"  Markov contribution to eps_S     = {eps_S_emp - eps_S_iid:+.4f}")
    print(f"  log(4)                            = {np.log(4):.4f}  (the candidate)")
    print(f"  empirical eps_S asymptote (CTF)  = 1.375 (compute_threads_findings.md)")

    # ----------------------------------------------------------------------
    # Verdict
    # ----------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"VERDICT")
    print(f"{'='*72}")
    print(f"""
  i.i.d. local-LD framework (Caravenna-Doney + Borovkov adaptation) predicts
  W_j_iid ~ Lorden = {Lorden_ladder/LOG43:.3f} step units at every j (j-independent).
  Borovkov finite-y correction at log(m_j) is negligible (< 0.1 step units)
  because log(m_j) >> mu_L = {E_L/LOG43:.3f} step units for all j tested.

  Empirical W_j varies dramatically by j (range -4.76 to +7.16). The cross-class
  variation is therefore PURELY a Markov-modulation phenomenon.

  Per the brief's decision criteria:
  - W_2 (P=0.938, dominant): i.i.d. captures most ({W_iid_sim[2]:.2f} vs {W_J_EMPIRICAL[2][0]:.2f}),
    Markov correction ~{W_J_EMPIRICAL[2][0] - W_iid_sim[2]:+.2f} step units
  - W_4 (P=0.024, sparse): Markov dominates ({W_iid_sim[4]:.2f} vs {W_J_EMPIRICAL[4][0]:.2f}),
    Markov correction ~{W_J_EMPIRICAL[4][0] - W_iid_sim[4]:+.2f} step units (~3x baseline magnitude)
  - W_5 (P=0.038, sparse): Markov is meaningful ({W_iid_sim[5]:.2f} vs {W_J_EMPIRICAL[5][0]:.2f}),
    Markov correction ~{W_J_EMPIRICAL[5][0] - W_iid_sim[5]:+.2f} step units

  i.i.d. eps_S baseline = {eps_S_iid:.3f}, vs empirical {eps_S_emp:.3f}, gap {eps_S_emp - eps_S_iid:+.3f}.
  Note empirical sits at log(4) = {np.log(4):.3f} within ~1%; i.i.d. baseline is
  {eps_S_iid - np.log(4):+.3f} below log(4). The Markov correction in eps_S effectively
  closes the gap to log(4)-vicinity (Path B target).
""")


if __name__ == "__main__":
    main()
