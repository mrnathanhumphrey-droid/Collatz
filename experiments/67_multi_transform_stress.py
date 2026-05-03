"""
Multi-transformation stress test of three Lagarias-class observable slices:
  1. w_q(q) — per-band Esscher tilt
  2. P(q|j) — conditional σ-band given absorbing attractor
  3. ⟨v|q,j⟩ — joint table

Operates on existing measurements; no new orbit generation.

Tests 1-5 per brief.
"""
import sys
import io
import numpy as np
from scipy.stats import norm
from scipy.special import logit
import polars as pl
from pathlib import Path

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG2 = np.log(2.0); LOG3 = np.log(3.0)
out_dir = Path("C:/Collatz/experiments_output")

# ============================================================
# DATA from existing Results
# ============================================================

# Result 26 / exp 61: 10-quantile w_q at N=2^36
W_Q_DATA = [
    # (q, z_q, E[v], w_q)
    (0.050, -1.645, 2.4976, -0.2621),
    (0.125, -1.150, 2.2730, -0.1637),
    (0.200, -0.842, 2.1876, -0.1187),
    (0.375, -0.319, 2.0710, -0.0486),
    (0.500,  0.000, 2.0115, -0.0082),
    (0.625, +0.319, 1.9657, +0.0254),
    (0.800, +0.842, 1.9043, +0.0744),
    (0.875, +1.150, 1.8747, +0.0998),
    (0.950, +1.645, 1.8290, +0.1416),
    (0.975, +1.960, 1.8079, +0.1620),
]

# Result 33 P(q|j) at N=2^32 (from P_q_given_j_log.txt and joint_qj_and_prime_test.md)
# 4 q-bands × 3 j-classes. P(q) is uniform = 0.25 (quartile bands).
P_Q_GIVEN_J = {
    # q : { j : P(q|j) }
    0.125: {2: 0.2323, 4: 0.5512, 5: 0.4724},
    0.375: {2: 0.2393, 4: 0.2517, 5: 0.2642},
    0.625: {2: 0.2616, 4: 0.1320, 5: 0.1678},
    0.875: {2: 0.2667, 4: 0.0651, 5: 0.0957},
}
P_Q_MARGINAL = {0.125: 0.25, 0.375: 0.25, 0.625: 0.25, 0.875: 0.25}

# joint ⟨v|q,j⟩ table (Result 33) at N=2^32
V_Q_J = {
    # q : { j : ⟨v|q,j⟩ }
    0.125: {2: 2.356, 4: 2.439, 5: 2.409},
    0.375: {2: 2.078, 4: 2.087, 5: 2.084},
    0.625: {2: 1.966, 4: 1.973, 5: 1.970},
    0.875: {2: 1.865, 4: 1.873, 5: 1.872},
}

# ⟨v|q⟩ marginal (bulk values from joint table — joint_qj_and_prime_test.md)
V_Q_MARGINAL = {0.125: 2.365, 0.375: 2.079, 0.625: 1.966, 0.875: 1.865}

# Per-j observables (W_j, ⟨v|j⟩, ⟨σ_S|j⟩) — from prime_vs_all log
PER_J = {
    # j : (W_j, ⟨v|j⟩, ⟨σ_S|j⟩)
    2: (+7.141, 2.057, 76.17),
    4: (-4.679, 2.251, 54.49),
    5: (+4.638, 2.199, 58.98),
}

# Result 20 Markov correction (j-specific step units)
# j=2: +3.69, j=4: −8.20, j=5: +1.14
MARKOV_CORR = {2: +3.69, 4: -8.20, 5: +1.14}

# log(m_j) = log((4^j - 1)/3)
def log_m_j(j):
    return float(np.log((4**j - 1) / 3))

# ============================================================
# TEST 1: KL divergences P(q|j) || P(q) per j
# ============================================================
def test1_kl():
    print("\n" + "="*70, flush=True)
    print("TEST 1: KL divergences P(q|j) || P(q) and structural correlations", flush=True)
    print("="*70, flush=True)

    j_list = sorted(PER_J.keys())
    q_list = sorted(P_Q_GIVEN_J.keys())

    rows = []
    for j in j_list:
        # D_KL(P(q|j) || P(q))
        d_kl = 0.0
        d_chi2 = 0.0
        d_hell = 0.0
        for q in q_list:
            p_qj = P_Q_GIVEN_J[q][j]
            p_q = P_Q_MARGINAL[q]
            if p_qj > 0:
                d_kl += p_qj * np.log(p_qj / p_q)
            d_chi2 += (p_qj - p_q)**2 / p_q
            d_hell += (np.sqrt(p_qj) - np.sqrt(p_q))**2
        d_hell = 0.5 * d_hell  # squared Hellinger

        W_j, v_j, sigma_S_j = PER_J[j]
        markov = MARKOV_CORR[j]
        log_mj = log_m_j(j)

        rows.append({
            'j': j, 'D_KL': d_kl, 'D_chi2': d_chi2, 'D_Hellinger': d_hell,
            'W_j': W_j, 'v_j': v_j, 'sigma_S_j': sigma_S_j, 'markov_corr': markov,
            'log_m_j': log_mj,
        })

    print(f"\n  {'j':>3}  {'D_KL':>9}  {'D_chi2':>9}  {'D_Hell':>9}  {'W_j':>9}  "
          f"{'⟨v|j⟩':>7}  {'⟨σ|j⟩':>9}  {'Markov':>8}  {'log(m_j)':>10}", flush=True)
    for r in rows:
        print(f"  {r['j']:>3}  {r['D_KL']:>9.4f}  {r['D_chi2']:>9.4f}  "
              f"{r['D_Hellinger']:>9.4f}  {r['W_j']:>+9.3f}  {r['v_j']:>7.3f}  "
              f"{r['sigma_S_j']:>9.2f}  {r['markov_corr']:>+8.3f}  {r['log_m_j']:>10.4f}", flush=True)

    # Correlations
    arrs = {k: np.array([r[k] for r in rows]) for k in
            ['D_KL', 'D_chi2', 'D_Hellinger', 'W_j', 'v_j', 'sigma_S_j', 'markov_corr', 'log_m_j']}

    print(f"\n  Pearson correlations of D_KL with j-specific observables:", flush=True)
    for k in ['W_j', 'v_j', 'sigma_S_j', 'markov_corr', 'log_m_j']:
        if arrs[k].std() > 1e-9:
            r = np.corrcoef(arrs['D_KL'], arrs[k])[0,1]
            print(f"    D_KL vs {k:>15}: r = {r:+.4f}", flush=True)

    # Note: only 3 j values, so correlations are essentially perfect with one degree of freedom
    print(f"\n  NOTE: only 3 j classes → 1 df. Perfect correlation possible by chance.", flush=True)
    print(f"        Report direction and magnitude, not significance.", flush=True)

    return rows


# ============================================================
# TEST 2: Log-ratio structure across j
# ============================================================
def test2_log_ratios():
    print("\n" + "="*70, flush=True)
    print("TEST 2: log P(q|j_a)/P(q|j_b) functional forms in q", flush=True)
    print("="*70, flush=True)

    q_list = sorted(P_Q_GIVEN_J.keys())
    q_arr = np.array(q_list)
    z_arr = np.array([norm.ppf(q) for q in q_list])
    logit_arr = np.array([np.log(q/(1-q)) for q in q_list])

    pairs = [(2,4), (4,5), (2,5)]
    rows = []
    for ja, jb in pairs:
        log_ratio = np.array([np.log(P_Q_GIVEN_J[q][ja] / P_Q_GIVEN_J[q][jb]) for q in q_list])
        print(f"\n  pair (j={ja}, j={jb}):", flush=True)
        for label, X in [('q', q_arr), ('z_q', z_arr), ('logit(q)', logit_arr)]:
            # Linear fit
            Xc = X - X.mean(); Yc = log_ratio - log_ratio.mean()
            slope = float((Xc*Yc).sum() / (Xc*Xc).sum())
            intercept = float(log_ratio.mean() - slope*X.mean())
            pred = intercept + slope*X
            ss_res = float(((log_ratio - pred)**2).sum())
            ss_tot = float(((log_ratio - log_ratio.mean())**2).sum())
            R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 1.0
            print(f"    fit log_ratio = {intercept:+.4f} + {slope:+.4f}·{label:>8}    R² = {R2:.4f}", flush=True)
            rows.append({'pair': f'{ja}_{jb}', 'transform': label, 'slope': slope,
                         'intercept': intercept, 'R2': R2})
        print(f"    raw values: {log_ratio.tolist()}", flush=True)

    return rows


# ============================================================
# TEST 3: Conserved quantity search across (q, j) cells
# ============================================================
def test3_invariants():
    print("\n" + "="*70, flush=True)
    print("TEST 3: Conserved quantity search across (q,j) cells", flush=True)
    print("="*70, flush=True)

    q_list = sorted(P_Q_GIVEN_J.keys())
    j_list = sorted(PER_J.keys())

    # Build E_band(q) lookup from W_Q_DATA at the 4 q-bands
    E_band_lookup = {}
    for q, z, Ev, w in W_Q_DATA:
        E_band_lookup[q] = Ev
    # 0.625 not in W_Q_DATA exactly but 0.625 IS in. Let me check:
    # W_Q_DATA has 0.05, 0.125, 0.20, 0.375, 0.50, 0.625, 0.80, 0.875, 0.95, 0.975
    # So 0.125, 0.375, 0.625, 0.875 all present.

    # w_q lookup at the 4 q-bands
    w_q_lookup = {q: w for q, z, Ev, w in W_Q_DATA}

    cells = []
    for q in q_list:
        for j in j_list:
            cells.append({
                'q': q, 'j': j,
                'E_band': E_band_lookup[q],
                'w_q': w_q_lookup[q],
                'P_qj': P_Q_GIVEN_J[q][j],
                'P_q': P_Q_MARGINAL[q],
                'v_qj': V_Q_J[q][j],
                'log_mj': log_m_j(j),
                'W_j': PER_J[j][0],
                'v_j': PER_J[j][1],
                'sigma_S_j': PER_J[j][2],
            })

    # 12 cells total
    candidates = {
        # raw
        'E_band·P(q|j)': lambda c: c['E_band'] * c['P_qj'],
        'w_q·log(m_j)': lambda c: c['w_q'] * c['log_mj'],
        'v_qj·sigma_S_j (per cell)': lambda c: c['v_qj'] * c['sigma_S_j'],
        'P(q|j)·W_j': lambda c: c['P_qj'] * c['W_j'],
        # sqrt
        'sqrt(E_band)·sqrt(P(q|j))': lambda c: np.sqrt(c['E_band'] * c['P_qj']),
        'sqrt(w_q^2 + log^2(m_j))': lambda c: np.sqrt(c['w_q']**2 + c['log_mj']**2),
        'sqrt(P(q|j)·P(q))': lambda c: np.sqrt(c['P_qj'] * c['P_q']),
        # log
        'log(E_band)−log(P(q|j))': lambda c: np.log(c['E_band']) - np.log(c['P_qj']),
        'log|w_q|−log(m_j)': lambda c: np.log(abs(c['w_q'])) - np.log(c['log_mj']),
        # E[v|q,j] · P(q|j)/P(q)
        'v_qj·P(q|j)/P(q)': lambda c: c['v_qj'] * c['P_qj'] / c['P_q'],
        # ratios
        'v_qj/E_band': lambda c: c['v_qj'] / c['E_band'],
        'P(q|j)/E_band': lambda c: c['P_qj'] / c['E_band'],
    }

    print(f"\n  Compute coefficient of variation (SD/|mean|) across {len(cells)} (q,j) cells:", flush=True)
    print(f"  Lower CV = more invariant. CV < 0.02 = candidate structural invariant.", flush=True)
    print(f"\n  {'candidate':<35}  {'mean':>10}  {'SD':>10}  {'CV':>8}", flush=True)

    rows = []
    for name, fn in candidates.items():
        vals = np.array([fn(c) for c in cells])
        mean = float(vals.mean())
        sd = float(vals.std())
        cv = sd / abs(mean) if abs(mean) > 1e-9 else float('inf')
        print(f"  {name:<35}  {mean:>+10.4f}  {sd:>10.4f}  {cv:>8.4f}", flush=True)
        rows.append({'candidate': name, 'mean': mean, 'sd': sd, 'cv': cv})

    # Best (lowest) CV
    rows_sorted = sorted(rows, key=lambda r: r['cv'])
    print(f"\n  Top 3 lowest CV:", flush=True)
    for r in rows_sorted[:3]:
        print(f"    {r['candidate']:<35}  CV = {r['cv']:.4f}", flush=True)

    # Verdict
    if rows_sorted[0]['cv'] < 0.02:
        print(f"\n  → STRUCTURAL INVARIANT: {rows_sorted[0]['candidate']} (CV={rows_sorted[0]['cv']:.4f})", flush=True)
    elif rows_sorted[0]['cv'] < 0.05:
        print(f"\n  → MARGINAL: lowest CV = {rows_sorted[0]['cv']:.4f}", flush=True)
    else:
        print(f"\n  → NULL: no clean invariant; lowest CV = {rows_sorted[0]['cv']:.4f}", flush=True)

    return rows


# ============================================================
# TEST 4: Residual structure under different transformations of w_q
# ============================================================
def test4_w_q_transforms():
    print("\n" + "="*70, flush=True)
    print("TEST 4: w_q(q) under different transformations", flush=True)
    print("="*70, flush=True)

    # Filter out q=0.5 where w_q ≈ 0 and z_q = 0 → ratios diverge
    data = [(q, z, Ev, w) for q, z, Ev, w in W_Q_DATA if abs(z) > 0.01]
    q = np.array([d[0] for d in data])
    z = np.array([d[1] for d in data])
    w = np.array([d[3] for d in data])

    transforms = {
        'raw: w vs z':           (z, w),
        'log: w vs log|z|·sgn':  (np.sign(z) * np.log(np.abs(z)+1e-9), w),
        'sqrt: w vs sgn·sqrt|z|': (np.sign(z) * np.sqrt(np.abs(z)), w),
        'logit: w vs log(q/(1-q))': (np.log(q/(1-q)), w),
        'inverse: log|w| vs log|z|': (np.log(np.abs(z)), np.log(np.abs(w))),
        'w/z ratio asymm: |w|/|z| vs sgn(z)': (np.sign(z), np.abs(w)/np.abs(z)),
    }

    rows = []
    for label, (X, Y) in transforms.items():
        Xc = X - X.mean(); Yc = Y - Y.mean()
        if (Xc*Xc).sum() > 1e-15:
            slope = float((Xc*Yc).sum() / (Xc*Xc).sum())
            intercept = float(Y.mean() - slope*X.mean())
            pred = intercept + slope*X
            ss_res = float(((Y - pred)**2).sum())
            ss_tot = float(((Y - Y.mean())**2).sum())
            R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 1.0
        else:
            slope, intercept, R2 = float('nan'), float('nan'), float('nan')
        print(f"  {label:<45}  slope={slope:+.4f}  R² = {R2:.4f}", flush=True)
        rows.append({'transform': label, 'slope': slope, 'intercept': intercept, 'R2': R2})

    # Asymmetry check: |w_q|/|z_q| split by sign(z)
    ratio = np.abs(w) / np.abs(z)
    lower = ratio[z < 0]
    upper = ratio[z > 0]
    print(f"\n  Asymmetry check: |w_q|/|z_q|", flush=True)
    print(f"    lower half (z<0): mean = {lower.mean():.4f}  range = [{lower.min():.4f}, {lower.max():.4f}]", flush=True)
    print(f"    upper half (z>0): mean = {upper.mean():.4f}  range = [{upper.min():.4f}, {upper.max():.4f}]", flush=True)
    print(f"    ratio (lower/upper means): {lower.mean()/upper.mean():.4f}", flush=True)
    print(f"    → asymmetry persists if ratio > 1.3 or < 0.77; raw is {lower.mean():.4f} vs {upper.mean():.4f}", flush=True)

    return rows


# ============================================================
# TEST 5: N-stability (we have w_q at multiple N from various exps)
# ============================================================
def test5_n_stability():
    print("\n" + "="*70, flush=True)
    print("TEST 5: N-stability of identified structures", flush=True)
    print("="*70, flush=True)

    # We don't have full P(q|j) at multiple N from existing data,
    # only at N=2^32. w_q at multiple N is also limited.
    # For the candidates that came up in Test 4 best transformation,
    # we can use the per-band E_band data from exp 60 etc.

    # Note: per Result 32, drift in joint moments occurs slowly with N.
    # We have w_q at 5 bands at N=2^36 (exp 60). And 10-band at N=2^36 (exp 61).
    # No multi-N w_q sweep at all 10 quantiles.

    # For the basic stability check: verify the (q,j) reduction from Result 33
    # (same q-distribution → reduction holds) at the q-bands where we have data.
    # Compare exp 33 N=2^32 P(q|j) max-z vs theoretical.

    print(f"\n  P(q|j) measurement at N=2^32 (Result 33 only N available):", flush=True)
    print(f"    max pairwise |z| across q ∈ [0.125, 0.875]: 211.6 (q=0.125), 30.85 (q=0.375),", flush=True)
    print(f"    38.23 (q=0.625), 59.05 (q=0.875). All > 30σ, all structural.", flush=True)
    print(f"    No multi-N P(q|j) data available; cannot verify N-stability of P(q|j) directly.", flush=True)

    print(f"\n  w_q stability across N (from exp 60 vs exp 61):", flush=True)
    # exp 60 5-band w_q at N=2^36: see joint table
    w_q_exp60 = {0.125: -0.1373, 0.375: -0.0214, 0.625: +0.0456, 0.875: +0.1577, 0.975: +0.2491}
    w_q_exp61 = {q: w for q, z, Ev, w in W_Q_DATA}
    for q in [0.125, 0.375, 0.625, 0.875, 0.975]:
        if q in w_q_exp60 and q in w_q_exp61:
            print(f"    q={q}: exp60 w_q={w_q_exp60[q]:+.4f}  exp61 w_q={w_q_exp61[q]:+.4f}  "
                  f"gap={w_q_exp60[q]-w_q_exp61[q]:+.4f}", flush=True)

    return []


def main():
    print("="*70, flush=True)
    print("Multi-transformation stress test of three Lagarias-class slices", flush=True)
    print("="*70, flush=True)
    print(f"\nData: w_q (Result 26, exp 61, N=2^36, 10 quantiles)", flush=True)
    print(f"      P(q|j) (Result 33, P_q_given_j.py, N=2^32, 4 q × 3 j)", flush=True)
    print(f"      ⟨v|q,j⟩ (Result 33 joint table, N=2^32)", flush=True)
    print(f"      W_j, ⟨v|j⟩, ⟨σ_S|j⟩ (Result 30, prime_vs_all)", flush=True)
    print(f"      Markov correction (Result 20)", flush=True)

    r1 = test1_kl()
    r2 = test2_log_ratios()
    r3 = test3_invariants()
    r4 = test4_w_q_transforms()
    r5 = test5_n_stability()

    # Save
    pl.DataFrame(r1).write_csv(out_dir / "67_test1_kl.csv")
    pl.DataFrame(r2).write_csv(out_dir / "67_test2_logratios.csv")
    pl.DataFrame(r3).write_csv(out_dir / "67_test3_invariants.csv")
    pl.DataFrame(r4).write_csv(out_dir / "67_test4_transforms.csv")
    print(f"\n[save] CSVs written to {out_dir}", flush=True)


if __name__ == "__main__":
    main()
