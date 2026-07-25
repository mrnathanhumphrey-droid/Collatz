"""
GATE SPLITTING (S2) -- the DEEP run: bigint accumulator + ESS/level monitor, past the int64 wall.

S1 validated the estimator at r<=16 (gamma to 0.1%, unbiased). Two builds were owed before going deep:
  (1) int64 vectorization walls at r~18 (modmul mod 3^{r+1} overflows once M^2 > 2^63). -> Python-bigint
      via numpy object arrays (arbitrary precision, no overflow).
  (2) particle degeneracy over 30+ levels is UNTESTED by the r<=16 gate. -> monitor+report ESS/level.

STAGE A (validation, fast): run the bigint estimator AND the int64 estimator on the SAME seed at r<=16.
  With no overflow at r<=16 they must give BIT-IDENTICAL gamma. That certifies the bigint path
  (S1 already chained build_nu ~ int64; so build_nu ~ int64 == bigint => bigint correct). No build_nu needed.

STAGE B (deep, ESS-monitored): run bigint splitting to R_DEEP, moderate N, R_REP replicas, channels m in MLIST.
  Report per level: survival q_s, survivors k_s, and unique-state fraction (ESS proxy) after resample.
  Reconstruct the TELESCOPED cumulative:  eps_hat_r = eps16 + 2*sum_m 4^-m (gamma_r(m) - gamma_16(m)),
  per replica, then mean +/- SE across replicas. Question: is eps_hat still RISING at R_DEEP, or has it TURNED?

  This is a SCOPED first pass (R_DEEP=30 brackets the predicted crossing #3 interval [27,31]); the decisive
  >=3sigma decision run (to r=35, N~5e5) is a second step, gated on ESS holding up here.

EPISTEMIC GUARD: the output is a STATISTICAL object with an error budget, NEVER welded to the exact eps ladder.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

# ---------------- config ----------------
R_VALID = 16          # stage A: bigint-vs-int64 identical-seed check
R_DEEP  = 40          # stage B: deep target -- a CURVE to 40 (shape is the observable, not the endpoint)
CURVE_R = [20, 25, 28, 31, 35, 40]   # reported curve points (telescoped gamma-differences)
MLIST   = [1, 2]      # channels (turnover dominated by m=1; m=2 is a 6% check)
N       = 400000      # particles per replica per channel
R_REP   = 20          # replicas -> SE
JMAX    = 80          # clip geometric (P[v>80] ~ 1e-24)
SEED    = 20260723
FP_MOD  = (1 << 61) - 1               # int64 fingerprint for fast uniq-frac (collision prob ~ N^2/2^61 negl.)

EPS16   = 4.685e-3    # exact/nu-validated eps_16 anchor (result_logperiodic.md chain)
# ----------------------------------------

MDEEP   = 3 ** (R_DEEP + 1)                     # Python int modulus (bigint)
INV2_D  = pow(2, -1, MDEEP)
INV2POW_D = np.array([pow(INV2_D, j, MDEEP) for j in range(JMAX + 1)], dtype=object)
POW3_D    = [3 ** s for s in range(R_DEEP + 2)]  # Python ints

# int64 tables for the r<=16 validation
MV      = 3 ** (R_VALID + 1)
INV2_V  = pow(2, -1, MV)
INV2POW_V = np.array([pow(INV2_V, j, MV) for j in range(JMAX + 1)], dtype=np.int64)
POW3_V    = [3 ** s for s in range(R_VALID + 2)]


def run_replica(rng, R, mlist, N, M, INV2POW, POW3, dtype, ess=False):
    """One replica. Returns {m: gamma[0..R]} and (if ess) {'q','k','uniq'} per level (m=mlist[0])."""
    out = {}
    diag = {'q': np.zeros(R + 1), 'k': np.zeros(R + 1, dtype=int), 'uniq': np.zeros(R + 1)}
    for mi, m in enumerate(mlist):
        fourm = 4 ** m                                   # test: (4^m accX - accP) == 0 mod 3^{s+1}
        accX = np.ones(N, dtype=dtype); wX = np.ones(N, dtype=dtype)
        accP = np.ones(N, dtype=dtype); wP = np.ones(N, dtype=dtype)
        # re-seed identically per channel so channels are independent draws but reproducible
        rloc = np.random.default_rng(rng.integers(0, 2**63 - 1))
        gamma = np.zeros(R + 1); p = 1.0; dead = False
        for s in range(1, R + 1):
            if dead:
                gamma[s] = 0.0; continue
            vX = np.clip(rloc.geometric(0.5, N), 1, JMAX)
            vP = np.clip(rloc.geometric(0.5, N), 1, JMAX)
            wX = (wX * INV2POW[vX]) % M
            wP = (wP * INV2POW[vP]) % M
            t3 = POW3[s]
            accX = (accX + (t3 * wX) % M) % M
            accP = (accP + (t3 * wP) % M) % M
            D = (fourm * accX - accP) % M
            survive = (D % POW3[s + 1]) == 0
            k = int(survive.sum())
            if k == 0:
                dead = True; gamma[s] = 0.0; continue
            p *= k / N
            gamma[s] = POW3[s] * p
            idx = np.nonzero(survive)[0]
            pick = idx[rloc.integers(0, k, N)]
            accX = accX[pick]; wX = wX[pick]; accP = accP[pick]; wP = wP[pick]
            if ess and mi == 0:
                diag['q'][s] = k / N; diag['k'][s] = k
                fp = (accX % FP_MOD).astype(np.int64)       # fast int64 fingerprint of the state
                diag['uniq'][s] = len(np.unique(fp)) / N
        out[m] = gamma
    return (out, diag) if ess else out


def main():
    t0 = time.time()
    print(f"# GATE SPLITTING (S2) deep  R_DEEP={R_DEEP} N={N} R_REP={R_REP} MLIST={MLIST}\n")

    # ---------- STAGE A: bigint == int64 at r<=16 (bit-identical on same seed) ----------
    print("## STAGE A -- bigint-vs-int64 identical-seed validation (r<=16)")
    sA = np.random.default_rng(SEED)
    seedcheck = sA.integers(0, 2**63 - 1)
    gi = run_replica(np.random.default_rng(seedcheck), R_VALID, MLIST, 20000, MV, INV2POW_V, POW3_V, np.int64)
    gb = run_replica(np.random.default_rng(seedcheck), R_VALID, MLIST, 20000, MDEEP, INV2POW_D, POW3_D, object)
    maxdiff = max(np.max(np.abs(gi[m] - gb[m])) for m in MLIST)
    print(f"   max|gamma_int64 - gamma_bigint| over r<=16, m={MLIST}: {maxdiff:.2e}")
    print(f"   => bigint path {'CERTIFIED (bit-identical to validated int64)' if maxdiff < 1e-9 else 'FAIL -- divergence!'}\n")
    if maxdiff >= 1e-9:
        print("   ABORT: bigint arithmetic disagrees with int64 -- not running deep."); return

    # ---------- STAGE B: deep ESS-monitored run ----------
    print(f"## STAGE B -- deep run to r={R_DEEP}, ESS monitored")
    ss = np.random.SeedSequence(SEED + 1)
    rngs = [np.random.default_rng(s) for s in ss.spawn(R_REP)]
    G = {m: np.zeros((R_REP, R_DEEP + 1)) for m in MLIST}
    uniq_acc = np.zeros((R_REP, R_DEEP + 1)); q_acc = np.zeros((R_REP, R_DEEP + 1))
    for i, rng in enumerate(rngs):
        rep, diag = run_replica(rng, R_DEEP, MLIST, N, MDEEP, INV2POW_D, POW3_D, object, ess=True)
        for m in MLIST:
            G[m][i] = rep[m]
        uniq_acc[i] = diag['uniq']; q_acc[i] = diag['q']
        print(f"   replica {i+1}/{R_REP}  ({time.time()-t0:.1f}s)")
    print()

    # ESS / degeneracy: FULL per-level curve through r=40 (guardrail 1)
    print("## ESS / DEGENERACY -- full per-level curve (m=%d; uniq-frac after resample; baseline ceiling q(1-e^-3)~0.317)" % MLIST[0])
    print(f"   {'r':>3} {'q_s':>8} {'uniq':>8}   {'r':>3} {'q_s':>8} {'uniq':>8}")
    uq = uniq_acc.mean(0); qq = q_acc.mean(0)
    rows = list(range(17, R_DEEP + 1))
    for a in range(0, len(rows), 2):
        s1 = rows[a]; line = f"   {s1:>3} {qq[s1]:>8.4f} {uq[s1]:>8.4f}"
        if a + 1 < len(rows):
            s2 = rows[a + 1]; line += f"   {s2:>3} {qq[s2]:>8.4f} {uq[s2]:>8.4f}"
        print(line)
    print(f"   uniq-frac: r17={uq[17]:.3f} -> r30={uq[30]:.3f} -> r40={uq[40]:.3f}  "
          f"[{'holds' if uq[40] > 0.10 else 'COLLAPSING -> bars widen, SNR revises DOWN'}]\n")

    # telescoped cumulative eps_hat per replica -- the CURVE (shape is the observable)
    eps_rep = np.zeros((R_REP, R_DEEP + 1))
    for i in range(R_REP):
        for r in range(16, R_DEEP + 1):
            dg = sum(4.0 ** -m * (G[m][i, r] - G[m][i, 16]) for m in MLIST)
            eps_rep[i, r] = EPS16 + 2 * dg
    eh = eps_rep.mean(0); es = eps_rep.std(0, ddof=1) / np.sqrt(R_REP)
    print("## TELESCOPED CUMULATIVE eps_hat_r = eps16 + 2 sum_m 4^-m (gamma_r(m)-gamma_16(m))  [STATISTICAL, not the exact ladder]")
    print(f"   anchor eps_16 = {EPS16:.4e} (exact/nu-validated); everything below is a splitting ESTIMATE")
    print(f"   {'r':>3} {'eps_hat':>12} {'SE':>11} {'(hat-eps16)/SE':>15}")
    for r in CURVE_R:
        z = (eh[r] - EPS16) / es[r] if es[r] > 0 else float('nan')
        star = '  <-- curve point'
        print(f"   {r:>3} {eh[r]:>12.4e} {es[r]:>11.3e} {z:>15.2f}{star}")
    print()

    # ---- revised pre-committed decision rule (peak-based, curve not endpoint) ----
    peak_r = int(16 + np.argmax(eh[16:R_DEEP + 1]))
    drop = eps_rep[:, peak_r] - eps_rep[:, R_DEEP]           # per-replica (correlated) drop peak->40
    drop_m = drop.mean(); drop_se = drop.std(ddof=1) / np.sqrt(R_REP)
    rise = eps_rep[:, R_DEEP] - EPS16                        # per-replica rise anchor->40
    rise_m = rise.mean(); rise_se = rise.std(ddof=1) / np.sqrt(R_REP)
    print(f"## DECISION (pre-committed, revised for the curve)")
    print(f"   peak of eps_hat at r={peak_r} (eps={eh[peak_r]:.4e}); drop peak->40 = {drop_m:.3e} +/- {drop_se:.3e} "
          f"({drop_m/drop_se if drop_se>0 else float('nan'):+.1f}sigma)")
    print(f"   rise anchor->40 = {rise_m:.3e} +/- {rise_se:.3e} ({rise_m/rise_se if rise_se>0 else float('nan'):+.1f}sigma)")
    if drop_m > 3 * drop_se and peak_r < R_DEEP:
        if peak_r <= 33:
            print(f"   => TURNOVER CONFIRMED, peak r={peak_r} in/near [27,31] => LOG3 STRUCTURE + 7/15 BOTH SUPPORTED.")
        else:
            print(f"   => TURNOVER CONFIRMED, peak r={peak_r} PAST ~33 => 7/15 SUPPORTED but LOG3 (x3 crossing) REFUTED.")
    elif rise_m > 3 * rise_se and peak_r >= R_DEEP - 1:
        print(f"   => NO TURNOVER through r=40 ({rise_m/rise_se:.1f}sigma rising) => 7/15 IN REAL TROUBLE.")
    else:
        print(f"   => INCONCLUSIVE (drop {drop_m/drop_se if drop_se>0 else 0:.1f}s, rise {rise_m/rise_se if rise_se>0 else 0:.1f}s; neither >=3sigma).")
    print("   NB: whatever this says, it is a SPLITTING ESTIMATE with an error budget -- NOT a proof; the exact ladder still stops at r=16.")
    print(f"\n# total extensions ~ {R_REP*len(MLIST)*R_DEEP*N:,}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
