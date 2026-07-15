"""
result_1b_halving_pgf.py -- Probe 1B: two-root-of-unity audit of the halving PGF.
Pre-reg: PRE_REG_1B_GEOM_PGF_2026_07_13.md.

Audits whether the EMPIRICAL v-measure (not idealized Geom(1/2)) keeps the two
constants that the sec-CONJ mechanism 7/45 = N(2-w)/(3^2*(1+4)) rests on:
    G_hat(-1) = E[(-1)^v]        target -1/3  (the "9": class mass -> 1:4)
    |G_hat(w)|^2 = E[w^v] norm^2 target  1/7  (the "7": N(2-w)=7),  w=e^{2pi i/3}
Geom(1/2) baseline: G(-1)=-1/3 exactly, |G(w)|^2=1/7 exactly.

Data: data/v_seq_N8388608.parquet -- long format (n, step_idx, v), 427,241,688
per-step 2-adic valuations v=v2(3n+1) over 2.79M odd coprime-to-3 starts.
Marginal used: UNCONDITIONAL per-step v-marginal. Bracket-stratified by step_idx
(descent-funnel position) per pre-reg step 3 (mandatory; Bohr aggregate-vs-funnel).

Decision (pre-reg 4.5): H_ROBUST iff |dG(-1)|<=0.01 AND |d|G(w)|^2|<=0.01 in EVERY
bracket; pooled-pass-but-bracket-fail => H_BRACKET_SPLIT; pooled-fail => H_DRIFT.
Exact w = (-1 + i*sqrt(3))/2. No aggregation before per-bracket report.
Not at stake: the rigorous 7/45 (R75/R76 exact-rational, separate).
"""
import sys, math, cmath
import numpy as np
import pyarrow.parquet as pq
sys.stdout.reconfigure(encoding="utf-8")
REPO = r"C:\Collatz"
DATA = REPO + r"\data\v_seq_N8388608.parquet"

LOG = []
def log(m=""):
    print(m); LOG.append(str(m))

# step_idx brackets (descent-funnel position). Half-open [lo, hi].
BRACKETS = [("step0", 0, 1), ("step1-2", 1, 3), ("step3-5", 3, 6),
            ("step6-9", 6, 10), ("step10+", 10, 1 << 30)]
VMAX = 64  # v is v2(3n+1) >= 1; cap histogram width (>VMAX folded, should be ~0)

W = complex(-0.5, math.sqrt(3) / 2.0)   # exact e^{2pi i/3}

def ghat_from_counts(counts):
    """counts[v] for v=0..VMAX-1 (v>=1 in data). Return N, mean, Ghat(-1), Ghat(w), |Ghat(w)|^2."""
    N = int(counts.sum())
    vs = np.arange(VMAX)
    mean = float((counts * vs).sum()) / N
    g_m1 = float((counts * ((-1.0) ** vs)).sum()) / N
    wv = np.array([W ** int(v) for v in vs], dtype=complex)
    g_w = complex((counts * wv).sum() / N)
    return N, mean, g_m1, g_w, abs(g_w) ** 2

def main():
    log("# PROBE 1B -- halving-PGF two-root-of-unity audit (empirical v-measure)")
    log(f"# data: {DATA}")
    log(f"# exact w = e^(2pi i/3) = {W:.15f}")
    log(f"# Geom(1/2) baseline: G(-1) = -1/3 = {-1/3:.10f},  |G(w)|^2 = 1/7 = {1/7:.10f}")
    log("")

    nb = len(BRACKETS)
    pooled = np.zeros(VMAX, dtype=np.int64)
    bcounts = [np.zeros(VMAX, dtype=np.int64) for _ in range(nb)]
    lo = np.array([b[1] for b in BRACKETS]); hi = np.array([b[2] for b in BRACKETS])

    pf = pq.ParquetFile(DATA)
    total = 0; overflow = 0
    for batch in pf.iter_batches(batch_size=4_000_000, columns=["step_idx", "v"]):
        si = batch.column("step_idx").to_numpy(zero_copy_only=False).astype(np.int64)
        v = batch.column("v").to_numpy(zero_copy_only=False).astype(np.int64)
        ov = v >= VMAX
        if ov.any():
            overflow += int(ov.sum()); v = np.minimum(v, VMAX - 1)
        pooled += np.bincount(v, minlength=VMAX)[:VMAX]
        for bi in range(nb):
            mask = (si >= lo[bi]) & (si < hi[bi])
            if mask.any():
                bcounts[bi] += np.bincount(v[mask], minlength=VMAX)[:VMAX]
        total += len(v)
    log(f"streamed rows: {total:,}  (overflow v>={VMAX}: {overflow})")
    log("")

    # smoke check: pooled mean ~ 2.102
    Np, meanp, gm1p, gwp, gw2p = ghat_from_counts(pooled)
    log("## smoke: pooled per-step v-marginal")
    log(f"   N={Np:,}  mean_v={meanp:.5f}  (expect ~2.102 per result_density_one_v2_bounds.md)")
    if abs(meanp - 2.102) > 0.05:
        log(f"   !! mean_v deviates from 2.102 by {meanp-2.102:+.4f} -- reconcile before trusting")
    log("")

    csv = ["bracket,N,mean_v,Ghat_m1,dG_m1,Re_Ghat_w,Im_Ghat_w,absGhat_w2,d_absG_w2"]
    def row(name, counts):
        N, mean, gm1, gw, gw2 = ghat_from_counts(counts)
        d_m1 = gm1 - (-1/3); d_w2 = gw2 - 1/7
        csv.append(f"{name},{N},{mean:.6f},{gm1:.8f},{d_m1:+.8f},{gw.real:.8f},"
                   f"{gw.imag:.8f},{gw2:.8f},{d_w2:+.8f}")
        return N, mean, gm1, d_m1, gw, gw2, d_w2

    log("## per-bracket + pooled  (targets: Ghat(-1)=-0.33333, |Ghat(w)|^2=0.142857)")
    log(f"   {'bracket':10s} {'N':>13s} {'mean_v':>7s} {'Ghat(-1)':>10s} {'dG(-1)':>9s} "
        f"{'|Ghat(w)|^2':>11s} {'d|G(w)|^2':>10s}  verdict")
    results = {}
    for (name, _, _), counts in list(zip(BRACKETS, bcounts)) + [(("POOLED", 0, 0), pooled)]:
        N, mean, gm1, d_m1, gw, gw2, d_w2 = row(name, counts)
        ok = (abs(d_m1) <= 0.01) and (abs(d_w2) <= 0.01)
        results[name] = (d_m1, d_w2, ok)
        log(f"   {name:10s} {N:>13,d} {mean:>7.4f} {gm1:>10.5f} {d_m1:>+9.5f} "
            f"{gw2:>11.6f} {d_w2:>+10.5f}  {'OK' if ok else 'DRIFT'}")
    log("")

    # decision
    pooled_ok = results["POOLED"][2]
    brackets_ok = all(results[b[0]][2] for b in BRACKETS)
    if pooled_ok and brackets_ok:
        verdict = "H_ROBUST"
    elif pooled_ok and not brackets_ok:
        verdict = "H_BRACKET_SPLIT"
    else:
        verdict = "H_DRIFT"
    log(f"## VERDICT: {verdict}")
    log(f"   pooled within 0.01: {pooled_ok};  every bracket within 0.01: {brackets_ok}")
    if verdict == "H_ROBUST":
        log("   => both constants survive the true v-measure. The N(2-w)/(3^2*5) mechanism")
        log("      is measure-ROBUST, not a Geom(1/2) idealization. Phase-2 target stands:")
        log("      derive R75/R76's 7 as N(2-w).")
    elif verdict == "H_BRACKET_SPLIT":
        log("   => holds pooled but splits by descent-funnel position. The mechanism is")
        log("      position-dependent; do NOT read the pooled pass as a clean robustness result.")
    else:
        log("   => the mechanism is Geom(1/2)-IDEALIZED; 7=N(2-w) is approximate, not the exact")
        log("      origin of the rigorous 7/45. Motivates a THEOREM_C_745 Geom(1/2)-substitution audit.")
    log("")

    with open(REPO + r"\result_1b_pgf_data.csv", "w", encoding="utf-8") as f:
        f.write("\n".join(csv) + "\n")
    with open(REPO + r"\result_1b_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("[wrote] result_1b_pgf_data.csv, result_1b_log.txt")
    return verdict, results, (Np, meanp)

if __name__ == "__main__":
    main()
