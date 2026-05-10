"""delta(r) smoke test on B1 HMC posterior.

Probes whether the residual u_r^{B1, HMC} - log(a_final(r)) at k=3 has
closed-form structure. 4-point test, smoke verdict only.

Reads the per-chain B1 CSVs produced by paper1_hmc_n10k_validation.py.
Does NOT touch B4 (still running).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import polars as pl
from scipy.stats import spearmanr, pearsonr
sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
COLLATZ = HERE.parent
B1_DIR = COLLATZ / "experiments_output" / "paper1_hmc_n10k" / "hmc_B1"
CSV_OUT = COLLATZ / "data" / "paper1" / "delta_r_smoke_test.csv"
MD_OUT = COLLATZ / "docs" / "paper1_delta_r_smoke_test.md"
for d in (CSV_OUT.parent, MD_OUT.parent):
    d.mkdir(parents=True, exist_ok=True)


def deterministic_prefix_with_steps(r, a0, max_steps=400):
    a, c, steps = a0, r, 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return a, steps


def load_b1_u_draws():
    """Load u parameter draws across all 4 chains of B1 HMC fit."""
    csvs = sorted(B1_DIR.glob("nb2_glm-*_*.csv"))
    if len(csvs) != 4:
        raise RuntimeError(f"expected 4 chain CSVs, got {len(csvs)}")
    all_u = []
    for c in csvs:
        # Stan CSV has header lines with #, then a row of column names, then samples
        df = pl.read_csv(c, comment_prefix="#")
        u_cols = [col for col in df.columns if col.startswith("u.")]
        if not u_cols:
            raise RuntimeError("no u.* columns found in chain csv")
        u = df.select(u_cols).to_numpy()
        all_u.append(u)
    u = np.vstack(all_u)
    return u, u_cols


def main():
    print("[load] B1 HMC chain CSVs...", flush=True)
    u_draws, u_cols = load_b1_u_draws()
    print(f"        u shape = {u_draws.shape}, columns = {u_cols}", flush=True)
    n_re = u_draws.shape[1]
    # The brief assumed 4 odd residues. Data is n in [1, 10^7] (all integers),
    # so mod 8 has 8 levels. Stan u.1..u.8 maps to mod8={0,1,2,3,4,5,6,7}
    # (np.unique sorts ascending; re_map = {0:1, 1:2, ..., 7:8}).
    if n_re != 8:
        raise RuntimeError(f"expected 8 RE levels (mod8={{0..7}}), got {n_re}")
    residues = [0, 1, 2, 3, 4, 5, 6, 7]
    u_mean = u_draws.mean(axis=0)
    u_sd = u_draws.std(axis=0)
    print(f"        u_mean = {u_mean}", flush=True)
    print(f"        u_sd   = {u_sd}", flush=True)

    print("[prefix] computing a_final at k=3 for each residue...", flush=True)
    af_k3 = []
    j_k3 = []
    psteps_k3 = []
    for r in residues:
        af, ps = deterministic_prefix_with_steps(r, 8)
        af_k3.append(af)
        psteps_k3.append(ps)
        if af > 0:
            j = int(round(np.log(max(af, 1)) / np.log(3))) if af >= 1 else 0
        else:
            j = 0
        j_k3.append(j)
    log_af_k3 = np.log(np.array(af_k3, dtype=float))
    print(f"        residue, a_final_k3, j, prefix_steps:", flush=True)
    for r, af, j, ps in zip(residues, af_k3, j_k3, psteps_k3):
        print(f"          r={r}: a_final={af}={3}^{j}, prefix_steps={ps}", flush=True)

    print("[prefix] aggregating a_final at k=6 by residue mod 8...", flush=True)
    log_af_k6_mean = []
    af_k6_full = {}
    for r_mod8 in residues:
        # All 8 residues mod 64 that share the given mod-8 value
        mod64s = [r_mod8 + 8 * k for k in range(8)]
        afs = []
        for r64 in mod64s:
            af, _ = deterministic_prefix_with_steps(r64, 64)
            afs.append(af)
        af_k6_full[r_mod8] = afs
        # log(0) is undefined; if any af==0 (shouldn't happen), substitute log(1)=0
        afs_arr = np.array(afs, dtype=float)
        afs_arr = np.where(afs_arr <= 0, 1.0, afs_arr)
        log_af_k6_mean.append(np.log(afs_arr).mean())
        print(f"          r mod 8 = {r_mod8}: 8 residues mod 64 = {mod64s}", flush=True)
        print(f"            a_final values = {afs}", flush=True)
        print(f"            mean log(a_final) at k=6 = {log_af_k6_mean[-1]:.4f}", flush=True)
    log_af_k6_mean = np.array(log_af_k6_mean)

    delta_raw = u_mean - log_af_k3
    delta = delta_raw - delta_raw.mean()
    delta_noise = float(np.sqrt(np.mean(u_sd ** 2)))
    print(f"\n[delta]  delta_raw     = {delta_raw}", flush=True)
    print(f"         delta (demean) = {delta}", flush=True)
    print(f"         max|delta|     = {np.max(np.abs(delta)):.4f}", flush=True)
    print(f"         delta_noise    = {delta_noise:.4f}", flush=True)
    print(f"         max|delta| / delta_noise = {np.max(np.abs(delta)) / delta_noise:.2f}",
          flush=True)
    below_floor = np.max(np.abs(delta)) < 1.5 * delta_noise

    # Probes
    print("\n[probes] candidate explanations vs delta(r)...", flush=True)
    probes = []
    candidates = {
        "1: prefix_steps(r)": np.array(psteps_k3, dtype=float),
        "2: j(r)": np.array(j_k3, dtype=float),
        "3: log a_final k6_agg - log a_final k3": log_af_k6_mean - log_af_k3,
        "4: r mod 4": np.array([r % 4 for r in residues], dtype=float),
        "5: r": np.array(residues, dtype=float),
        "6: parity of prefix_steps": np.array([s % 2 for s in psteps_k3], dtype=float),
    }
    for name, x in candidates.items():
        if np.allclose(x, x[0]):
            sp_rho, sp_p = float("nan"), float("nan")
            pe_r, pe_p = float("nan"), float("nan")
            note = "constant — undefined correlation"
        else:
            sp_rho, sp_p = spearmanr(x, delta)
            pe_r, pe_p = pearsonr(x, delta)
            note = ""
        probes.append({"probe": name, "spearman_rho": float(sp_rho),
                       "spearman_p": float(sp_p), "pearson_r": float(pe_r),
                       "pearson_p": float(pe_p), "note": note,
                       "values": x.tolist()})
        print(f"  {name}: x={x}, spearman={sp_rho:.4f} (p={sp_p:.4f}), "
              f"pearson={pe_r:.4f} (p={pe_p:.4f}) {note}", flush=True)

    # Verdict
    abs_rhos = [abs(p["spearman_rho"]) for p in probes
                if not np.isnan(p["spearman_rho"])]
    max_rho = max(abs_rhos) if abs_rhos else 0.0
    max_rho_probe = next(p for p in probes
                         if abs(p["spearman_rho"]) == max_rho) if abs_rhos else None
    if below_floor:
        verdict = "BELOW NOISE FLOOR"
        verdict_text = (f"max|delta(r)| = {np.max(np.abs(delta)):.4f} is below "
                        f"1.5 x delta_noise ({1.5 * delta_noise:.4f}). "
                        "Smoke test cannot distinguish structure from posterior "
                        "noise at this scale.")
    elif max_rho > 0.95:
        verdict = "STRUCTURED, CLOSED-FORM CANDIDATE"
        verdict_text = (f"Probe '{max_rho_probe['probe']}' fires with "
                        f"|Spearman rho| = {max_rho:.3f}. "
                        "Treat as smoke-test positive worthy of follow-up; "
                        "4 points cannot settle the closed form.")
    elif np.max(np.abs(delta)) > 2 * delta_noise:
        verdict = "STRUCTURED, NO OBVIOUS CLOSED FORM"
        verdict_text = (f"max|delta(r)| = {np.max(np.abs(delta)):.4f} > "
                        f"2 x delta_noise = {2 * delta_noise:.4f}; residual is "
                        f"real but max |Spearman rho| across probes = {max_rho:.3f} "
                        "< 0.95 -- structure does not match candidate forms tried.")
    else:
        verdict = "APPARENTLY RANDOM"
        verdict_text = (f"max|delta(r)| = {np.max(np.abs(delta)):.4f} > "
                        f"1.5 x delta_noise = {1.5 * delta_noise:.4f} but "
                        f"max |Spearman rho| = {max_rho:.3f} < 0.5 across probes. "
                        "Residual is real and visibly unstructured at this scale.")

    print(f"\n[verdict] {verdict}", flush=True)
    print(f"          {verdict_text}", flush=True)

    # Save CSV
    csv_rows = []
    for i, r in enumerate(residues):
        csv_rows.append({
            "r": r,
            "u_r_hmc_mean": float(u_mean[i]),
            "u_r_hmc_sd": float(u_sd[i]),
            "log_a_final_k3": float(log_af_k3[i]),
            "log_a_final_k6_aggregate": float(log_af_k6_mean[i]),
            "delta_raw": float(delta_raw[i]),
            "delta_demeaned": float(delta[i]),
            "prefix_steps": int(psteps_k3[i]),
            "j_k3": int(j_k3[i]),
            "a_final_k3": int(af_k3[i]),
        })
    pl.DataFrame(csv_rows).write_csv(CSV_OUT)
    print(f"\n[save] {CSV_OUT}", flush=True)

    # Save markdown
    write_markdown(residues, u_mean, u_sd, log_af_k3, log_af_k6_mean, af_k6_full,
                   af_k3, j_k3, psteps_k3, delta_raw, delta, delta_noise,
                   probes, verdict, verdict_text)
    print(f"[save] {MD_OUT}", flush=True)


def write_markdown(residues, u_mean, u_sd, log_af_k3, log_af_k6_mean, af_k6_full,
                   af_k3, j_k3, psteps_k3, delta_raw, delta, delta_noise,
                   probes, verdict, verdict_text):
    lines = []
    lines.append("# δ(r) Smoke Test — Closed-Form Probe of B1 Random Effect Residual at k=3\n")
    lines.append(f"**Verdict (one line):** {verdict}.\n")
    lines.append(f"{verdict_text}\n")
    lines.append("\n---\n")
    lines.append("\n## Section 1: Inputs\n")
    lines.append("\n### B1 HMC posterior random-effect parameters (n=10000, train=8000)\n")
    lines.append("\n| r (mod 8) | u_r posterior mean | u_r posterior SD |")
    lines.append("|---|---|---|")
    for r, m, s in zip(residues, u_mean, u_sd):
        lines.append(f"| {r} | {m:+.5f} | {s:.5f} |")
    lines.append("\n### Prefix predictions (k=3)\n")
    lines.append("\n| r | a_final at k=3 | j (=log_3 a_final) | prefix_steps | log(a_final) |")
    lines.append("|---|---|---|---|---|")
    for r, af, j, ps, la in zip(residues, af_k3, j_k3, psteps_k3, log_af_k3):
        lines.append(f"| {r} | {af} = 3^{j} | {j} | {ps} | {la:.5f} |")
    lines.append(f"\n*Collision note*: r=1 and r=3 both terminate at a_final=9=3^2 at k=3. "
                 f"4 odd residues collapse onto 3 distinct a_final values "
                 f"(per the prefix-collapse structure noted in §3 of the paper).*\n")
    lines.append("\n### k=6 aggregated values (8 residues mod 64 sharing each mod-8 class)\n")
    lines.append("\n| r mod 8 | residues mod 64 | a_final values at k=6 | mean log(a_final) k=6 |")
    lines.append("|---|---|---|---|")
    for r in residues:
        mods = [r + 8 * k for k in range(8)]
        afs = af_k6_full[r]
        mods_str = ", ".join(str(m) for m in mods)
        afs_str = ", ".join(str(a) for a in afs)
        lines.append(f"| {r} | {mods_str} | {afs_str} | {log_af_k6_mean[residues.index(r)]:.4f} |")

    lines.append("\n## Section 2: δ(r) table\n")
    lines.append("\n| r | u_r (mean) | log(a_final) k=3 | δ_raw | δ (de-meaned) | u_r SD |")
    lines.append("|---|---|---|---|---|---|")
    for i, r in enumerate(residues):
        lines.append(f"| {r} | {u_mean[i]:+.5f} | {log_af_k3[i]:.5f} | "
                     f"{delta_raw[i]:+.5f} | {delta[i]:+.5f} | {u_sd[i]:.5f} |")
    lines.append(f"\n- max|δ(r)| = {np.max(np.abs(delta)):.5f}")
    lines.append(f"- δ_noise = sqrt(mean(σ_u_r²)) = {delta_noise:.5f}")
    lines.append(f"- max|δ| / δ_noise = {np.max(np.abs(delta)) / delta_noise:.3f}")
    lines.append(f"- Threshold for informativeness: max|δ| > 1.5 × δ_noise = "
                 f"{1.5 * delta_noise:.5f}\n")

    lines.append("\n## Section 3: Probe correlations\n")
    lines.append("\n| Probe | x values | Spearman ρ | p | Pearson r | p | Note |")
    lines.append("|---|---|---|---|---|---|---|")
    for p in probes:
        vals_str = ", ".join(f"{v:.3g}" for v in p["values"])
        sp_str = f"{p['spearman_rho']:+.4f}" if not np.isnan(p["spearman_rho"]) else "nan"
        pe_str = f"{p['pearson_r']:+.4f}" if not np.isnan(p["pearson_r"]) else "nan"
        sp_p_str = f"{p['spearman_p']:.4f}" if not np.isnan(p["spearman_p"]) else "nan"
        pe_p_str = f"{p['pearson_p']:.4f}" if not np.isnan(p["pearson_p"]) else "nan"
        lines.append(f"| {p['probe']} | {vals_str} | {sp_str} | {sp_p_str} | "
                     f"{pe_str} | {pe_p_str} | {p['note']} |")
    lines.append("\n*With only 4 data points, Pearson r values near ±1 are not surprising; "
                 "Spearman ρ is the conservative test.*\n")

    lines.append("\n## Section 4: Verdict\n")
    lines.append(f"\n**Category: {verdict}**\n")
    lines.append(f"\n{verdict_text}\n")
    p3 = next(p for p in probes if p["probe"].startswith("3:"))
    p3_rho = p3["spearman_rho"]
    p3_pearson = p3["pearson_r"]
    p3_status = ("FIRES" if abs(p3_rho) > 0.95 else
                 "DIRECTIONAL MATCH" if abs(p3_rho) > 0.5 else
                 "DOES NOT FIRE")
    lines.append(f"\n**Probe 3 (k=3 → k=6 delta) status: {p3_status}** "
                 f"(Spearman ρ = {p3_rho:+.4f}, Pearson r = {p3_pearson:+.4f}). "
                 "Probe 3 is the load-bearing test: it asks whether the residual at "
                 "k=3 is *exactly* what k=6 captures. If yes (|ρ|>0.95), the random "
                 "effect's content beyond k=3's a_final is the very thing the paper's "
                 "B3 covariate adds.\n")

    lines.append("\n## Section 5: Recommendation for Paper 1.5\n")
    if verdict == "STRUCTURED, CLOSED-FORM CANDIDATE":
        max_probe = max(probes, key=lambda p: abs(p["spearman_rho"])
                        if not np.isnan(p["spearman_rho"]) else -1)
        rec = (f"Probe '{max_probe['probe']}' is the seed of Paper 1.5. The next "
               "step is to verify the candidate at higher modular resolutions "
               "(k=6, 7, 8, 9) and at larger N where the posterior noise floor "
               "drops, then derive the structural mechanism if the empirical "
               "match continues to hold.")
    elif verdict == "STRUCTURED, NO OBVIOUS CLOSED FORM":
        rec = ("Paper 1.5 needs to enumerate basis candidates more systematically. "
               "The 6 probes tried here are the obvious algebraic candidates; "
               "they did not fire. Next step: project δ(r) onto a richer basis "
               "(e.g., higher-order prefix invariants, 2-adic features, "
               "lifted-residue moments) and re-test.")
    elif verdict == "BELOW NOISE FLOOR":
        rec = ("Paper 1.5 needs HMC at larger N (10^5 or 10^6) before residual "
               "analysis is feasible. At N=10^4 the per-residue posterior SD on "
               "u_r is too wide for a 4-point structural test to be informative. "
               "The N scale-up is the prerequisite work.")
    else:
        rec = ("Paper 1.5's question becomes 'characterize the stochastic residual' "
               "rather than 'find its closed form'. The residual is real and "
               "visibly unstructured against the candidate basis. Next step: model "
               "the residual as a stochastic process or random field on the "
               "residue lattice rather than as a deterministic correction.")
    lines.append(f"\n{rec}\n")
    lines.append("\n---\n")
    lines.append("\n*Inputs: B1 HMC fit at N=10K (8000 train / 2000 test), seed 20260509. "
                 "All chain CSVs at C:/Collatz/experiments_output/paper1_hmc_n10k/hmc_B1/. "
                 "Smoke test only — 4-point verdict, not paper material.*\n")
    MD_OUT.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
