# Probe D1 — cap-sector recon (Track D's input)

**Date:** 2026-07-18  Dense/direct at L=2,3 (sanctioned); block SpMM subspace iteration at L=4 (Lambda A10,
233k states). Goal: name the effective-model basis and deliver its judges for the three-body reading
(partner + condensing complex pair).

## Bookkeeping rider (Wilson's catch, entered in STATE)
Two labeled sequences, not to be conflated. ρ₂=0.34682666, ρ₃=0.333236, ρ₄=0.33349990132.
- **detuning-vs-c₀** [c₀(L)−ρ]: −2.911e-3 (L2), +9.958e-5 (L3), −1.666e-4 (L4)
- **distance-to-1/3** [1/3−ρ]: −1.349e-2 (L2), +9.73e-5 (L3), −1.666e-4 (L4)  (coincide at L4: c₀(4)≈1/3)

## D1-A — sector content of the top modes (the basis-namer)
Gauge-character mass Σ|f_k|² per sector k (Fourier in log₂ a), dense eigenvectors.

| | partner | leading pair | pair |·|, arg |
|---|---|---|---|
| L=2 (D=6) | 0.3468267 | 0.02232+0.23513j | 0.23619, 1.4762 rad |
| L=3 (D=18) | 0.3332363 | 0.23764+0.18303j | 0.29996, 0.6563 rad |

- **Partner is k=0-LEADING** — but only **40% (L2) → 32% (L3)** of its mass sits in k=0; the rest spreads (L3: k=9,3,15). **The k=0 concentration DECREASES with L.**
- **Pair is k=+1-LEADING** (its conjugate lives at k=−1) — L3 k=1 mass 0.274, then a real tail k=10,16,4. So the pair names a **k=±1** sector, **but it is broad, not pure** (leading sector only ~27%).
- **⚠️ REPORT-LOUDLY caveats for the pen (before deriving):**
  1. The top modes are **sector-LEADING, not sector-concentrated**, and the concentration **weakens with L** (partner k=0: 0.40→0.32). A pure {k=0, ±1} 3×3 captures the leading sectors but a **growing majority of mass lives in other k**.
  2. The subdominant is a **CLUSTER of near-degenerate complex pairs**, not one (L3: |λ|=0.29995 AND 0.29794; L4: 0.32042+0.07524j AND 0.32022+0.07525j). The "pair" is a doublet.
  3. Single-eigenvector sector content is NOT k↔−k symmetric (maxdev 0.15–0.24) — expected (the +imag vector sits at +k; the pair-TOTAL symmetry is trivially forced and not an independent test). The real readout is **which k: k=±1.**
- **Basis named:** the effective model is **{k=0 (partner), k=+1, k=−1}** at leading order — with the above corrections.

## D1-B — coupling tables (raw material, dumped)
`outputs/d1b_NR_tables_q3_L{2,3}.txt` for k∈{0, 1, and the pair's second sector} + k=0. R_k(s) on the D-shell
(R_0 the real autocorrelation R_0(0)=1/3+…, R_1/R_10 complex); N_κ(e′,γ,γ′) twisted counts (k=0 = integer
counts =18; the depth-selection rule is visible, e.g. γ=3→γ′=1 carries N only at k=0). Tables first, per protocol.

## D1-C — the pair at L=4 (the condensation's third point + the phase)
Block subspace iteration (block 6, doubles). Partner **0.333500** (= G4's 0.33349990132 — cross-validated by
an independent method). Leading pair **0.320423 + 0.075243j**, with a near-degenerate twin 0.320223+0.075252j.

| L | partner | \|pair\| | arg(pair) rad | ratio \|pair\|/partner | gap partner−\|pair\| |
|---|---|---|---|---|---|
| 2 | 0.346827 | 0.236186 | 1.47615 | 0.68099 | 0.11064 |
| 3 | 0.333236 | 0.299955 | 0.65631 | 0.90013 | 0.03328 |
| 4 | 0.333500 | 0.329139 | 0.23065 | **0.98692** | **0.00436** |

**THE CONDENSATION IS CONFIRMED.** The complex pair coalesces onto the partner along all three axes: modulus
ratio → 1 (0.68 → 0.90 → 0.987), phase → 0 (1.476 → 0.656 → 0.231, rotating onto the real axis), and gap → 0
(0.111 → 0.033 → 0.0044, closing ~×0.30 then ×0.13). An exceptional point is forming as L→∞ — the partner and
the tower's leading complex pair are merging. (The measured ratios supersede the cited 0.87/0.97; the trend —
condensation — is exactly as posed.)

**Hypothesis under test — does the partner's detuning sign track the pair's phase?** Phase is **monotone
decreasing** (1.476, 0.656, 0.231); detuning sign **oscillates** (−, +, −). A simple monotone phase→sign map
does **NOT** hold across these three points. So: the condensation strongly supports the three-body interference
reading (a pair merging with a real mode is exactly what produces an oscillating detuning), but the **specific
sign law — the braid — is not read off the phase directly; it is the pen's to derive** from the effective
model, and these three phases + three detunings judge it. Reported, not forced.

## D1-D — SKIPPED (reported, per spec)
c₀-side eigenvector / overlap₄ / g₄. At L=4 the partner (0.33350) exceeds c₀ (≈1/3), so the c₀-side right
eigenvector is subdominant in the full operator AND the overlap/g₄ need F2's exact biorthogonal convention;
deferred to avoid a mislabeled EP witness. A–C are the spine.

## Headline for the pen
The three-body reading is **materially supported**: a complex pair (k=±1) is condensing onto the partner
(k=0) — ratio→1, phase→0, gap→0 across L=2,3,4 — and that condensation is the natural source of the partner's
braid. Two things the pen must build in: the modes are **sector-broad** (leading sector ~30%, weakening with L)
and the pair is a **near-degenerate doublet**, so the {k=0,±1} model is a leading reduction with known
corrections. The braid's **sign law is not a direct function of the phase** — it must fall out of the derived
effective model, which the three (phase, detuning) points then judge.

Probes `probes/d1ab_sectors.py`, `probes/d1c_pair.py`; dumps `outputs/d1a_sectors.json`,
`outputs/d1b_NR_tables_q3_L{2,3}.txt`; log `logs/d1c_pair_L4_log.txt`.
