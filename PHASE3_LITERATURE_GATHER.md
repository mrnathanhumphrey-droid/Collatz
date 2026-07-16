# Phase 3 — Literature Gather (for the spectral-gap brief)

**Purpose.** Map the on-disk literature onto the proof brief's lemmas L1–L4, with exact citations, and state honestly what plugs in vs. what is genuinely ours. Mined 2026-07-16 from `Bourgain-Konyagin/` and `references/Q-sweep/` (three parallel deep reads). Companion to `PHASE3_SPECTRAL_GAP_BRIEF.md`.

## Summary table

| Lemma | Best source(s) | Status |
|---|---|---|
| **L1** well-posedness (finite matrix per level) | trivial in ultrametric (clopen cylinders) | plugs in |
| **L2** Perron: λ₁=1/3 simple, dominant | Ruelle, *Method of Transfer Operators*, Notices AMS 49 (2002), p.891–892 (Perron–Frobenius); Prop. 1 p.889 | **plugs in directly** — need only primitivity of the level matrix |
| **L4** geometric decay of the correlation sum, uniform in k | Solomyak, *Notes on Bernoulli convolutions*, Def. 4.4 + Thm. 4.5; **Siegel diss. eq. 2.180 (Parseval)** | **plugs in for the DIAGONAL** (`(Σp_v²)^k = 3^{-k}`, exactly geometric, k-uniform automatic). Cross term is L3. |
| **L3** the gap: `d = ord_q(2) ≥ 3 ⟹` strict off-diagonal contraction, **uniform in k** | **NONE.** Boundary side (`d=2`) = Konyagin small-subgroup non-cancellation. Positive side = **novel, ours.** | **the real work** |
| our recursion / functional equation | **Siegel diss. Prop. 2.18, eq. 2.173–2.174** = our `μ̂(ξ)=Σ_v 2^{-v}e(ξ2^{-v})μ̂(q2^{-v}ξ)` | **already proven (cite it)** |

## L2 — Perron (plugs in)

Our level-k transfer operator is a **finite nonnegative matrix** (clopen ultrametric cylinders ⟹ locally-constant functions ⟹ finite-dimensional; essential spectral radius = 0, so "quasi-compactness" is automatic and needs no Doeblin–Fortet/Lasota–Yorke inequality). Ruelle's Perron–Frobenius statement (*"L positivity-preserving + mixing ⟹ e^{P(A)} is a simple eigenvalue with no other eigenvalue of equal modulus"*, Notices p.891–892) then gives `λ₁ = Σ_v p_v² = 1/3` simple and strictly dominant, **provided the level matrix is primitive** — the one thing to verify (the qx+1 mixing/overlap structure should give it). `C_q ≥ 1` is forced by Cauchy–Schwarz (already in R8).

## L4 — geometric decay (diagonal plugs in; cross is L3)

**Solomyak, Def. 4.4 + Thm. 4.5.** The cylinder correlation sum `μ_k := (μ×μ){(ω,τ): |ω∧τ| = k}` for a **product measure** is *exactly* `(Σ_i p_i²)^k` — geometric, ratio `Σ p_v² = λ₁ = 1/3`, **uniform in k automatically** (the ratio is constant). This is the DIAGONAL term, verbatim in our geometry (ultrametric `|ω∧τ|` = q-adic level).

⚠️ **Correction to a tempting oversimplification:** our measure is **not** strongly separated. The IFS images `2^{-v}(1+qZ_q)` coincide when `v ≡ v' (mod d)`, so there **are** overlaps — precisely the collisions `2^{-S} ≡ 2^{-S'} mod q^k` — and they carry the entire `cross(k)` term. Solomyak's clean geometric decay covers only the diagonal; the overlap/cross decay at rate `r_q` is exactly L3 and is **not** handed to us by the Bernoulli-convolution machinery. (Supporting, weaker: Peres–Schlag–Solomyak Prop. 4.1(i), subadditive correlation `a_{m+n} ≤ C a_m a_{n-k}` — but its proof uses Plancherel/overlap geometry we'd strip out.)

**The Pisot/Salem phase boundary of Bernoulli convolutions is VACUOUS here** (it is a single-ratio-on-ℤ arithmetic resonance; ours is q-adic with disjoint-or-equal balls). So our `q=3` boundary is a *different*, genuinely arithmetic phenomenon — not the Pisot obstruction. Useful conclusion: the diagonal gap `λ₁=1/3` is structural, unconditional.

## L3 — the gap (THE WORK; no literature route for the positive side)

**Boundary side — why the gap CLOSES at q=3 — is anchored, two independent ways:**
1. **Konyagin, small-subgroup non-cancellation.** Thm 1.8 (p.21): `|H| ≪ log q ⟹ max_a |S(a,H)| ≥ v|H|` (no cancellation). Worked `|G|=2` example (p.6): `S(1,{1,−1}) = 2cos(2π/q) = |H| + O(q^{-2})`. At `q=3` (`d=2`, `H={1,−1}`) the leading mode is preserved to `O(q^{-2})` — **exactly our `λ₂ = λ₁`, gap closed.**
2. **Siegel diss. eq. 4.191 (p.289–290).** The closed form of the *non-archimedean* transform `χ̂_q` carries a factor `1/(q−3)`, **singular at q∈{1,3}**; Siegel flags (heuristically) that this singularity is the finitely-many-orbit-classes / Collatz case. Independent sighting that q=3 is the critical value.

**Positive side — the gap OPENS for `d ≥ 3` (including small d=3,4) — has NO literature route:**
- **Konyagin exponential-sum bounds are VACUOUS for small d.** Every strict `|S(a,H)| < |H|` needs `|H| > √q` (Thm 1.7), `|H| > q^δ` (Thm 3.3, Bourgain–Konyagin), or additive-energy `|H| ≥ q^c` (Thm 2.1 Garcia–Voloch `N_2(b) ≤ 4|H|^{2/3}`, Thm 2.2/2.8). All need `|H|` a positive power of q; at `d = 3, 4` they say nothing. Konyagin Thm 1.8 says the *opposite* for small subgroups.
- **The natural qualitative argument is FALSE.** "H proper ⟹ not additively closed ⟹ contraction" **fails at q=5**: `2` is a primitive root mod 5, so `H = F_5^*` is the *whole* group, yet `r_5 ≈ 0.62 < 1` has a gap. So the gap is **not** a subgroup-additive-closure property. The distinguisher is `d = 2` vs `d ≥ 3` *specifically* (q=3 and q=5 both have full `H`; only q=3 loses the gap).
- **⇒ L3's positive side must come from DIRECT spectral analysis of the concrete cascade operator** (`probe_25`'s `M`, gate-validated), not a generic subgroup theorem. This is the genuine mathematical content and it is **novel** — no source in the corpus supplies it, and Siegel explicitly flags exactly this (`‖π_k‖²` decay / spectral bound) as **open** (diss. p.92–93: "it remains to be seen whether recursive formulae can be derived … for a general p-Hydra map"; the `3^{-k}` decay he attributes to Tao, q=3 only).

## What Siegel already owns (cite, don't re-derive)

- **Prop. 2.18 (eq. 2.173–2.174):** the functional equation `φ_H(t) = (1/p)Σ_j e^{-2πi{b_j t/d_j}} φ_H(a_j t/d_j)` — **identical to our μ̂ recursion**, two-branch form. For `H=T₃` it is Tao's `φ₃`.
- **eq. 2.180 (Parseval):** `q^{-n} Σ_k |φ_H(k/q^n)|² = Σ_k P(χ_H ≡ k mod q^n)²` — the right side **is our `‖π_k‖²`**. Our L²-object = his averaged `|φ_H|²`.
- He does **not** prove the decay, has **no** transfer operator, and does **not** have our sum-to-zero conservation identity `Σ_j M(η₀+j q^k)=0`. Those are ours.

## Net assessment

- **The RPF framework (L1, L2, L4-diagonal) is standard and citable** — Ruelle Perron–Frobenius + Solomyak geometric decay + Siegel's exact recursion/Parseval. In the ultrametric exact-contraction setting the usual hard hypotheses (Hölder potentials, bounded distortion, transversality, Pisot non-resonance) are trivial or vacuous.
- **The crux (L3 positive side: `d ≥ 3 ⟹` uniform off-diagonal gap) is genuinely novel and is the whole mathematical content of Result 1.** Literature supplies the *boundary* (`d=2 ⟹` no gap, Konyagin + Siegel) but not the *interior*. It must be proved by direct analysis of the concrete operator `M`, and the value `r_q` (algebraic, no closed form, R28) is a by-product of that spectrum.
- **Positioning:** Siegel independently reached our recursion and Parseval identity and *flagged our remaining step as the open problem*. We are doing the flagged-open part, with an operator he didn't build — not duplicating.

## Files
On disk: `Bourgain-Konyagin/Konyagin_Lectures.pdf` (load-bearing for L3 boundary), `references/Q-sweep/Ruelle_dynamical_zeta_transfer_operators.pdf` (L2), `references/Q-sweep/Solomyak_Bernoulli_notes.pdf` + `PeresSchlagSolomyak_sixty_years_bernoulli.pdf` (L4), `references/Q-sweep/Siegel2024_pq_adic_Collatz_consolidated.pdf` (recursion Prop. 2.18, Parseval eq. 2.180, q=3 pole eq. 4.191). Cached text: scratchpad `{C2024,FP2,C2019}.txt`.
