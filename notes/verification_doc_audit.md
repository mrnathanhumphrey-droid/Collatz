# Audit of `milicevic_banks_verification.md`

**Date:** 2026-05-04
**Auditor:** independent review (skeptic mode)
**Target:** `C:\Collatz\milicevic_banks_verification.md`

---

## Per-claim verdicts

### CLAIM 1 — F̂ vs 1̂ notation (§C2). **BUG (cosmetic / framing), but bottom line survives.**

Empirical at r = 3 (script `C:\Collatz\audit_claim1_claim4.py`):

| a | \|F̂(3a)\| | \|1̂(3a)\| |
|---|---|---|
| 1 | 27.0000 | 7.4598 |
| 4 | 27.0000 | 1.9297 |
| 7 | 27.0000 | 1.1906 |
| 10 | 27.0000 | 0.9432 |
| 13 | 27.0000 | 0.8675 |
| 16 | 27.0000 | 0.9040 |
| 19 | 27.0000 | 1.0797 |
| 22 | 27.0000 | 1.5760 |
| 25 | 27.0000 | 3.7553 |

3√q = 27.0000 exactly. **|F̂(3a)| = 3√q is saturated; |1̂(3a)| varies 0.87–7.46.** They are different objects.

§C2 of the doc (lines 304–309) writes "the primal sum |Σ_{n ∈ window} χ(n) e(c·n)|, which equals our 1̂(3·) on the support" then "we already know |1̂(3a)| = 3√q exactly (Theorem 78.3)." That second clause **conflates 1̂ and F̂.** Theorem 78.3 (per `result_78_extended.md` line 95) is about **F̂**, not 1̂.

**Impact:** The "Resolution" paragraph (line 307) recovers — it correctly distinguishes worst-case primal vs Plancherel-saturated dual. So the *argument* lands on the right structural mismatch, but the bridge sentence "we already know |1̂(3a)| = 3√q" is plain wrong. The §C2 conclusion ("Milićević bounds worst-case primal, we have exact value at specific dual frequency") still holds because Milićević applied to the primal would bound `Σ_n χ(n) e(c n)` ≪ q^{1/2−η}; that primal sum, by Theorem 78.3 / the smooth-completion identity, *does* equal something of size 3√q at our frequency (since |F̂(3a)| = 3√q after the Plancherel split). The notation just got swapped one layer too early. Bottom-line outcome 3 is independently established by §A2.2 (Theorem 3 condition (iv) failure) anyway.

### CLAIM 2 — Theorem 3 condition (iv) failure (§A2.2). **CONFIRMED.**

- Definition of ρ_p(y) at `milicevic.txt:599`: "[y+α] ∈ I_{1+α}[λ − ρ_p(y)](Z_p), where ρ_p(y) equals p if ord_p y = 0 and 0 otherwise." Verbatim.
- ord_3(1) = 0 ⇒ ρ_3(1) = 3. Confirmed.
- Theorem 3 statement at `milicevic.txt:2094–2100`: requires λ̃ = min(κ − ρ_p(y), λ) > 0. Verbatim (line 2100).
- κ = 1, ρ_3(1) = 3, λ = ∞ ⇒ λ̃ = min(−2, ∞) = −2 < 0. **Condition (iv) fails.**
- Sanity flag: `milicevic.txt:2085–2090` explicitly notes conditions (i)–(ii) "can occasionally be somewhat relaxed... relevant for p ∈ {2,3} only." That's a relaxation note for (i)–(ii), NOT for (iv). The doc's claim that (iv) is the *hard* failure is correct.
- Whether y = 1 (pure log) is "covered separately" — Milićević's Lemma 13 (line 862-ish, cited by doc at 30) gives the F-class membership for primitive characters with κ_1 = 1 + ι'(2). That puts the Dirichlet-character case into F-class with y = 1 directly; there's no y = 1 carve-out. The doc's framework-application is the right reading.

### CLAIM 3 — Exponent pair arithmetic (§A3.1). **CONFIRMED (with one self-correction the doc already flagged).**

Recomputed via Fraction arithmetic:
- B(0,1) = (1/2, 1/2) ✓
- AB(0,1) = (1/6, 2/3) ✓
- A²B(0,1) = (1/14, 11/14)
- A³B(0,1) = (1/30, 13/15)
- BA³B(0,1) = (11/30, 8/15)
- ABA³B(0,1) = **(11/82, 57/82)** ✓ matches `milicevic.txt:184`
- θ = (11/82 + 57/82)/2 − 1/4 = 68/164 − 41/164 = **27/164 ≈ 0.16463** ✓ matches `milicevic.txt:3571, 3594`

Doc correctly self-flags an early-typo (line 128–129) and lands at θ = 27/164.

### CLAIM 4 — Saddle-class linear phases at r = 3 (§C3 caveat 9). **CONFIRMED.**

Empirical (J_{3,1,4} = 3 ✓):
- L(1) = 0, mod 81 = **0** ✓
- L(4) = 15/2 (as Fraction), mod 81 = **48** = 15·41 mod 81 ✓
- L(7) = 60 (as Fraction; integer because v_3(L(1+9)) ≥ 2), mod 81 = **60** ✓

So:
- P_a(0) = 0 ✓
- P_a(1) = 3 − C_a · 48 mod 81 = **3 − 48·C_a mod 81 = 3 − (15/2)·C_a mod 81** ✓ (the 15/2 is correct as a fraction; reducing mod 81 gives 48; both phrasings are equivalent)
- P_a(2) = 6 − C_a · 60 mod 81 = **6 − 60·C_a mod 81** ✓

Cross-check vs `result_78_extended.md` table at r = 3: all 9 (C_a, s*, P_a(s*)) triples reproduce exactly. **Linear-in-C_a structure holds.** Saddle-class partition path is well-posed at r = 3.

Side note on the doc's "L̃ = L(4)/3 = 5/2" phrasing (line 47 / Brief): As a Fraction, L(4) = 15/2 and L(4)/3 = 5/2. ✓ But the actual code (`path_B_explicit_phase.py:80`) computes `L_tilde = L4_mod // 3 = 48 // 3 = 16`. Both are right — 5/2 mod 27 = 5·14 = 70 mod 27 = 16. Notation slightly schizophrenic across the doc but mathematically consistent.

### CLAIM 5 — B-S G is in summation variable (§B2). **CONFIRMED.**

`banks_shparlinski.txt:175–202`:
- Line 177: S(M,N;G) = Σ_{n=M+1}^{M+N} χ(n) · e(G(n))
- Line 197–199: G ∈ ℝ[x], deg G ≤ C

G is unambiguously a polynomial in the summation variable n. Our P_a(s) is polynomial in the saddle parameter s, and ψ(a) is the Plancherel-dual phase — neither sits in B-S's slot. Doc's "structural mismatch, not parameter issue" framing is correct. There is no reinterpretation that handles dual-side cubic phases inside B-S's S(M,N;G) shape; the dual sum doesn't even have a `χ(n) e(G(n))` form (the kernel 1̂(3a) is not χ at a polynomial argument).

---

## Bottom-line assessment

The doc's headline conclusion **survives audit**. Outcome 3 (neither framework cleanly closes |S_partial| ≪ q^{1/2−δ}) is supported by:

1. **Hard structural failure:** Theorem 3 (iv) λ̃ = −2 < 0 (CLAIM 2 confirmed). Independent of the §C2 notation slip.
2. **Wrong-object failure:** B-S's S(M,N;G) shape doesn't match our dual-side bilinear (CLAIM 5 confirmed).
3. **Wrong-length failure:** Exponent-pair bound at length q^{1/2} gives no sub-trivial saving on |S_partial| (§A3.2; arithmetic checks).

The §C2 notation bug (CLAIM 1) is a *framing slip*, not a load-bearing error. The argument it's supposed to support — that Milićević bounds the worst-case primal and we have a frequency-saturated dual — is correct; the doc just stated it through a wrong intermediate identification (1̂ ≡ F̂). Fixing the slip strengthens the doc, doesn't reverse it.

CLAIMS 3 and 4 are clean. The path-forward §C3-9 (saddle-class partition with linear-in-a phases at r = 3) is mathematically well-posed and the doc's own caveat that this only buys logarithmic savings at r = 3 (not polynomial) is honest.

**No bug found that opens a closed path.**

---

## Recommendation

**(a) Trust the doc and proceed**, with a small revision pass:

1. Fix §C2 lines 304–307: replace "we already know |1̂(3a)| = 3√q exactly (Theorem 78.3)" with "we already know |F̂(3a)| = 3√q exactly (Theorem 78.3); after the Plancherel split, the primal short character sum has size matching this saturation." The argument lands the same; the citation needs to point at F̂, not 1̂.
2. Optional: in Caveat 9 (line 350), surface the explicit values L(4) ≡ 48 mod 81 and L(7) ≡ 60 mod 81 alongside the fraction forms — disambiguates the L̃ = 5/2 vs 16 thing.

The substantive mathematical conclusion (Heath-Brown / Cochrane-style cubic-character-sum bound on (Z/3^r)* is the right next tool, not Milićević or Banks-Shparlinski) stands.
