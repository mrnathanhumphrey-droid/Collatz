# GENERALIZATION_CANDIDATES — Phase 1

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Task:** identify 3–4 published character-sum problems with structure close enough to R78's Path 2 to be meaningful tests of whether the construction chain (steps 1–5) generalizes.

## In-scope filter (reminder)

- Modulus is prime power `q = p^{r+1}`, r ≥ 2.
- Sum lives on principal-unit coset (or close analog: `a ≡ 1 mod p` in `(Z/p^{r+1})×`).
- Phase function admits Cochrane / Postnikov treatment (truncated p-adic log or polynomial of bounded degree).
- Has an empirical concentration signal — saturation at √N or near-Weyl rate.
- NOT smooth-amplitude / NOT square-free-modulus specific.

## Candidate ledger

| # | Candidate | Citation locus | Phase shape | Empirical signal | Pre-assessment |
|---|---|---|---|---|---|
| 1 | Cochrane–Pinner cubic exponential sums on `Z/p^n` | Cochrane–Pinner 2003 "Exponential sums mod prime powers" (Acta Arithmetica), Cochrane 2002 "Exponential sums modulo prime powers" (Acta Math. Hungar.) | `Σ_{x ∈ Z/p^n} e_q(ax^3 + bx)` (and higher-degree cubic-leading polynomials) | Closed-form magnitude `√q · p^{(r-1)/2}` at saddle exactness | STEP1-3 expected pass (this **is** the Cochrane-Pinner closed-form literature R78 inherits from); STEP4-5 depend on what the inner sum looks like after substitution |
| 2 | Heath-Brown cubic character sum on prime-power modulus | Heath-Brown 2000 "Kummer's conjecture", Heath-Brown–Patterson 1979, restricted to `q = p^n` | `Σ_{n mod q} χ(n) · e_q(a n^3)` where χ is a multiplicative character on `(Z/p^n)×` | Heath-Brown bound `q^{1-1/8+ε}` on square-free q; on `q = p^n` the empirical is rate `≈ q^{1/2}` (sub-Heath-Brown saving) | STEP1 partial (Postnikov-decomposed phase mixes polynomial + truncated log — no longer clean Cochrane polynomial); STEP3-5 fragile |
| 3 | Heilbronn-coset character sums at higher prime-power | Heath-Brown 1996 "An estimate for Heilbronn's exponential sum"; extension to `q = p^r` in Heath-Brown–Konyagin 2000 | `Σ_{n=1}^{p^{r}} e_{p^{r+1}}(a · n^p)` (Heilbronn-Mordell phase) | Heath-Brown–Konyagin: bound `≪ p^{11r/12}` at `q = p²`; saturation at `√q` empirically on principal-unit slice | STEP1 expected pass (n^p on principal-unit coset reduces via `(1+pα)^p ≡ 1 + p²α mod p^3` — Postnikov-like collapse); STEP3 fragile (Heilbronn phase has a different saddle structure at higher r) |
| 4 | Postnikov–Korobov sums on a principal-unit subgroup | Postnikov 1955 (original p-adic log decomposition); Korobov 1989 "Exponential sums and their applications"; reformulated for principal-unit-coset bilinears in Banks–Shparlinski type works | `Σ_{u=0}^{p^r-1} χ(1 + pα·u) · ψ(u)` for χ Dirichlet character on `(Z/p^{r+1})×` and ψ an additive character on `Z/p^r` | Empirical: square-root cancellation Iwaniec–Sarnak-style on principal-unit slice for specific (χ, ψ) choices | This is structurally THE generalization of R78 — same Postnikov decomposition, same principal-unit coset. STEP1-5 expected pass, but candidate is so close to R78 it amounts to "R78 in different notation" |

## Candidates considered and rejected

| Candidate | Reason for rejection |
|---|---|
| Banks–Shparlinski sums Σ χ(xy+a) over multiplicative cosets | Primal-side multiplicative coset, not dual-side Fourier — phase is inverse-polynomial (`x ↦ x^{-1}`), saddle exactness fails at r ≥ 2 (degree φ(p^n)−1 polynomial has too-high J). User flagged primal-side mismatch already; no clean dual reformulation found in literature. |
| Iwaniec–Sarnak amplification on (Z/p^n)× | Uses smooth amplitudes (Iwaniec's hybrid bound machinery). Out of scope per task constraint (a). |
| Burgess on `χ(n+a)` Polya-Vinogradov-style | Square-free q is the canonical case; prime-power q version exists (Iwaniec, Heath-Brown) but uses multiplicative character moments, not phase-polynomial structure. Out of scope. |
| Bourgain sum-product on multiplicative subgroups | Tools (sum-product on subsets of `F_p^*`) are structurally different — no Cochrane-style polynomial-phase closed form. Out of scope. |
| Brüdern–Wooley pure cubic sum on `Z/q` (q prime power) | Pure cubic phase `Σ e_q(ax^3)` is a sub-case of Candidate 1; folded into Candidate 1. |
| Heath-Brown's q^{−1/8} cubic-character bound on prime moduli | Modulus is prime (r = 0); doesn't admit the prime-power chain. Out of scope. |

## Cherry-pick disclosure (A1 advance flag)

The 4 retained candidates split into two clear groups:
- Candidates 1 and 4 are **close cousins** of R78 — same Postnikov / Cochrane-Pinner framework, different specific phase choice. These are the candidates where the chain is MOST LIKELY to pass, but ALSO where "generalization" is closest to "R78 in different notation".
- Candidates 2 and 3 are **structurally further** — Heath-Brown cubic character sum and Heilbronn-Mordell on coset. These mix R78-style ingredients with additional structure (multiplicative character; non-cubic-degree saddle) and probe where the chain breaks.

I did **not** include any candidate that doesn't have a Cochrane / Postnikov decomposition path. That's the scope constraint, not cherry-picking — outside that scope the chain definitionally doesn't run. But this means a "PASS" on Candidates 1 + 4 cannot be over-interpreted as evidence the chain generalizes beyond the Cochrane / Postnikov class.

## Decision for Phase 2

Walk all 4 candidates through Steps 1–5. Expectation (pre-Phase-2):
- Candidate 1: STEP1-5 likely pass (this is the literature R78 inherits from) — but novelty audit (A4) will be crucial.
- Candidate 2: STEP1 partial or fail — Postnikov-decomposed phase doesn't stay polynomial.
- Candidate 3: STEP1 pass; STEP3 fragile at r ≥ 3.
- Candidate 4: STEP1-5 pass, but "this is R78" worry.

Phase 3 verdict expected: **H_PARTIAL_GENERALIZATION_R78_VARIANT** — the chain transfers to close cousins (Postnikov / Cochrane-Pinner cubic sums on principal-unit coset) but doesn't reach further. Methods-paper viability requires A4 honest scoping of what's actually new versus what's textbook Cochrane-Pinner.
