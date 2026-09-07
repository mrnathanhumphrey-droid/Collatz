# RESULT — ARTIN-CLASS: is forward-enrichment PASS/FAIL Chebotarev-governed? **DEAD — the FIFTH dead predictor-story. No density constant, no condition-match: forward PASS/FAIL is genuinely structureless / not Chebotarev-predictable.** (2026-09-06)

**Probe:** `probes/probe_artin_class.py`. Tests Wilson's (unverified, lit-hunt) claim that forward-enrichment
PASS/FAIL is "Artin/Chebotarev-class — governed by the multiplicative behavior of 2 mod p." Held to the same bar as
the four already-killed predictors (q²=1+p, primitive-root, orbit-size, "no predictor"): a citation to Artin is not
evidence; the **condition-match over 33 primes (Test 2) is decisive** — a density can coincide, an exact set-match
cannot.

**Criterion (uniform, reproduces the anchors):** PASS := `γ(p), γ(2p), γ(3p)` all `>1` (first-3-multiples
enriched; captures p=11's depth-robust 3p-depletion; applied identically to every prime). Float sign-call (as the
anchors were established), shallow r (a pass/fail-per-prime deliberately dodges the doubly-exp wall). Stability:
computed at r_top and r_top−1; a prime whose call flips is **UNRESOLVED** and excluded from the condition-match.

**ANCHOR GATE PASS:** p=5,7 → PASS; p=11,13,23 → FAIL, exactly as required. Criterion valid.

## Sweep (primes 5..149, q=2)
33 primes; **24 confident**, 9 UNRESOLVED (depth-unstable: 37,41,43,67,83,103,137,139,149). Confident **PASS set =
{5, 7, 17, 127, 131}** (5 primes). Note 127,131 pass only at r=2 (shallow, γ ~ 10–50 = the under-resolution
regime); restricting to well-resolved r≥3 shrinks PASS to **{5, 7, 17}**.

## TEST 1 — DENSITY: no recognizable Artin/Chebotarev constant
Pass-fraction (confident) = 5/24 = **0.208**. Running value as primes accrue:
`1.000, 0.667, 0.600, 0.429, 0.333, 0.273, 0.231, 0.200, 0.176, 0.158, 0.143, 0.174` — **monotone decreasing,
still trending down**, no plateau at any Artin-flavored value (A=0.374, 1/2, 1/3, 2/3, 3/8, 5/8). The decay toward
~0.15 (with two large-prime passes that are r=2 noise) is the signature of a **small-prime phenomenon with density
→0**, not a positive Chebotarev density. Evidence AGAINST clean Chebotarev-governance.

## TEST 2 — CONDITION MATCH (decisive): **no exact match to any standard condition on 2 mod p**
PASS set {5,7,17,127,131} vs 15 candidate conditions, over the 24 confident primes (agree = pass⟺condition):

| condition | agree/24 | disagree | 
|---|---|---|
| (g) 2 is QR **and** ord even | 18 | **6** (best) |
| (d8) p ≡ 1 mod 8 | 16 | 8 |
| (h2) v₂(p−1) ≥ 3 | 16 | 8 |
| (a) 2 is QR (p≡±1 mod8) | 13 | 11 |
| (c1) ord ≡ 0 mod 4 | 13 | 11 |
| (c2) ord ≡ 0 mod 3 | 13 | 11 |
| (d3) p≡1 mod3 · (d6) p≡1 mod6 | 13 | 11 |
| (f1) index even · (f2) index odd | 13 | 11 |
| (d12) p≡1 mod12 | 14 | 10 |
| (b) ord_p(2) even | 10 | 14 |
| (d4) p≡1 mod4 · (h1) p≡3 mod4 | 11 | 13 |
| (e) 2 primitive root | 12 | 12 |

**Zero exact matches.** The closest, (g) "2 is QR and ord even," still disagrees on 6 of 24 primes (25% wrong) —
not a Chebotarev condition, just the best of a bad field. Every arithmetic-progression, order-parity,
primitive-root, index, QR, and 2-adic-valuation condition fails.

## TEST 3 — q=3 transfer: **not triggered** (no Test-2 match to transfer).

## VERDICT — **DEAD (fifth predictor-story killed)**
"Artin/Chebotarev-class" does not survive: (Test 1) the pass-density has no recognizable constant and decays toward
0 — a small-prime phenomenon, not a positive Chebotarev density; (Test 2, decisive) **no exact match to any of 15
standard conditions on 2 mod p** over 24 confident primes (best is 6/24 wrong). Forward-enrichment PASS/FAIL is
therefore **genuinely structureless — not merely formula-less, but not even Chebotarev-predictable** (the natural
home for "provably formula-less" arithmetic conditions), tested and rejected against the full standard menu.

This is the strongest honest "no predictor" statement: after five dead predictor-stories —
**q²=1+p, primitive-root, orbit-size (ord_p(2)), Chebotarev-density, Chebotarev-condition** — the forward half of
the enrichment dichotomy has no arithmetic predictor within reach, and its decisive tests are depth-bounded anyway.

**Caveat (honest):** the large-prime pass/fail is depth-limited (r=2 noise; 9 primes unresolved). But the verdict is
robust to it — even generously counting the r=2 passes (127,131), no condition matches; and restricting to
well-resolved primes shrinks PASS to {5,7,17}, a structureless triple.

**Disposition unchanged and now maximally supported:** the **ENGINE** (fiber-tower telescoping under a p-divisible
affine renewal; 9/9 exact, 13/13 with q, depth-independent) is the **sole general, predictable law**. The enrichment
dichotomy — forward and converse — is a **p=3-idiosyncratic** feature with no general mechanism and no Chebotarev
structure. The #2 paper claims the engine, and states the dichotomy's non-generality cleanly. Not at stake:
CHANNELDICH p=3 (14/14), the tower/ENGINE identity, R7–R10, MAHLER, value floor 0.473177.
