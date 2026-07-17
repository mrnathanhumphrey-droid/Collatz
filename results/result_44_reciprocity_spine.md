# Result 44 (qx+1 paper) — the order-reciprocity SPINE does NOT transmit force. H_RIG FAILS: the 2-adic side speaks only in powers of 2; d_q=ord_q(2) generically doesn't. Claim A true, Claim B false.

**Date:** 2026-07-16. **Verdict: ✗ H_RIG REFUTED (as pre-registered). The (2,3) order-2/order-2 "reciprocity" is a POWER-OF-2 coincidence, not a coupling law. The pipeline is a true DESCRIPTION (four faces of Tao's measure), NOT a mechanically-coupled rig.**

**Headline: the spine posited ord_8(3)=2 (Chang) ↔ ord_3(2)=2 (Nathan) as reciprocal order-2 facts that, if a law, would make the four tools one rig and hand r_q a structural handle from the 2-adic side. It fails, and the reason is hard and structural, not "small-prime scatter": (Z/2^s)^* ≅ Z/2 × Z/2^{s−2}, so EVERY 2-adic element order is a power of 2. Nathan's d_q=ord_q(2) is generically NOT a power of 2 (d_7=3, d_11=10, d_13=12, d_19=18), so the order-matching reciprocity ord_{2^{s_q}}(q)=d_q is PROVABLY IMPOSSIBLE there. And ord_8(q)=2 is generic (q^2≡1 mod 8 for every odd q), so the mod-8 persistence constant is a FLAT 1/4 for all q — d_q-blind. The (2,3),(2,5),(2,17) matches are exactly the q where d_q=2,4,8 happen to be powers of 2 (reinforced at (2,3) by the Catalan pair 2^2−1=3, 3^2−1=8). Claim A [see the whole object] TRUE; Claim B [seeing ⇒ landing] FALSE. 0-for-10 on "and therefore." The machine working.**

Probe: `probe_44_reciprocity_spine.py`. Log: `result_44_reciprocity_spine_log.txt`. Runtime: instant (cheap arithmetic, the guard-mandated pre-check).

## The three structural checks

**(A) The reciprocity `ord_{2^{s_q}}(q) = d_q = ord_q(2)` exists only when `d_q` is a power of 2.**

| q | d_q = ord_q(2) | d_q pow-of-2? | ord_{2^s}(q), s=2..6 | s_q with ord=d_q? |
|---|---|---|---|---|
| 3 | 2 | ✓ | 2,2,4,8,16 | s=2 |
| 5 | 4 | ✓ | 1,2,4,8,16 | s=4 |
| **7** | **3** | **✗** | 2,2,2,4,8 | **NONE (impossible)** |
| **11** | **10** | **✗** | 2,2,4,8,16 | **NONE** |
| **13** | **12** | **✗** | 1,2,4,8,16 | **NONE** |
| 17 | 8 | ✓ | 1,1,1,2,4 | s=7 |
| **19** | **18** | **✗** | 2,2,4,8,16 | **NONE** |

`(Z/2^s)^*` has only power-of-2 element orders, so the 2-adic side can only ever realize a power-of-2 order. `d_q` is a power of 2 at q=3,5,17 (rare) and not at q=7,11,13,19 (generic). Where it isn't, the reciprocity **cannot exist** — no computation, a group-theoretic impossibility.

**(B) `ord_8(3)=2` is generic, not a q=3 reciprocal fact.** `q^2 ≡ 1 (mod 8)` for **every** odd q (all odd squares are 1 mod 8, since `(Z/8)^* ≅ Z/2 × Z/2` has exponent 2). So `ord_8(q)=2` for every odd `q ≢ 1 mod 8` (q=3,5,7,11,13,19 all give 2; only q=17≡1 gives 1). The "ord_8(3)=2" that looked like Chang's reciprocal to `ord_3(2)=2` is the universal odd²≡1 mod 8 fact — it carries no q-specific information.

**(C) The mod-8 persistence constant is flat in d_q.** `Pr[persistent]_q = Σ_k 2^{−k}·frac(μ: q^k·μ ≡ 7 mod 8)`. Because `q^k mod 8` has period ≤ 2 for every odd q, exactly one μ-class is persistent per k → `Pr = Σ_k 2^{−k}·(1/4) = 1/4` for **every** odd q. Computed: identical value for q=3,5,7,11,13,17,19. **Flat — it does not track `d_q`.** The brief's pass bar ("four constants on one monotone curve in d_q, hitting 1/4 at q=3") fails: the curve is a horizontal line, i.e. the constant is independent of d_q.

## Verdict on H_RIG — refuted, three ways, one reason

Every reading of the reciprocity fails, and all for the same structural cause — **the 2-adic multiplicative group speaks only in powers of 2, and `d_q = ord_q(2)` generically does not:**
- order-matching (`ord_{2^{s_q}}(q)=d_q`): impossible for non-power-of-2 d_q (q=7,11,13,19);
- fixed-modulus (mod 8): the constant is flat 1/4, d_q-blind;
- deeper 2-adic (mod 2^s): orders are still powers of 2, never a general d_q.

The apparent reciprocity at (2,3) — and at (2,5), (2,17) — is precisely the sparse set where `d_q ∈ {2,4,8}` happens to be a power of 2. At (2,3) it is additionally the **Catalan pair** (`2^2−1=3`, `3^2−1=8=2^3`, the only consecutive perfect powers), which is *why* the smallest instance looked like a law. It is not one.

## What this settles (the honest fork, brief §4)

- **Claim A — "I see the whole object" — TRUE.** The pipeline is real: Chang (2-adic input), Siegel (transform), Tao & Nathan (q-adic output) are four faces of one measure (Tao's Syracuse random variable on Z/3^n = our π). R42/R43 stand; the shared object is in the files.
- **Claim B — "seeing it whole ⇒ I can land it" — FALSE.** The order-reciprocity spine does **not** transmit force. The tools are aimed at one animal from four angles, but their loads do **not** sum through this spine, because the coupling was a power-of-2 coincidence, not a mechanism. `r_q` gains **no** structural handle from the 2-adic side via reciprocity.
- **The pipeline is a DESCRIPTION, not an active RIG.** L3 stays where it was: the q-adic output second-moment rate, to be bounded on its own terms (Kahane–Salem–Zygmund / the R42 renewal generating function `A(z)`), not derived through a 2-adic reciprocal law that does not exist.

In the metaphor: the net + reciprocal spear do **not** load Siegel's ground-to-air line. The mammoth keeps its altitude — for the same reason it always has: the escape axis is transcendental and no combination of ground tools climbs it via this spine. Claim A true, Claim B false, the 0-for-9 trap caught one more time (→ 0-for-10 on "and therefore"). The machine working — stop feeling, start knowing.

## Not at stake
R1–R43. This tests whether the order-reciprocity is a coupling law (it is not); it changes no r_q value, and it does not touch the (real) pipeline description or L3's statement.

_Reporting discipline: fired as the §5-guard-mandated structural pre-check (cheap arithmetic, read-only), not a heavy compute. H_RIG's failure was pre-registered by the brief itself; the added value is the SHARPER mechanism — a group-theoretic impossibility (power-of-2 orders in (Z/2^s)^*), decisive at q=7,11,13,19 — replacing the vaguer "small-prime scatter." The (C) constant's "vary" display label is a truncation artifact (sum to k=39 ≠ exactly 1/4); the values are identical across q = flat. No fit; exact orders and exact rational constant._
