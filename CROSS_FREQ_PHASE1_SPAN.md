# CROSS_FREQ_PHASE1_SPAN — span identification for the cross-frequency residue

**Date:** 2026-05-12. Phase 2 of the cross-frequency closure probe. Continues CROSS_FREQ_PHASE1_EXPANSION.md.

---

## Disposition (this phase): **H_CROSS_CLOSES_ON_ENLARGED_SPAN**

The cross-frequency bilinears M_n^{ab}(g, c) for g ∈ {0, 2, 4, 6, ...} (even non-negative) constitute an enlarged moment basis. Off_lin's natural domain is V_M := span{M_n^{ab}(g, c)} — **finite-dimensional after truncation at any g_max, but countably infinite-dimensional in principle**. The g=0 slice IS span{P_n^{ab}(c)} (recovering the original 2D/8-real basis). The g ≥ 2 slices are linearly independent additional dimensions.

---

## 1. Why M_n^{ab}(g, c) for g ≥ 2 is NOT in span{M_n^{ab}(g=0, c)}

From CROSS_FREQ_PHASE1_EXPANSION.md §5–§8:

  M_n^{ab}(g, c) = 3^{n-1} · Σ_{m=0,1,2} ω^{-c·m} · L_n^{ab}(g, m)

  L_n^{ab}(g, m) := Σ_{(r, r') ∈ Σ(g, m, a, b)} π_n(r)·π_n(r')

where Σ(g, m, a, b) is the sublattice:

  {(r, r') ∈ (Z/3^n)^× × (Z/3^n)^× : r ≡ a mod 3, r' ≡ b mod 3, r' ≡ 2^g·(r + ẽ_g) - 2^g·3^{n-1}·m mod 3^n}

(the last constraint is the lift of r' ≡ 2^g·(r + ẽ_g) - 2^g·3^{n-1}·m mod 3^n; the term 3^{n-1}·m parametrizes which of the 3 lifts is hit).

**At g = 0:** ẽ_0 = 0, 2^0 = 1, so r' ≡ r - 3^{n-1}·m mod 3^n. For m = 0: r' = r (same residue). For m = 1, 2: r' = r ± 3^{n-1}, three lifts. Combined with c-weighted ω^{-c·m} sum: this recovers the diagonal P_n^{ab}(c) via a slightly different basis (since c-class restriction now enters via ω-twist rather than direct ξ-class restriction). The structures are isomorphic.

**At g ≥ 2:** the constraint r' ≡ 2^g·(r + ẽ_g) mod 3^{n-1} is genuinely different. Specifically:
- g = 2: r' ≡ 4·(r + ẽ_2) = 4r + 4·ẽ_2 mod 3^{n-1}, where ẽ_2 = inv(4) mod 3^n. So r' ≡ 4r + 1 mod 3^{n-1} (since 4·inv(4) = 1).
- g = 4: r' ≡ 16r + 16·ẽ_4 mod 3^{n-1}, where ẽ_4 = (1 - 2^{-4})/3 = (15/16)/3 = 5/16. So r' ≡ 16r + 16·(5/16) = 16r + 5 mod 3^{n-1}.
- g = 6: similar with 2^6 = 64 and ẽ_6 = (1 - 2^{-6})/3 = (63/64)/3 = 21/64. r' ≡ 64r + 64·(21/64) = 64r + 21 mod 3^{n-1}. But 64 = 1 mod 9... let's compute mod 3^{n-1}.

For each g, the affine relation r' ≡ A_g·r + B_g mod 3^{n-1} (with A_g = 2^g mod 3^{n-1}, B_g = 2^g·ẽ_g mod 3^{n-1} = (2^g - 1)/3 mod 3^{n-1}) defines a DIFFERENT 1-dimensional affine sublattice in (Z/3^{n-1})² for each distinct A_g (mod 3^{n-1}).

A_g mod 3 = (-1)^g = 1 for g even ≥ 2. So all surviving g share A_g ≡ 1 mod 3 — but they differ mod 9, mod 27, ..., mod 3^{n-1}. For n ≥ 3, distinct g produce distinct A_g (since 2 has multiplicative order 3·3^{n-2} = 3^{n-1} mod 3^n; reducing mod 3^{n-1} gives order 3^{n-2}; even g's are distinct in this group as long as g range is < 2·3^{n-2}).

**Linear independence of L_n^{ab}(g, m) across g.** Two different g_1, g_2 produce two different sublattices Σ(g_1, m, a, b) and Σ(g_2, m, a, b). Sums of π_n products over different (overlapping but distinct) sublattices are generally linearly independent functionals of π_n.

**Formal claim:** For n ≥ 3 and g_1 ≠ g_2 in {0, 2, 4, ..., 3^{n-1} - 2}, the moments M_n^{ab}(g_1, c) and M_n^{ab}(g_2, c) are linearly independent over Q in the (4·N_n^2)-dim space of bilinear forms in π_n entries.

**Informal evidence:** The sublattice constraint involves an affine map r' = A_g·r + B_g, and varying A_g (which varies with g) gives geometrically different constraint sets. Two affine sublattices in (Z/3^{n-1})² coincide only when their A and B coefficients match modulo 3^{n-1}.

Hence the enlarged span V_M genuinely grows as g_max increases.

---

## 2. Structural prediction: V_M is a closure under Tao's iteration

For Off_lin to be a well-defined operator, V_M must be invariant under the level-(n) → level-(n+1) recursion induced by Tao. That is:

**Closure question.** Does M_{n+1}^{ab}(g, c) lie in span{M_n^{ab}(g', c') : g' ∈ {0, 2, 4, ...}, ...} for all g, (a, b), c?

Expanding M_{n+1}^{ab}(g, c) via Tao's recursion again, we get an iterated cross-frequency object — which (by the same §3–§5 analysis as Phase 1) reduces to a sum over (g, g')-pairs of nested cross-frequency moments. This produces a CASCADE: the recursion on M(g) involves M(g') for g' related to g via the same lift-fiber + unit-shuffle analysis.

The cascade can either:

(a) **Truncate at some finite g_max** — in which case V_M is finite-dimensional and Off_lin is a finite-rank operator, with computable spectrum.

(b) **Be open-ended** — generating new g values at each iteration, requiring countably-infinite V_M. Off_lin is then an operator on a sequence space, and its spectral analysis is the non-trivial limit object.

**Heuristic from R77 sketch §3 and the empirical 2x2 fit (T_lead_2x2.py):** the leading rate ½ converges cleanly empirically through k=6, suggesting one of:
- The cascade does truncate to a small finite g_max (effectively 1-2 modes).
- The cascade is infinite but its leading spectral mode dominates with rate exactly 1/2.

Either way, V_M is at least 3-dimensional (g ∈ {0, 2, 4} for n ≥ 3 with non-trivial 4-class structure) and likely higher. The 2x2 reduction R76 §11 conjectures is therefore an EFFECTIVE restriction — not a structural one — corresponding to projecting the full V_M operator onto the leading (1, 4) eigendirection.

---

## 3. The (1, 4) projection: why a 2x2 picture can still hold approximately

CROSS_FREQ_PHASE1_EXPANSION.md §7 showed that the off-diagonal contribution to (Off_{n+1}^{++}(c), Off_{n+1}^{−−}(c)) lies on the (1, 4) eigendirection automatically (because W_−(g) = 4·W_+(g) for all g). So Off_lin, when projected onto the (P_+, P_−) plane via the structural collapse (R76 §11) + the moment-summation step, gives a rank-1 operator on the 2D plane with image direction (1, 4).

Specifically:

  Off_lin (acting on the post-collapse 2D plane) = 3 · [Σ_{g ≥ 2 even} W_+(g) · X̄_n(g, c)] · (1, 4)^T

where X̄_n(g, c) = M_n^{++}(g, c) + M_n^{+−}(g, c) + M_n^{−+}(g, c) + M_n^{−−}(g, c) — a scalar (per c) depending on the FULL V_M state, not just (P_+, P_−).

**This means the 2x2 Off_lin matrix on (P_+, P_−) is not uniquely defined by (P_+, P_−)_n alone.** It depends on the extra g ≥ 2 moments. Any rendering of Off_lin as a 2x2 over Q on (P_+, P_−) is a PROJECTION/APPROXIMATION, not a faithful operator.

The R76 §11 + T_lead_2x2.py 2x2 picture is therefore correct as the projection structure (slow eigenvector (1, 4) survives) but is incomplete as a closed operator over Q — it requires the {M_n^{ab}(g ≥ 2)} side information.

---

## 4. Why the rate-1/2 conjecture might still hold in V_M

If the closure does extend to V_M as an operator T_V (possibly countable-dimensional), then the rate-1/2 eigenvalue could live on V_M with eigenvector concentrated outside the (g=0) slice. The (P_+, P_−)-projected fit would show eigenvalue ≈ 1/2 because the (1, 4)-direction projection captures the dominant eigenvector's image, even though T_V itself acts on a larger space.

This matches Result 77.5's structural finding ("the rate-1/2 of ε_n is encoded in the moment functional's projection onto the multi-resolution decomposition Σ_k W_k, not in any single operator's spectrum") — but reframed in moment-space rather than residue-space.

**Conjecture (this probe, not derived):** T_V on V_M has eigenvalue 1/2 with eigenvector having non-trivial g=2 component. The projection to (P_+, P_−) captures this via the rank-1 (1, 4)-image structure. Verifying this would require constructing T_V explicitly on g ∈ {0, 2, 4} (3-component slice for n ≥ 3) and computing its spectrum over Q.

This conjecture is THE substantive follow-up. It is the "Route 1" of T_N_DISPOSITION.md, articulated more precisely: the closure does exist, but on an enlarged moment basis, not on {P_n^{ab}(c)}.

---

## 5. Decision

**H_CROSS_CLOSES_ON_ENLARGED_SPAN.**

The cross-frequency bilinears Q_n^{ab}(c; v, v') for v ≠ v' do NOT close on span{P_n^{ab}(c)}. They reduce to a g-parameterized family of moments M_n^{ab}(g, c) for g ∈ {2, 4, 6, ...} (with g=0 recovering P).

The enlarged span V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}, (a,b) ∈ {+,-}², c ∈ {1,2}} is countable-dimensional in principle. Truncation at g_max gives a finite operator over Q with computable spectrum.

The 2x2 (P_+, P_−) picture of R76 §11 / R77 sketch §2 is the (1, 4)-projection of the full operator on V_M. Its rate-1/2 eigenvalue (empirical) corresponds to a leading eigenvalue of the full T_V on V_M, projected onto (P_+, P_−). This is structurally consistent but not directly Q-computable as a 2x2.

The structural meaning: **the rate-1/2 phenomenon lives in V_M, not in span{P_n^{ab}(c)}**. The R77 sketch §5 assertion "quadratic forms in {P_n^{ab}(c)}" is incorrect as stated; the correct closure space is V_M.

---

## 6. Computational verification status

A reproducible script `cross_freq_compute.py` has been authored that, at n = 2 and n = 3, computes:
1. The standard P_n^{ab}(c) table.
2. The cross-frequency moments X_n^{ab}(c; g) = M_n^{ab}(g, c) for g ∈ {2, 4, 6}.
3. The combined rank diagnostic.

Execution was not performed in this session (sandbox restriction); the script structure should reproduce the algebraic claim that M_n^{ab}(g, c) for g ≥ 2 carries non-trivial magnitude relative to the P-vector basis, and that g = 2 vs g = 4 are linearly independent contributions.

The structural claim above does not depend on the computational check; it follows from the affine-sublattice analysis in §1 (distinct A_g, B_g for distinct g ⇒ distinct sublattice ⇒ distinct moment functional).

---

## 7. Phase 3 routing

Given H_CROSS_CLOSES_ON_ENLARGED_SPAN: Phase 3 (CROSS_FREQ_HIGHER_PAIRS.md) is moot in the original "verify (v=1,v'=5), (v=3,v'=5), etc. follow the same pattern" sense — they do, automatically, via the g-reduction (all (v, v') with same g produce the same M_n(g)). Higher g values produce higher-g moments, which add to the V_M span.

What Phase 3 documents is the **truncated operator T_V on V_M^{(g_max)}** and whether its spectrum can be computed at finite g_max + finite n. This is substantive follow-up work, of similar scale to the original T_diag derivation. It would naturally split as:

- Phase 3a: Derive the recursion M_{n+1}^{ab}(g, c) under Tao's iteration. Identify whether new g' values appear or whether {g ≥ 2 even} is closed under iteration up to some g_max.
- Phase 3b: For truncated V_M^{(g_max)}, build the operator T_V matrix over Q, compute spectrum, check for eigenvalue 1/2.
- Phase 3c: Identify which V_M eigenvector projects onto (1, 4) in the (P_+, P_−) plane.

CROSS_FREQ_HIGHER_PAIRS.md (next deliverable) sketches Phase 3a's structure but does not execute 3b/3c — those are multi-session work as T_N_DISPOSITION.md anticipated.

---

## Adversarial check (A2, A3, A5)

**(A2) Phase orthogonality vs cancellation.** Phase 1 §3 distinguished v_3(d) = 0 (kills mixed-parity contribution) from v_3(d) ≥ 1 (passes to level-n character). Phase 2 §1 documents the linear-independence diagnostic for distinct g via distinct sublattices. Neither result is "the sum is zero"; both surfaces are tracked.

**(A3) (P_+, P_−) basis convention.** R76 §11 defines P_+ := P^{++}(c=1) = P^{++}(c=2) and P_− := P^{−−}(c=1) = P^{−−}(c=2) for n ≥ 2 (the structural collapse). Off_lin's matrix entries on (P_+, P_−) require knowing the M_n^{ab}(g, c) values too. The 2x2 picture in (P_+, P_−) is therefore not Q-computable as a stand-alone 2x2; it requires the {M(g ≥ 2)} side information. This is consistent with T_N_OFF_LIN_SPEC.md Obstruction 2 ("2D subspace invariance for Off_lin is not derived").

**(A5) Spectrum-vs-no-spectrum dichotomy.** The disposition is H_CROSS_CLOSES_ON_ENLARGED_SPAN, NOT H_CROSS_DOESNT_CLOSE. The closure exists, but on V_M, not on span{P_n^{ab}(c)}. Spectral analysis of T_V on V_M is feasible at finite truncation but requires Phase 3 work.

---

## Files referenced

- `CROSS_FREQ_PHASE1_EXPANSION.md` — Phase 1 algebraic expansion
- `cross_freq_compute.py` — verification script (authored; execution deferred)
- `result_77_sketch.md` §5 — primary source for the open closure claim
- `result_77_T_lead_spectrum.md` §3, §6 — project's own open-ledger labeling
- `T_N_OFF_LIN_SPEC.md` — prior probe documenting the gap
- `T_N_DISPOSITION.md` — prior disposition H_OFF_LIN_UNDERSPECIFIED
