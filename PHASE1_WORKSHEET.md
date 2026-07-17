# PHASE 1 WORKSHEET — the five substrate lemmas
## Status: proofs written (Nathan); G1 adversarial check downgraded to by-inspection; numerical hygiene confirmed (probe_G1_smoke)

Setting: q odd prime, T_v(x) = (qx+1)·2^{−v} on ℤ/q^N (2 invertible),
p_v = 2^{−v}, K_N = Σ_v p_v (T_v)_*, lift(μ)(y) = μ(y mod q^k)/q,
π = stationary law of the IFS on ℤ_q, π_k its reduction mod q^k,
d_k := π_k − lift(π_{k−1}),  X_k := 3^k‖π_k‖²,  c_k := 3^k‖d_k‖².

── LEMMA 1 (FORGET). T_v(x) mod q^{k+1} depends only on x mod q^k.
PROOF. x′ = x + aq^k ⟹ qx′+1 = qx+1 + aq^{k+1} ≡ qx+1 (mod q^{k+1});
multiplying by the unit 2^{−v} preserves this. ∎
COROLLARY (descent + SPREAD). T_v induces T̂_v : ℤ/q^k → ℤ/q^{k+1},
and T̂_v is INJECTIVE (x ≠ x′ mod q^k ⟹ qx ≠ qx′ mod q^{k+1}), with
image 2^{−v}(1 + qℤ/q^{k+1}) = {y : y ≡ 2^{−v} mod q}.
[The pair forget-the-top-digit / spread-the-rest is the engine of
everything below.]

── LEMMA 2 (ONE-STEP). μ·K_{k+1} = π_{k+1} for EVERY probability
measure μ on ℤ/q^{k+1} with proj_k μ = π_k. In particular
π_{k+1} = lift(π_k)·K_{k+1}.
PROOF. Stationarity on ℤ_q reduces mod q^{k+1} (reduction commutes
with T_v), so π_{k+1}K_{k+1} = π_{k+1}. By Lemma 1, μK_{k+1} is a
function of proj_k μ only. proj_k π_{k+1} = π_k = proj_k μ, hence
μK_{k+1} = π_{k+1}K_{k+1} = π_{k+1}. ∎

── LEMMA 3 (INTERTWINE). For EVERY (signed) measure μ on ℤ/q^{k−1}:
      lift²(μ)·K_{k+1} = lift( lift(μ)·K_k ).
PROOF. Fix v and y ∈ ℤ/q^{k+1}; both sides vanish unless
y ≡ 2^{−v} (mod q) [Corollary: image condition on the left; the mod-q^k
solvability condition on the right — the SAME condition]. When it
holds: LEFT — T̂_v is injective, the unique preimage is
x_v(y) = (2^v y − 1)/q ∈ ℤ/q^k, contributing p_v·lift(μ)(x_v) =
p_v·μ(x̄)/q where x̄ := x_v mod q^{k−1}. RIGHT — the mod-q^k equation
pins x only mod q^{k−1} (= x̄), giving q lifts, each carrying
lift(μ)(x)= μ(x̄)/q; the outer lift contributes 1/q:
(1/q)·p_v·q·μ(x̄)/q = p_v·μ(x̄)/q. Equal term-by-term in v. ∎
[Verified numerically on random μ at 1e-17 before proving — the
proof explains the number, per protocol.]

── LEMMA 4 (REFINE). d_{k+1} = lift(d_k)·K_{k+1} — exactly, no
projection.
PROOF. lift(d_k)K_{k+1} = lift(π_k)K_{k+1} − lift(lift(π_{k−1}))K_{k+1}
= π_{k+1} [L2] − lift(lift(π_{k−1})K_k) [L3] = π_{k+1} − lift(π_k)
[L2 at level k] = d_{k+1}. ∎

── LEMMA 5 (PYTHAGORAS / RENEWAL). d_k ⊥ lift(π_{k−1}), and
      X_k = (3/q)·X_{k−1} + c_k        (exactly, every k ≥ 2).
PROOF. proj(d_k) = π_{k−1} − π_{k−1} = 0, so d_k has zero fiber-sums;
lift(π_{k−1}) is fiber-constant; such vectors are orthogonal. Hence
‖π_k‖² = ‖lift π_{k−1}‖² + ‖d_k‖² = ‖π_{k−1}‖²/q + ‖d_k‖² (lift
scales ‖·‖² by 1/q). Multiply by 3^k. ∎
[This is R7's Pythagoras / R42's renewal, now derived, not measured.]

── G1 STATUS (adversarial s≥2 / Wieferich check): DOWNGRADED TO
BY-INSPECTION. All five proofs are elementary algebra in q and never
touch τ, s, or d — they are s-blind and hold verbatim at q = 1093.
The Wieferich frontier enters only where the PHASES enter (Phase 3's
character sum), not in the substrate. Optional smoke test at moderate
q for hygiene; not load-bearing.
[CONFIRMED 2026-07-16: probe_G1_smoke verified all five identities to
machine precision at q=5,7,11 and L1 pointwise at q=1093 — see below.]

── PHASE 1 REMAINDER: none. Five lemmas stated and proved. Deliverable
for STATE: this worksheet + one dated entry. Phase 2 (entrance exam,
attempt three) now stands on theorems: prove the triple mode-collision
{1/3, 1/q, r_q/3} → 1/3 is FORCED at d = 2 and FORBIDDEN at d ≥ 3.

---
## Claude review note (2026-07-16)
Read all five proofs as reviewer before banking. Verdicts: L1 exact (q·q^k=q^{k+1}≡0); L2 correct
(T_v preimage pinned only mod q^k ⇒ μK depends on proj_k μ = fiber sums); L3 term-by-term in v with
both 1/q factors accounted; L4 chains L2/L3 correctly; L5 uses ‖lift μ‖²=‖μ‖²/q + zero-fiber-sum ⊥
fiber-constant. All correct. Numerical hygiene: probe_G1_smoke (see STATE entry).
