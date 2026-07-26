# RESULT — PROBE BRIDGE: both gates fail; the root cause is the dlog coordinate (2026-07-26)

**Probe:** `probe_bridge.py` (+ dlog inspection). Wilson's run order: (1) additive recursion on the FORWARD measure,
(2) the Gauss-sum bridge. He asked to know immediately if either fails.

## Sanity — the forward measure is correct
`stationary_trunc(3,k)` (the q-sweep forward-Syracuse chain) reproduces R66's max exactly: 0.14283, 0.06362, 0.03133,
0.01671, 0.00924 for k=2..6 = R66 to 1e-4. **So the forward chain = R66's/your W measure — the measure is right.**

## Gate 1 — the additive recursion STILL fails on the forward measure
`μ̂_n(ξ) = Σ_a 2^{−a} e(ξ2^{−a}/3^n) μ̂_{n−1}(ξ2^{−a} mod 3^{n−1})`, μ̂ = fft(forward), n=3..6: **rel ~1.33–1.41.**
So the perpetuity-vs-forward offset was *not* the (only) issue — it fails on the forward measure too. This is pure
additive DFT (no dlog), so the failure is in the **level bookkeeping** of the recursion: `stationary_trunc(3,n)` is
the *fixed point* at level n, and the fixed-point self-consistency is at a single level with a 3× frequency scaling
(`μ̂(ξ)=Σ_a 2^{−a}e(…)μ̂(3ξ2^{−a}/…)`), not a clean level-n-from-level-(n−1) relation. The n↔n−1 form as written
doesn't hold on the stationary measures.

## Gate 2 — the bridge fails, and the reason is concrete: `|τ| ≠ √q`
`ρ̂(a) = (1/τ(χ̄_a)) Σ_t χ̄_a(t) μ̂(t)`: rel ~1.8–2.1, and the Gauss-sum modulus is wrong — `|τ|` = 4.90, 9.85, 20.69
at k=3,4,5 vs `√(3^k)` = 5.20, 9.00, 15.59. A primitive character mod 3^k has `|τ|=3^{k/2}`; mine doesn't, so my
`χ_a` is **not a proper multiplicative character**.

## The root cause — R10's "dlog" is NOT the group discrete log
Inspecting `R10.dlog_table`: it's a **permutation of {0,…,3^k−1}** — 3^k distinct values. But `(ℤ/3^k)*` is cyclic of
order `φ(3^k) = 2·3^{k−1}`, and its discrete log lands in `ℤ/(2·3^{k−1})`, not `ℤ/3^k`. **So the codebase's dlog is a
bespoke 3^k-point coordinate (the one where ×4 acts as +1 — the channel coordinate from CHANNEL_ID), not the
group-theoretic discrete log.** Consequences:
- The channel object `fft(ρ)` is a well-defined additive DFT over *this* 3^k coordinate (MAXMODE2's U² identity holds
  there) — but it is **not** the multiplicative character transform of `(ℤ/3^k)*`.
- The Gauss-sum bridge (standard `χ = ` additive-character-sum identity) is a statement about `(ℤ/3^k)*` characters
  living in `ℤ/(2·3^{k−1})`. It does **not** connect to `fft(ρ)` in R10's 3^k coordinate — the two "multiplicative"
  objects are different, which is why `|τ|` is wrong and the bridge misses.

## Verdict — the bridge as written doesn't connect the banked channel object; convention is yours to pin
Both your gates fail, and the diagnostics are clean:
1. The recursion's **level indexing** (fixed-point single-level self-consistency vs n↔n−1) needs restating.
2. The bridge assumes standard `(ℤ/3^k)*` characters (index space `ℤ/(2·3^{k−1})`), but the channels live in R10's
   3^k dlog coordinate. **The additive↔multiplicative bridge only closes if the channels are recomputed in the
   group's discrete-log coordinate (`ℤ/(2·3^{k−1})`), or if the R10↔standard-dlog map is supplied.**

Per the reconstruction guardrail (three failed reproductions now — perpetuity, forward, bridge), I stop here rather
than guess the coordinate. What unblocks it, concretely: the standard discrete log of `(ℤ/3^k)*` w.r.t. a fixed
generator (index space `ℤ/(2·3^{k−1})`, the `ℤ/2 × ℤ/3^{k−1}` split — the `ℤ/2` being your `(−1)^a` parity). With
that coordinate, `ρ̂` becomes a genuine multiplicative transform and the Gauss-sum bridge is a machine-precision
identity. It's a real object mismatch, not a cleverness gap — the hunt found the exact seam where the two coordinates
fail to line up.

**Not at stake:** the channels/MAXMODE2 (self-consistent in R10's coordinate), MEAN1, HIERARCHY, CHANNEL_ID, R1–R30.
Cheap (stationary_trunc, 0.1s).
