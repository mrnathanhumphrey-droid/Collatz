# RESULT — ORBIT-GATE: the ⟨2⟩-orbit block-diagonalization conjecture is FALSIFIED; the "+1" breaks ×2-equivariance (2026-07-28)

**Probe:** `probes/probe_orbit_gate.py`. Gate on Wilson's YOLO ⟨2⟩-graded unification conjecture — the proposed single
mechanism that would make the denominator theorem (`2^M−1`) and the c̃_q law (`2^{−ord}`) two corollaries of one statement:
*"the 2-adic renewal's transfer operator commutes with ×2, hence block-diagonalizes over the character group graded by the
⟨2⟩-orbit of each character (χ↦χ²), each orbit of size L contributing a `2^L−1` denominator."* Tested empirically at
G=(ℤ/q^k)\*, q=7 and q=13, k=2 (k=1 is degenerate: `qr+1≡1`).

## The gate — NEITHER grading block-diagonalizes K
Transfer operator K built in the character basis (`K_char = F† K F`); off-block `|K_char|²` mass measured under each
candidate partition of Ĝ, vs a random-partition null of matched block sizes:
```
 q=7,  k=2 (φ=42):   χ→χ² orbits: off-block 0.929   ×2-eigenclasses: 0.857   random null: 0.824
 q=13, k=2 (φ=156):  χ→χ² orbits: off-block 0.981   ×2-eigenclasses: 0.981   random null: 0.961
```
**Off-block mass is at the random-null level for both partitions** — K has *no* block structure under any ⟨2⟩-grading. (Too
high for block-triangular either, which would leave one triangle ~0.4–0.5.) **Candidate A (χ→χ²) and B (χ(2)-eigenclass) both FAIL.**

## Root cause (gated, not asserted): the +1 breaks ×2-equivariance
Any ⟨2⟩-orbit block-diagonalization requires `[K, P₂]=0` (P₂ = the permutation r↦2r). Direct check:
```
 ||K P₂ − P₂ K||_max   WITH the +1 (real qr+1 map):  0.500   |   WITHOUT it (pure halving r·2^{−v}):  0.000e+00
```
- The real qx+1 operator does **not** commute with ×2 (‖[K,P₂]‖=0.5). The pure-halving surrogate `r↦r·2^{−v}` (drop the whole
  affine `qr+1`) commutes **exactly** (0.0).
- Mechanism: from `r` the target is `(qr+1)·2^{−v}`; from `2r` it is `(2qr+1)·2^{−v} ≠ 2·(qr+1)·2^{−v} = (2qr+2)·2^{−v}`. The
  offending term is the constant **+1**. **The very "+1" that makes this the 3x+1 problem is exactly what kills the
  ⟨2⟩-equivariance the mechanism assumed.** It cannot be dropped — `qr+1` is what keeps the image in the units; `qr` alone
  leaves G. So the mechanism holds only for the pure-halving (non-Collatz) map.

## Verdict
- **The ⟨2⟩-orbit block-diagonalization conjecture is falsified** — not just χ↦χ², but the whole family (any grading needs
  `[K,×2]=0`, which fails by 0.5). **Corpse, cleanly diagnosed.**
- The empirical `2^{−ord}` (c̃_q) and `2^M−1` (denominator theorem) scalings are **still real** — but they are **not** explained
  by this mechanism. So they remain what the pre-YOLO gate said: a **rhyme / conceptual ⟨2⟩-in-a-multiplicative-group sibling**,
  **not** one theorem. `result_CTILDE_EXTEND.md`'s corrected framing stands; nothing here promotes it.
- **Falsifier #1 (the scary branch) is MOOT.** The "single full orbit ⟹ S_∞ rational" prediction rested on the block
  mechanism, which is dead. So there is **no reversal** of the irrationality lean — GARSIA (dimension-1/not-Garsia-singular),
  DENOM, and the i=20 no-crossing are **untouched**. 7/15 does not come back.
- **What survives:** the conceptual framing (both results are about `⟨2⟩` in a multiplicative group). A genuine unification
  would require a symmetry the **+1 respects** — the affine map `r↦qr+1` is ⟨2⟩-equivariant on *no* nontrivial grading here, so
  it is not the naive ⟨2⟩ action. That symmetry (if any) is unknown; this gate rules out the obvious candidate.

## Net
Wilson's own pre-YOLO caution was right: rhyme, not theorem. The YOLO attempt to promote it via character-orbit
block-diagonalization **fails the gate** — the `+1` breaks the ×2-commutation the mechanism needed. Either way we now know: no
false girder banked, and the irrationality lean is safe. **Not at stake:** c̃_q data, denominator theorem, GARSIA, DENOM,
SOLSTICE, R1–R30.
