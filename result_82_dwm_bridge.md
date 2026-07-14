# Result 82 — The DWM Chirp Bridge (disposition)

**Date:** 2026-07-14. **Verdict: H_BRIDGE_PARTIAL — structural common-form established; the EVIDENTIAL bridge is UNTESTABLE ON EXISTING DATA (DWM only at n=3 ↔ r≤2; R81 evidential floor is r≥3).**

Probe `probe_82_dwm_chirp_bridge.py`; data `result_82_data.csv`; log `result_82_log.txt`.

## Step 1 — common-form table

| field | char modulus D | multiplier κ | orbit variable | sign | measure |
|---|---|---|---|---|---|
| F̂ (Probe 81) | 3^{r+1} | c_{ℓ,ε} (a unit) | 2^{2u} in ⟨4⟩=⟨2²⟩ (index 2) | + | flat counting |
| DWM (step j) | 3^{n−2j+2} (effective) | 3^{2j−2} → unit 2^{−b} | 2^{−v} in ⟨2⟩ | − | Geom(½) 2^{−v} |

- **Modulus reduction verified** (`e_{3^n}(3^{2j−2}·ξ·w)=e_{3^{n−2j+2}}(ξ·w)`): PASS. This is the load-bearing fact: κ=3^{2j−2} is a *power of 3*, so it reduces the modulus rather than acting as a unit multiplier like F̂'s c.
- **Orbit / measure distinction verified** (index[⟨2⟩:⟨4⟩]=2): PASS. F̂ lives on the index-2 even-power sub-orbit ⟨4⟩ with flat counting measure (this is *what forces* the a≡1 mod 3 support and ord(4)=3^r); DWM walks the full ⟨2⟩ with Geom(½) weight 2^(-v). **Same character species, different sub-orbit + measure** — exactly the pre-reg's §2 framing.
- **Sign:** F̂ exponent +2u, DWM exponent −(…). Inversion 2^{−1} is an automorphism of ⟨2⟩ and the exponent sign is the additive-character negation ξ→−ξ. Recorded and verified, not assumed. No sign flip was *needed* to make forms agree (unlike the superseded 2026-05 DWM attempt).

Step-1 outcome: the two objects **can** be put in a common form (H_BRIDGE does not fail at step 1). The functional species matches; the differences (sub-orbit, measure, j-graded modulus) are exactly the documented content of the bridge.

## Step 2′ — the j-graded n↔r map (PRE-FIRE AMENDMENT, 2026-07-14)

The pre-reg guessed `n=r+1`. The modulus reduction above forces the finer, **j-graded** correspondence:

&nbsp;&nbsp;&nbsp;&nbsp;**r = n − 2j + 1.**

| n | j=1 | j=2 | j=3 |
|---|---|---|---|
| 3 | r=2 | r=0 | r=-2 |
| 4 | r=3 ✓evid | r=1 | r=-1 |
| 5 | r=4 ✓evid | r=2 | r=0 |
| 6 | r=5 ✓evid | r=3 ✓evid | r=1 |

**DWM data exists only at n=3**, whose steps land at r=2 (j=1) and r=0 (j=2) — **both below R81's evidential floor r≥3**. So the evidential bridge is **untestable on existing data**.

**The amendment changed the extension recommendation.** Under the naïve `n=r+1`, one would say "extend DWM to n=4 to reach r=3." The j-graded map shows n=4 makes *only the j=1 leg* evidential — and j=1 is the **exceptional step**: R3_DARK_SUBSPACE_STRUCTURAL.md, `x_1 = 2^{−b}` is a unit mod 3 (the unique W↔W^⊥ mixing event), while `x_{j≥2} = 3^{2j−2}·2^{−b} ≡ 0 mod 9`. Testing a claim about a **3²-phase object** on the one step where the mod-9 twist does not apply yields a near-uninterpretable number. **An evidential bridge requires n≥5 with j≥2** (n=6 puts both legs at r≥3), i.e. the state_count=162–486 regime (~6–160 h) — to be scoped as its own compute, with its own new Syracuse-side measurement and pre-registration, not drifted into from here.

**Step 2′ is itself a structural result, not bookkeeping.** `r = n−2j+1` says the DWM Kraus phase carries **3-adic depth n−2j+1 at step j** — the twist does not merely sit at 3²; it burns two units of 3-adic depth per level and **exhausts at j = (n+1)/2**. This upgrades R3_DARK_SUBSPACE_STRUCTURAL.md's `x_{j≥2} ≡ 0 mod 9` from a binary on/off threshold to the **j=2 instance of a level-graded depth sequence**: the mod-9 darkness is one term of `depth = n−2j+1`, not a wall. That is a statement about the trajectory's phase geometry, and it is what the modulus-reduction observation actually means.

## Step 5 — free triangulation: UNREACHABLE on existing data

Predicting the four reductions from the R81 J₄ indices independently requires overlapping (r≥3 evidential, matching-n DWM) data. No such overlap exists (DWM only at n=3 ↔ r≤2). Not reachable; deferred to the n≥5 extension.

## NON-EVIDENTIAL — FORMULA-LEVEL DIAGNOSTIC ONLY

Re-expressing the archived DWM operator as the common-form weighted chirp orbit sum (transcription check: PASS) and running the n=3 moments:

| V_MAX | moment | reduction | predicted | measured | ratio |
|---|---|---|---|---|---|
| 12 | G1 | sum_entries | +1.08427e-01 | 1.0783e-01 | 1.0055 |
| 12 | G2 | sum_entries | +6.08142e-01 | 6.0890e-01 | 0.9988 |
| 12 | G2 | tr_pi | +5.35700e-02 | 5.3570e-02 | 1.0000 |
| 12 | G2 | delta_1 | +5.74167e-02 | 5.7420e-02 | 0.9999 |
| 12 | G2 | vac_pi | +4.77296e-03 | 4.7750e-03 | 0.9996 |
| 16 | G1 | sum_entries | +1.07831e-01 | 1.0783e-01 | 1.0000 |
| 16 | G2 | sum_entries | +6.08879e-01 | 6.0890e-01 | 1.0000 |
| 16 | G2 | tr_pi | +5.35722e-02 | 5.3570e-02 | 1.0000 |
| 16 | G2 | delta_1 | +5.74203e-02 | 5.7420e-02 | 1.0000 |
| 16 | G2 | vac_pi | +4.77548e-03 | 4.7750e-03 | 1.0001 |
| 20 | G1 | sum_entries | +1.07820e-01 | 1.0783e-01 | 0.9999 |
| 20 | G2 | sum_entries | +6.08922e-01 | 6.0890e-01 | 1.0000 |
| 20 | G2 | tr_pi | +5.35746e-02 | 5.3570e-02 | 1.0001 |
| 20 | G2 | delta_1 | +5.74226e-02 | 5.7420e-02 | 1.0000 |
| 20 | G2 | vac_pi | +4.77574e-03 | 4.7750e-03 | 1.0002 |

**This is NOT evidence for H_BRIDGE.** It confirms only that the DWM operator was transcribed correctly — that the DWM machinery literally *is* a weighted chirp orbit sum (part (B) of the bridge, at the formula level). Whether that chirp is *structurally* R81's (part (A), evidentially) is exactly what n=3 cannot decide — the R81 lesson is that a thing can have the right shape and the wrong structure. The verdict does not cite this table.

## §6 — Igusa re-read against ĝ(a): two of three barriers are now stale

Re-reading IGUSA_DISPOSITION.md's three categorical barriers against the **R81 incomplete-sum object ĝ(a)** (not the R78 complete sum):

- **Barrier 1 (trivial substrate — R78 D=0 ⇒ g(u)≡c mod 3 ⇒ Z=1): STALE.** That was about the *complete* sum. ĝ(a) is the *incomplete* sum: flat nonzero magnitude 3√q (Th 78.3) and a certified nonzero 3-adic-analytic phase (Probe 81). The substrate is genuinely nontrivial; barrier 1 does not apply.
- **Barrier 2/3 (positive-irrational target vs negative-rational Igusa poles): INTRINSIC, still lethal.** ĝ(a)'s exponent scale is log₃4 = 2·log₃2, positive irrational (§ closing paragraph).

Net: the Igusa arc stays closed, but the *reason* is sharper — **two of the three barriers are stale (barrier 1 dies against ĝ(a); the substrate is now nontrivial), and the closure rests entirely on the one intrinsic barrier.** "Two of three stale, the third intrinsic" is a more useful map than "NO_FIT."

---

## The category obstruction (closing paragraph)

A positive-irrational exponent cannot be a pole of any machinery whose output is a finite set of rational poles. The load-bearing irrational is a single one: **log₃2** — the incommensurability of the 2-adic and 3-adic scales, which is the actual root of the obstruction and traces straight back to arc 4's archimedean finding. It is what surfaces as the R81 exponential chirp e_q(c·4^j) = e_q(c·exp₃(j·log₃4)): the chirp's scale **log₃4 = 2·log₃2 is the same irrational, not a second one.** Two arcs have proved this one number cannot be a rational pole from two sides: **arc 4 (ADELIC)** found log₃2 is not a Tate local-factor pole (Tate poles sit at s=0,1) and named the missing category as a non-Tate Mellin object; **arc 9 (IGUSA)** proved log₃2 cannot be an Igusa local-zeta pole (rationality + Monodromy + Bernstein–Sato force pole real-parts into the negative rationals). These are not independent negatives — they are **one categorical obstruction, stated once: the c=7/45 rate is encoded in an object whose natural singularity is the positive irrational log₃2 (a transfer-operator spectral radius or a branch cut), and no p-adic oscillatory-integral / algebraic-polynomial-Mellin machinery — Tate, Igusa, adelic — can produce it, because that entire class emits only finitely many rational poles.** (The T_lead gap log₃(45/43) is a *conjectural* third instance — almost certainly irrational, since 45/43 is not a power of 3, but its irrationality is not established here, so it is named as conjecture, not as a co-equal.) That is why Collatz resists this whole family of machinery, and the load-bearing sentence is the single one about log₃2.

_Reporting discipline: H_BRIDGE_PARTIAL is reported as fired. The n=3 reproduction is walled as non-evidential and not cited by the verdict. The j-graded amendment (step 2′) was made pre-fire and is what makes the untestability precise and what changed the extension recommendation. r=2 was not substituted for evidence._
