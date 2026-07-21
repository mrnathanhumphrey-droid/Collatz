# Probe R7 — the channel engine (the ⟨4⟩-orbit law) — **FULL GATE PASS (k=2,3,4,5)**

**Date:** 2026-07-21  Exact rationals, minutes. Probe `probes/probe_engine_R7.py`. Gates Wilson's closed-form
off-diagonal **engine**: the ledger reorganized from |v−v′| bins (R6-C) into the ⟨4ᵐ⟩-orbit law

> **OffDiag_k = (2/3) Σ_{m≥1} 4^{−m} C_k(m)**,  C_k periodic in m mod 3^{k−1} = ord_{3^k}(4),
> **C_k(m) = Σ_{a,a′} μ_{k−1}(a) μ_{k−1}(a′) · c_{3^k}( 4^{−m}(1+3a) − (1+3a′) )**  (c = Ramanujan sum).

**The engine was built from first principles, then checked against the frozen increments — not read off the
S-table.** The channel weight (2/3)4^{−m} is derived: even gap g=2m, both orderings, Σ_{v′≥1} 2·2^{−(v′+2m)}2^{−v′}
= (2/3)4^{−m}; the inner unit factor 2^{−v′} is a unit mod 3^k so c_{3^k}(2^{−v′}δ)=c_{3^k}(δ) — the correlation
depends on the gap-index m alone. The measure hierarchy μ_k is the renewal X = 2^{−v}(1+3·a), a∼μ_{k−1},
v∼Geom(½), built exactly (folded valuations); |supp μ_k| = 2, 6, 18, 54 for k=1..4; μ₂ reproduces R6's pi2.

## R7-A — the engine at k=2 and k=3: **GATE PASS**
| quantity | engine | pre-registered | verdict |
|---|---|---|---|
| C₂(m), m=1,2,3 | **(−1, −1, +2)** | (−1, −1, +2) | **PASS** |
| OffDiag₂ | **−4/21** | −4/21 | **PASS** |
| C₃(1) | **4/49** | (forced) | — |
| m=1 channel (2/3)4⁻¹·C₃(1) | **+2/147** | R6-C gap-2 = +2/147 | **PASS — sign-flip DERIVED** |
| OffDiag₃ | **−2980/203889** | −2980/203889 (frozen incr.) | **PASS** |

**The gap-2 sign flip −1/6 → +2/147 (k=2→3) is now derived, not observed.** Both are the single m=1 channel
(2/3)4⁻¹·C_k(1) = (1/6)C_k(1): at k=2, C₂(1)=−1 → −1/6; at k=3, C₃(1)=+4/49 → +2/147. **The flip is C_k(1)
crossing zero** as the measure deepens (μ₁→μ₂ at the same orbit position m=1). The mechanism has a name: the
orbit-correlation character sum at m=1 changes sign with k.

## R7-B — the C-tables (exact, period 3^{k−1}) — the limit law's raw material and its judge
| k | P=3^{k−1} | C_k(m), m=1…min(P,9) |
|---|---|---|
| 2 | 3 | −1, −1, **+2** |
| 3 | 9 | 4/49, −26/49, −5/7, 22/49, 22/49, −5/7, −26/49, 4/49, **10/7** |
| 4 | 27 | 2588172/94264681, 1505640/94264681, −19158/67963, 17644728/94264681, −20771232/94264681, … |
| 5 | 81 | 595374766986694714604534052/49461751815032265013619722849, … (81 classes) |

**Structural facts read off the tables (for the pen's limit law C_∞):**
1. **Palindrome:** C_k(r) = C_k(3^{k−1}−r) — verified exactly at k=3 (C(1)=C(8)=4/49, C(2)=C(7)=−26/49,
   C(3)=C(6)=−5/7, C(4)=C(5)=22/49). This is the m→−m (gap-sign) symmetry: c is even, orbit-correlation is
   real-symmetric. The tail law only needs the first half plus the self-paired DC class.
2. **DC/self-orbit class is positive and large:** C_k(3^{k−1}) — the m≡0 class — is **+2** (k=2), **10/7** (k=3);
   it is the ⟨4⟩-orbit's "diagonal" term (4^{−m}≡1) and the only class that is not conjugate-cancelled. It is
   what the negative interference channels must overcome — and don't (net negative through k=3, then the balance
   tips as the weight 1/(4^P−1) of this class collapses).
3. **Mersenne denominators arrive as steered:** the class weights carry 1/(1−4^{−P}) = 4^P/(4^P−1), the (4^P−1)
   family — 63=4³−1 (k=2), 262143=4⁹−1 (k=3, = denominator scale of the 94264681 = 9707²… the k=4 shell). 21 =
   63/3 is −4/21's denominator; representation (Ramanujan/Plancherel) and mechanism (⟨4⟩-orbit geometry) meet.

## R7-C — the running ledger (engine vs frozen; target −1/5)
| k | engine OffDiag_k | float | frozen S_k−S_{k−1} | match | running Σ | vs −1/5 |
|---|---|---|---|---|---|---|
| 2 | −4/21 | −0.190476 | −4/21 | ✅ | −0.1904762 | +0.009524 |
| 3 | −2980/203889 | −0.014616 | −2980/203889 | ✅ | −0.2050920 | −0.005092 |
| 4 | +5699915795296300/… | +0.002640 | (frozen) | ✅ | −0.2024523 | −0.002452 |
| 5 | +6958280182844849…/… | +0.001301 | (frozen) | ✅ | −0.2011517 | −0.001152 |

**Engine = frozen at every k = 2,3,4,5.** Σ_{k≥2} OffDiag → −1/5 (the constant's target); Σ_{k≥3} → −1/105.
The overshoot (−,−,+,+) is now inside the engine: k≥4 channels turn net positive because the DC class weight
1/(4^P−1) collapses super-geometrically while the negative bulk channels (C(2),C(3),…) also shrink — the
positive DC term stops being outrun.

## Guard — odd-gap channels vanish (conjugate-kill; verified, not assumed)
Odd gaps g=1,3,5,7 give character sum **0 at both k=2 and k=3** (all four, exactly). Confirms the weight ledger's
"odd gaps 4/9 conjugate-killed": the interference lives entirely in the even-gap 2/9-weight channel system, as
the engine assumes. The diagonal 1/3 (replication) + odd 4/9 (killed) + even 2/9 (live) = 1 ledger is complete.

## Status
**R7 FULL GATE PASS.** The channel engine — derived independently (channel weight (2/3)4^{−m}, twisted
⟨4ᵐ⟩-orbit character sum C_k(m) of μ_{k−1}) — reproduces OffDiag_k for **k=2,3,4,5** exactly, matching the frozen
S_k−S_{k−1}. R6-C's gap-2 sign flip (−1/6→+2/147) is **derived** as C_k(1) crossing zero (C₂(1)=−1, C₃(1)=4/49).
The C-tables (R7-B) are the raw material for the pen's limit law C_∞(m): palindromic, with a positive DC/self-orbit
class of collapsing weight 1/(4^P−1) (Mersenne) racing the negative bulk channels — the machine behind the −,−,+,+
overshoot and the −1/5 total. **Still owed (pen):** the stationary profile C_∞(m) and Σ_k (2/3)4^{−m}C_k(m) → −1/5
in closed form (the L→∞ limit of Theorem S). No fitting; exact rationals; the ledger is finitely computable term
by term, forever, and now every term through k=5 is engine-derived.
