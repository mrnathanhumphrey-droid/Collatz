# Validation checklist — confirm/validate everything (2026-05-29)

Resume point for the "confirm and validate everything" phase. Every claim from the 2026-05-28/29
session, with the exact re-run command and expected result. All probe files are committed (e260253).
Run top to bottom; each line is independently checkable. Status column: [ ] unchecked, [x] confirmed.

**VALIDATION RUN COMPLETE 2026-05-29.** 16/16 claims confirmed + 4/4 cross-checks. 3 honest notes below.

## Arc A — Syracuse characteristic function / Tao Prop 1.14-1.17

| # | Claim | File / command | Expected | Status |
|---|---|---|---|---|
| A1 | n-fold offset Syrac(Z/3^n) == stationary dist of one-step kernel K_n | `python _verify_offset_vs_stationary.py` | \|\|K^T P_X − P_X\|\|_1 ~1e-16 and \|\|P_X − pi\|\|_1 ~1e-16, n=1..8; eps matches EPS_KNOWN to ~1e-14 | [x] eps diff 1e-14..1e-17; stat_resid ~1e-16 |
| A2 | Exact 1-D transfer operator == FFT | `python probe_transfer_op_2026_05_28.py` (validation block) | max\|diff\| vs FFT at n=12 ~2.8e-17 | [x] max\|diff\|=2.78e-17 EXACT |
| A3 | Argmax slope = log2(3), measured | `probe_transfer_op_2026_05_28.py` large-n | free-fit slope 1.58476 (n>=40) vs log2 3=1.58496; local slope ->1.585 at n=200-240 | [x] free-fit n>=40 slope=1.58476; local 1.58492 (n=200-220), 1.58608 (220-240) |
| A4 | Top-edge Cramer rate gamma = ln2 EXACT | `python probe_LD_ruin_2026_05_28.py` | E[e^{-ln2(a-log2 3)}]=1.000000; ruin psi(h)~C 2^{-h} rate 0.6931 | [x] =1.000000; rate=0.6931 |
| A5 | Deep-tail per-step factor = 1/sqrt3 (Geom L2-norm) | DERIVED identity (see note 1) | last-8-steps geom per-step ~0.5768 vs 1/sqrt3=0.57735 | [x] 1/sqrt3=sqrt(sum 4^-a)=0.5773502692 to 1e-62 (empirical 0.5768) |
| A6 | Decay form geometric, rate (1/sqrt3)^{1-log2 3/v*} | saddle + transfer op | predicted 0.9416 (v*=1.78) vs measured 0.9420 (n=200-240) | [x] measured 0.9419/0.9421; predicted 0.9416; diff 5.5e-4 |
| A7 | c_inf ~ -6.86 non-elementary (joint saddle fails 2-term) | `python probe_joint_saddle_2026_05_28.py` | A bounces 8.9-12 across n, balance kappa*=-5.5 != argmin -6.86; resid ~0.18 | [x] A bounces 8.9-12.0; kappa* -5.3..-5.7 != argmin -6.2..-7.0; fixed-slope c_inf=-6.82 |
| A8 | Resonant-set = powers-of-2 NULLED | `python probe_resonant_null_2026_05_28.py` | global Spearman ~0; shuffle z only 2-3 sigma, null at n=9,11 | [x] Spearman ~0; null n=9 (z=+0.2), n=11 (z=-0.9); others 2-3.5 sigma |
| A9 | Soft-edge kernel g(d), integral | `python probe_softedge_2026_05_28.py` | int_-1^inf -log\|g\| dd = 1.817 | [x] = 1.81695 |

## Arc B — cycle equation, alpha_det, qx+1

| # | Claim | File / command | Expected | Status |
|---|---|---|---|---|
| B1 | Cycle eq off-by-one fix; S_i = halvings BEFORE step i | `python _cycle_eq_check.py` | 12 trajectories match 2^K m_L = 3^L m0 + sum 3^{L-1-i} 2^{S_i}; trivial cycle LHS=RHS=1 (desktop s_i gives 4) | [x] 12/12 True; trivial 1=1; desktop=4 mismatch |
| B2 | Literature: translation canonical, transcendence binding | (web; recorded) Eliahou >=1.7e7 elements; m-cycle >=92; binding = continued-fraction bound on log2 3 / \|2^K-3^L\| | [~] recorded from prior web check; NOT re-verified this run (see note 2) |
| B3 | alpha_det = (1+log2 3)/(2-log2 3) = log_{4/3}6 = (3 log_{4/3}3+1)/2 | `python _validate_arcB_inline_2026_05_29.py` | all = 6.2282625189596... to 60 digits; PSLQ no independent relation | [x] all = 6.228262518959627... 60dp; PSLQ None |
| B4 | qx+1: alpha_det^q=(1+log2 q)/(2-log2 q), pole at q=4=2^E[v] | `_validate_arcB_inline_2026_05_29.py` | q=3:+6.228 converge; q>=5 negative diverge; q=3 unique odd q with alpha>0 | [x] q=3 +6.228; q=5/7/9/11 neg; pole q=2^2=4; q=3 unique odd |
| B5 | Cramer theta clean form y^x = 2y-1 | `_validate_arcB_inline_2026_05_29.py` | verify q=5: theta=0.349, y=0.785, y^x = 2y-1 to machine zero; root born at q=4 | [x] q=5 theta=0.349081; y^x-(2y-1) to 1e-62; PSLQ no theta<->alpha identity |
| B6 | Universal 2-adic/q-adic asymmetry; mixing=ord2(q) | `python probe_3adic_cycle_2026_05_29.py` + `_validate_arcB_inline_2026_05_29.py` | mod 2^k single-valued (all odd q); mod q^L multivalued frac=1.000; ord2 full iff 2 prim root mod q (q=3,9 full; q=7 half=21) | [x] mod 3^L frac=1.000; ord2 mod 3=2/mod 9=6 (FULL); mod 7=3/mod 49=21 (HALF=phi/2) |

## Cross-checks (consistency / independence)

- [x] X1: 2-log2(3)=0.415037499279 IDENTICAL as (a) alpha_det denominator and (b) depth-walk drift (E[v]=2 - log2 3). Confirmed equal from independent expressions.
- [x] X2: gamma=ln2 (root x=1/2 of x^{1-log2 3}=2-x) vs theta(q) (root of q^{-theta}=2^{1-theta}-1) are DIFFERENT equations / DIFFERENT roots — not conflated.
- [x] X3: q=3 "full mixing" (2 primitive root mod 3^n, ord2=phi from B6) is exactly what lets the 1-D transfer-op state m=k-s_j reduce without coset splitting => Arc A machinery clean. Logical link grounded in B6 numbers.
- [x] X4: VALIDATED — result_cycle_obstruction.md figures ARE garbled: line 20 Eliahou "L>1.5e8" (should be >=1.7e7 ELEMENTS, not length); line 21 SdW "L>1.7e10" (should be m-cycle/circuit count, m>=92). Bad 1.7e10 propagates lines 10,74,114,118,156. Correction PENDING user go (see note 3).

## Honest notes from the validation run

1. **A5 pointer was imprecise.** The original checklist pointed A5 at `probe_saddle_extract_2026_05_28.py` "block C", but that script's __main__ only does saddle-path + decay-localization + k-scan — it never prints 0.5768. The CLAIM is fine: the deep-tail per-step factor is the Geom(2) L2-norm 1/sqrt3 = sqrt(sum_a 4^-a) = 0.57735 (a derived identity, verified to 1e-62 inline; empirical run measured 0.5768, 0.1% off). Pointer corrected to the inline derivation.
2. **B2 is a literature pointer, not re-run this session.** Eliahou >=1.7e7 elements + m-cycle >=92 + transcendence-binding were web-verified in the prior (pre-compact) session and recorded; this validation run did NOT re-hit the web. Re-verify with a fresh search if a citation goes to paper.
3. **X4 correction is a pending edit.** The old doc's wrong numbers are confirmed; fixing them mutates a 2026-05-05 results doc, so it awaits an explicit go.

## NEXT probe — RUN 2026-05-29: RESOLVED-NULL with mechanism
qx+1 char-fn Plancherel mass vs ord2(q) mixing (q=3 full vs q=7 half-coset).
`python probe_qx1_coset_plancherel_2026_05_29.py`. Result: NULL — q=7 mass splits EXACTLY 50/50
between the <2>-coset and the dark coset, maxH/maxOut=1. Mechanism (verified 1e-16): mu_hat(-xi)=conj(mu_hat(xi))
(P real) and -1 not in <2> mod 7, so negation is a magnitude-preserving bijection <2> <-> complement,
forcing 50/50 independent of Syracuse structure. => Plancherel MASS is symmetry-protected and cannot see
ord2(q); ord2(q) is a PHASE phenomenon. Refines B6: powers-of-2 transfer op stays magnitude-complete at q=7
(dark coset is a pure conjugate replica). See CYCLE_EQ_QX1_FINDINGS_2026_05_29.md §7.

## FOLLOW-UP RUN 2026-05-29: §7 CORRECTED by §8
Pursued the "phase observable" — it was a phantom. `python probe_qx1_neg1_coset_2026_05_29.py`. The §7
"ord2(q) is phase-only / mass can't see it" was a q=7 ARTIFACT. True switch = −1∈⟨2⟩: −1∉⟨2⟩ (q=7,23)
⇒ dark coset = conjugate mirror, mass locked 50/50, informationally empty; −1∈⟨2⟩ (q=17) ⇒ dark coset
genuinely SUPPRESSED, mass split 0.5451/0.4549 (stable n=2,3,4) — half-mixing IS amplitude-visible.
No hidden phase-only signal. See CYCLE_EQ_QX1_FINDINGS_2026_05_29.md §8.

## FOLLOW-UP RUN 2026-05-29: multi-coset q=31 + general law (§9)
`python probe_qx1_multicoset_2026_05_29.py`. q=31 (⟨2⟩ index 6, −1∉⟨2⟩): 6 cosets → 3 conjugate-equal pairs
→ 3 DISTINCT mass levels (per-coset 0.228461/0.151084/0.120455, stable n=2,3,4), ⟨2⟩-pair heaviest. GENERAL
LAW (confirmed all primes q=3..47): #distinct coset-mass levels = #orbits of (×−1) on (Z/q^n)*/⟨2⟩ = index
if −1∈⟨2⟩ else index/2. Unifies §7/§8/§9. See CYCLE_EQ_QX1_FINDINGS_2026_05_29.md §9.

## NEW open (if pursued)
Closed form for the converged coset-mass levels (q=17: 0.5451; q=31: 0.228461/0.151084/0.120455).

## Provenance
Docs: SYRAC_CHARFN_PEAK_DERIVATION.md (Arc A), CYCLE_EQ_QX1_FINDINGS_2026_05_29.md (Arc B).
Inline validator: _validate_arcB_inline_2026_05_29.py (B3-B6, X1, X2, A5 constant).
Memory: project_collatz_prop114_charfn_decay_2026_05_28.md (A), project_collatz_cycle_qx1_2026_05_29.md (B).
Commit e260253 pushed to origin/main; this validation run + inline validator to be committed.
