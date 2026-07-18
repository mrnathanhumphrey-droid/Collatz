# PATH2 — Family-level extension of R78.4-78.6 to general prime p ≥ 3

**Status:** Phase 1 deliverable. Family-level statements; q=3 specifics flagged.

## Notation (locked in PRECISE_ASK.md style)

For prime p ≥ 3, integer r ≥ 1, c ∈ (Z/p^{r+1})*:
- `q := p^{r+1}`, `period := p^r`, `N := p^{r-1}`
- `M := p^{r+1}` (same as q; both names used in different docs)
- `f_p(u) := e_q(c · (1+p)^u)`, period p^r in u
- `F̂_p(ξ) := Σ_{u=0}^{q-1} e_q(c·(1+p)^u − ξ·u)` (full-period DFT)
- `1̂(ξ) := Σ_{u=0}^{N-1} e_q(ξ·u)` (Dirichlet kernel, length N)

The R78.4-78.6 statements at q=3 (`result_78_extended.md` lines 9-27):

> **(78.4):** F̂(3a) = 3 · e_q(1) · G(a), G(a) = Σ_{s=0}^{period-1} e_q(P_a(s)), P_a(s) = 3s − C_a·L(1+3s), L(1+3s) = Σ_{j=1}^J (-1)^{j-1}/j·(3s)^j, C_a = a·L̃^{-1} mod 3^r, L̃ = L(4)/3.
>
> **(78.5):** Bijection a ↔ C_a on {a ≡ 1 mod 3 in Z/3^r}.
>
> **(78.6):** At r=3, ψ(a) = G(a)/√q = e_q(P_a(s*(C_a))), s*(C_a) = (C_a−1)/3 mod 3.

## Family-level analogs (Phase 1 derivation)

### Setup at general p

The principal-unit structure used:
- (1+p) generates a cyclic subgroup of (Z/p^{r+1})* of order p^r (for p ≥ 3 odd; fails at p=2 — FHAT verification doc §4).
- So f_p(u) = e_q(c·(1+p)^u) has period p^r in u.

The Cochrane Prop 4 truncated p-adic log applies at any odd prime — it's a formal identity in `Z_p[[X]]`:
- `log(1+px) = Σ_{j≥1} (-1)^{j-1}/j · (px)^j`
- Series converges p-adically (each term has v_p ≥ j − v_p(j) ≥ j − log_p(j))
- Truncated to J terms with J = `J_{p,1,r+1}` = max j such that `j − v_p(j) < r+1`. Then L_p(1+px) mod p^{r+1} captures the full log mod p^{r+1}.

**Computing J_p:** For p=3, J grows roughly linearly with r. For p≥5, J grows faster (v_p(j) is rarer). Specifically:
- p=3, r=2 → J=3 (J_{3,1,3}=3)
- p=3, r=3 → J=3
- p=5, r=2 → J=3 (since 5−v_5(5)=4 ≥ 3, but 4−v_5(4)=4 ≥ 3; check: 3−v_5(3)=3, want < 3? No, want < r+1=3. So 3−0=3 NOT < 3, fails; 2−v_5(2)=2 < 3, so J=2. Actually: J = max j with j − v_p(j) < r+1, so 2−0=2 < 3 ✓, 3−0=3 NOT < 3 ✗; J=2.)
- Let me recompute via the code in path_B_explicit_phase.py: `J_for_p3(m): j+1 with (j+1)−v_3(j+1) ≥ m`. So J is max j such that for ALL k ≤ j, k−v_p(k) < m. For p=5, m=3: j=1 (1−0=1 < 3 ✓), j=2 (2−0=2 < 3 ✓), j=3 (3−0=3 NOT < 3, fails). So J=2.

(This matters for the explicit polynomial form. **It is NOT a structural barrier** — J just sets the polynomial degree.)

### Family Theorem 78.4_p (Explicit Gauss-sum factorization)

> For prime p ≥ 3, r ≥ 1, c ∈ (Z/p^{r+1})*, a ≡ c (mod p) in Z/p^r:
>
> **F̂_p(p·a) = p · e_q(c) · G_p(a)** where **G_p(a) = Σ_{s=0}^{p^r − 1} e_q(P_a(s))**
>
> with `P_a(s) = p·c·s − C_a · L_p(1+ps)`, `L_p(1+ps) = Σ_{j=1}^{J_p} (-1)^{j-1}/j · (ps)^j` mod q, `C_a = a · L̃_p^{-1}` mod p^r, `L̃_p = L_p(1+p)/p` (the unit mod p^r after stripping single p-factor from L_p(1+p)).

**Why this generalizes p-blindly:**

1. **Cochrane Prop 4 setup.** Prop 4 holds for any odd prime: a character χ on principal units 1+p·Z_p is determined by χ(1+p) (cyclic generator), and χ(1+ps) = e_q(−C_χ · L_p(1+ps)) where C_χ depends linearly on χ. The structural statement is p-blind.

2. **Variable substitution.** F̂_p(p·a) sums e_q(c(1+p)^u − p·a·u). Reindex u → s via (1+p)^u = 1+ps mod q for some s = s(u) (the principal-unit parametrization). The Jacobian of u ↔ s in this parametrization gives a factor of p (in q=3 case, the factor 3 in front of G(a); same factor p appears at family level).

3. **Character pairing.** The χ_a character defined by χ_a(1+p) = e_period(−a) (so χ_a(4) at p=3 becomes χ_a(1+p) at family p) pairs against c · (1+p)^u to give the C_a coefficient via Cochrane: C_a · L_p(1+p) ≡ p·a mod q, hence C_a · p · L̃_p ≡ p·a, hence C_a ≡ a · L̃_p^{-1} mod p^r. Argument is p-blind.

4. **The factor e_q(c) in front:** when v = 1+ps, c·v = c + cps, so e_q(c·v) = e_q(c) · e_q(cps); the e_q(c) factors out of the sum. At q=3, c=1, this gives e_q(1) (matching R78.4 statement). General c: e_q(c).

**Where q=3-specifics MIGHT appear:**

- (Q1) **The Jacobian factor of p.** The change-of-variables u → s on the principal-unit subgroup gives multiplicity exactly p (group quotient `(Z/p^{r+1})*` / `{u : (1+p)^u ≡ 1}`). This is p-blind, follows from cyclic-group order = p^r.
- (Q2) **Truncation level J_p.** Different p give different J. **NOT a structural obstruction** — same form, different coefficient count.
- (Q3) **L̃_p invertibility.** L̃_p = L_p(1+p)/p mod p^r. From the series, L_p(1+p) = p − p²/2 + p³/3 − ..., so L̃_p = 1 − p/2 + p²/3 − ... ≡ 1 (mod p) for p odd (since p/2 mod p is well-defined for p odd: 2^{-1} mod p exists). Hence L̃_p is a UNIT mod p^r, and L̃_p^{-1} exists. **For p=2: 2^{-1} undefined, structural failure** — confirmed by FHAT doc §4.

**Verdict on T78.4_p:** Generalizes cleanly for p ≥ 3. NO q=3 specifics in the statement itself. J_p replaces J=3, otherwise verbatim substitution.

### Family Theorem 78.5_p (Bijection on support)

> The map `a ↔ C_a` is a bijection on `{a ∈ Z/p^r : a ≡ c (mod p)}` to itself (for fixed c ≡ 1 case: from {a ≡ 1 mod p} to itself).

**Why:** Multiplication by L̃_p^{-1} is multiplication by a unit mod p^r, hence a bijection of (Z/p^r); restricting to {a ≡ 1 mod p} preserves the set since L̃_p ≡ 1 mod p implies L̃_p^{-1} ≡ 1 mod p, so multiplication by L̃_p^{-1} sends 1-mod-p elements to 1-mod-p elements.

**Argument is p-blind.** No q=3 specifics.

### Family Theorem 78.6_p (Saddle-point closed form)

> Define `s*(C_a) = (C_a − 1)/p mod p` (integer division — for C_a ≡ 1 mod p which holds by T78.5_p, this is well-defined).
>
> At r=2 (small r baseline): `|G_p(a)| = p^{(r+1)/2}` exactly (Plancherel + uniform support); phase has additional saddle-multiplicity factor from p simultaneous saddles giving P_a = 0 (analogous to q=3 r=2's e^{iπ/6} from 3 saddles).
>
> At r=3 (medium r): conjectured `G_p(a) = p^{(r+1)/2} · e_q(P_a(s*(C_a)))` exact (analog of R78.6 at q=3 r=3 with J=3).
>
> At r ≥ 4: Hensel correction needed (analog of R79b's full ψ_true vs ψ_lead delocalization).

**Derivation of saddle:** dP_a/ds = pc − C_a · dL_p(1+ps)/ds = pc − C_a · (p − p²s + p³s² − ...) = pc − C_a·p + C_a·p²·s − ... mod q. Setting ≡ 0 mod q: at leading order p·(c−C_a) + C_a·p²·s ≡ 0 mod p², which gives s* ≡ (C_a − c)/(p·C_a) mod p, and using C_a ≡ c mod p (T78.5_p restricted to coset c), s* ≡ (C_a − c)/(p·c) mod p. **For c=1: s* = (C_a − 1)/p mod p.** Argument is p-blind.

**Why r=3 is the saddle-exact case:** At r=3, J_p = J such that J − v_p(J) < r+1 = 4 first fails. For p=3, J=3 (path_B_saddle_point.py J_for_p3 gives J=3 at m=4). For p≥5, J=2 or 3 depending. The saddle integration via stationary-phase in p-adic arithmetic terminates exactly when J ≤ degree of dP_a/ds modulo p^2 contributions. **Will need empirical check** since J_p differs from J=3.

**Where q=3-specifics MIGHT appear:**

- (Q4) **Exact r=3 saddle.** R78.6 states "exact at r=3 with J=3". For p=5 r=3 J=2 (computed above), the saddle prediction MIGHT fail because the polynomial degree of P_a(s) is lower. Need empirical check.
- (Q5) **r=2 phase factor.** At q=3 r=2, R78.6 notes 3 simultaneous saddles all giving P_a(s)=0, yielding Gaussian-integration factor e^{iπ/6}. At general p r=2, we expect p simultaneous saddles giving an analogous factor (some root of unity). Phase prediction needs adjustment, but **magnitude prediction √q is exact** (Plancherel-driven, structurally p-blind).
- (Q6) **Hensel correction at r ≥ 4.** R79b documents this at p=3. Family-level analog should hold structurally (Hensel lifting is p-blind). The specific class-correlated deviation (j=0 anomalous, j≥1 regular) **needs empirical verification at general p** — the j=0 anomaly might be q=3-specific or might generalize.

**Verdict on T78.6_p:** Magnitude generalizes p-blindly (already verified by FHAT doc — `|G_p(a)| = p^{(r+1)/2}` to machine precision at 33 cells). Phase formula at r=3 is **CONJECTURAL at family level** until verified empirically. Hensel correction structure at r ≥ 4 is **OPEN** at family level.

## Summary table

| Theorem | Family-level statement | Status | q=3-specifics in extension? |
|---|---|---|---|
| 78.4_p | F̂_p(p·a) = p·e_q(c)·G_p(a), G_p in closed form via L_p | Rigorous (Cochrane Prop 4 p-blind) | None — J_p replaces J |
| 78.5_p | a ↔ C_a bijection on coset {a ≡ c mod p} | Rigorous (unit multiplication on Z/p^r) | None |
| 78.6_p magnitude | \|G_p(a)\| = p^{(r+1)/2} on support | Verified empirically (FHAT doc, 33 cells) | None |
| 78.6_p phase at r=3 | G_p(a) = p^{(r+1)/2}·e_q(P_a(s*)) exact | Conjectural — needs empirical check at p≥5 since J_p differs | (Q4) — saddle exactness depends on J_p |
| 78.6_p phase at r=2 | Magnitude √q exact; phase needs Gaussian factor | Conjectural — analog of q=3 e^{iπ/6} | (Q5) — root-of-unity factor |
| 78.6_p Hensel at r ≥ 4 | ψ_true delocalizes from ψ_lead; class structure | Open at family level | (Q6) — j=0 anomaly p-specific? |

## Phase 2 verification target

Empirical check cells (locked in pre-reg): `(p, r) ∈ {(3, 2), (3, 3), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2), (11, 3)}`.

Per cell, verify:
- C1: a ↔ C_a is a bijection (exhaustive enumeration; collision count = 0)
- C2: |G_p(a)| = p^{(r+1)/2} for all a in support (FHAT doc already covers, re-verify for completeness)
- C3: At r=2 — magnitude only; phase factor logged
- C4: At r=3 — full saddle prediction G_p(a) = p^{(r+1)/2}·e_q(P_a(s*(C_a))) compared per a

**STOP condition:** If C1 or C2 fails → H_CLOSES with structural obstruction at the failing cell.

## Files

- `PATH2_PRE_REGISTRATION.md` — locked pre-reg
- `PATH2_FAMILY_EXTENSION.md` — this file
- `PATH2_FAMILY_EXTENSION_VERIFICATION.md` (next) + CSV
