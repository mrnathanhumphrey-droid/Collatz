# ADELIC_J — Saloff-Coste 2001 (random walks + invariant diffusions on groups)

**Source:** C:/tmp/adelic/Saloff_Coste_Notices_Probability_Groups.txt. Notices AMS survey.

## Verbatim Theorem 2 (Varopoulos, p. ~4)

> "Theorem 2. Assume that there exists a positive constant c such that V(n) ≥ c n^d for all n. Then there are positive constants C_1 and c_1 such that
>    φ(n) ≤ C_1 n^{-d/2} and I(n) ≥ c_1 n^{1 - 1/d} for all n."

(V(n) = volume of ball of radius n in the Cayley graph; φ(n) = probability random walk returns to origin at time n; I(n) = isoperimetric profile.)

## Verbatim Theorem 3 (Gromov+Varopoulos, p. ~5)

> "Theorem 3. For a finitely generated group G, the following are equivalent properties: (1) V(n) ≍ n^d; (2) I(n) ≍ n^{1 - 1/d}; (3) φ(n) ≍ n^{-d/2}; (4) G contains a nilpotent subgroup N of finite index, and d = Σ_i i · r_i, where r_i is the torsion-free rank of the abelian group N_i / N_{i+1}, and (N_i) is the lower central series of N."

## Verbatim Theorem 4 (Varopoulos / Coulhon-Saloff-Coste, p. ~5)

> "Theorem 4. Fix α ∈ [0, 1]. Assume that there exists a positive constant c such that log V(n) ≥ c n^α for all n. Then there are positive constants c_1 and c_2 such that
>    log φ(n) ≤ −c_1 n^{α/(α+2)} and I(n) ≥ c_2 n / [log n]^{1/α} for all n."

## Hypotheses isolated

- **h1 (GROUP):** G a finitely generated group with a symmetric generating set S; Cayley graph (G, S).
- **h2 (RANDOM WALK):** Simple symmetric random walk on Cayley graph; equivalently SYMMETRIC probability measure μ on G (μ(g) = μ(g^{-1})).
- **h3 (REVERSIBILITY):** The walk is *reversible* w.r.t. counting measure — this is the standard setup for the heat-kernel methods.
- **h4 (VOLUME GROWTH):** V(n) (number of group elements with word-length ≤ n) satisfies V(n) ≥ c n^d (polynomial) or log V(n) ≥ c n^α (exponential).
- **CONCLUSION:** φ(n) ≤ C n^{-d/2} (return probability decay); isoperimetric profile.

## Hypothesis × input check

| Hyp | Syracuse |
|---|---|
| h1 (finitely generated group) | (Z/3^n)* IS a finite group, hence trivially finitely generated. ✓ But Saloff-Coste's setting is for *infinite* finitely generated groups; the asymptotics are at n → ∞ in the *walk-step* variable, on a fixed group. We'd want to apply it to the *inverse limit* lim (Z/3^n)* = ℤ_3* — which IS infinite and finitely generated (topologically by 2). |
| h2 (symmetric μ) | FAILED — Tao's recursion is fundamentally NOT symmetric. The step r ↦ (3r+1)/2^v has no inverse in the natural sense; the *inverse Syracuse map* exists but is multi-valued (2-to-1 fibration). Saloff-Coste's framework explicitly requires symmetric μ. |
| h3 (reversibility) | FAILED for the same reason — Tao K_k is the *forward* transition matrix; its reversibility ratio K_k(r → r') · π(r) vs K_k(r' → r) · π(r') is NOT 1 in general for Tao. |
| h4 (volume growth) | (Z/3^n)* has size 2·3^{n-1} which grows exponentially in n; ℤ_3* topologically has continuous Haar mass. Not directly the V(n) of Saloff-Coste's framework. |

## Disposition for J

**NO_FIT (categorical, REVERSIBILITY).**

This is the same disposition as in FG_DISPOSITION (Furstenberg-Guivarc'h probe): Saloff-Coste's heat-kernel framework requires SYMMETRIC random walks, and the Tao recursion is fundamentally non-reversible. The forward Syracuse map has no symmetric counterpart in the natural category.

Already pre-classified as low prior (3%) per ADELIC_PRE_REGISTRATION.md.

**Adelic factorization tag:** **GLOBAL_BUT_PLACE_BLIND** — the random-walk framework is on an abstract Cayley graph; no place-by-place factorization.

## Mode H circular fingerprint

None: h2/h3 fail categorically.
