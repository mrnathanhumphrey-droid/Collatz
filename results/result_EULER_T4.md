# RESULT — EULER-T4: NO structured product form for S_∞. Space closed by exclusion — A EXCLUDED, B VACUOUS, C RE-RUN⊕EXCLUDED, all backstopped by MAHLER. (2026-08-28)

**Probe:** `probes/probe_t4.py` (Wilson-penned kill-first spec). Machinery reused verbatim (`exact_Re`
ladder via the certified EPS json; `build_M` from `probe_25_transfer_operator_Aprime.py`). Guardrails
G1–G5 confirmed verbatim (see `notes/CORPUS_LEDGER.md` and the two extraction passes).

**Question (sharpened, kill-first):** a trivial `S_∞=(7/15)·Π(1+δ_i)` with `δ_i` back-solved is VACUOUS
(the R9-D/R10 tautology trap). T4 excludes a **structured** product — one whose factors are
*independently* determined (periodic orbits / places / cyclotomic data) and carry information the tail
**sum** `Σ Λ_i` does not.

## Anchors (exact)
`S_1=2/3`, `S_2=10/21`; `Λ_1=−2/21`; `T_0=1/3`; tail `= S_∞/2 − T_0 = Σ_{i≥1}Λ_i`. Dead leading mode
`7/15` (`⟺ Σ_{r≥1}Λ_r=−1/10`, char-ledger R10); true `S_∞≈0.475 > 7/15`.

## NULL A — dynamical / Ruelle-zeta Euler product → **EXCLUDED**
CLAIM: `S_∞ = ` residue at the leading `1/3` pole of `det(1−zM)`, Euler product `Π_γ(1−z^{n_γ}A_γ)^{−1}`
over primitive orbits. Three independent kills:
1. **Framing NO_FIT (pre-disposed):** `S_n = Σ|μ̂_n(ξ)|²` is a **Plancherel sum, not a dynamical trace
   `tr(Lⁿ)`** — `FAURE_IJ_HYPOTHESES:16,20` "Syracuse: not in this category." No orbit product exists to
   take a residue of.
2. **Leading residue = the dead value.** On the certified operator `M` (`build_M`), the generating
   function `⟨1|(1−zM)^{−1}|v₀⟩ = Σ‖π_k‖²z^k` has `λ_1 → 1/3` (Perron; measured `L=1: 0.555556`,
   `L=2: 0.346827`). Its leading residue = the leading `3^{−k}` mode = **`7/15`, not the tail `0.475`.**
3. **Subdominant pole dissolves.** `|λ_2|/λ_1` climbs `0.6000 (L=1) → 0.9916 (L=2) → 1` — the q=3
   exceptional point (`biov→0`, cond# `2.5e14→2.4e17`, `L=4` witness inaccessible, G5). The discrete
   subdominant factor does **not survive `L→∞`**; it becomes the branch cut / continuous spectrum
   (README:196, INTERLEVEL:57–59, D1_T_M:54). No discrete Euler factor in the limit.

## NULL B — `⟨2⟩`-cascade cyclotomic product → **VACUOUS** (the sharpest, most concrete kill)
The `Λ_i` **denominators** factor *exactly* into the cyclotomic `⟨2⟩`-cofactor primes — e.g.
`den(Λ_4) = 7²·19²·73²·163·2593·135433·262657²·87211²·71119·97685839·272010961` (all cofactors of
`2^{2·3^{i−1}}−1`). That factorization is **pure number theory — vacuous.** The value lives in the
**numerators**, and the numerators carry **FOREIGN primes with no cascade structure**:
```
   Λ_1 num = −2                      = −2
   Λ_2 num = −1490                   = −2·5·149
   Λ_3 num =  2849957897648150       = 2·5²·163·883·396022747
   Λ_4 num =  3479…371450            = 2·3²·5²·31249·197209152271·4633083256267931·270785938212059525530169
```
Foreign numerator primes (NOT in the cofactor set): `{2, 5, 149, 883, 31249, 396022747, 197209152271,
4633083256267931, 270785938212059525530169}`. A product over the cascade primes `{7,19,73,163,…}`
could **never** reproduce `396022747` or `270785938212059525530169`. So the tail is **not** a product
over the cascade: denominators-only factoring is a rewriting, and the numerators — where the value is —
carry no independent cascade product. **VACUOUS.** (And it is not the disposed cyclotomic-7, G2.)

## NULL C — infinite Mahler product of `μ̂` → **RE-RUN ⊕ EXCLUDED**
1. **Premise fails (G4).** The FE `μ̂_{n+1}(ξ)=Σ_v 2^{−v}e(·)μ̂_n(ξ2^{−v})` is a **weighted sum** over `v`,
   not a single-argument rescaling `f(ξ)=a(ξ)f(qξ)`. Iterating two steps gives a **double sum over
   valuation paths** `Σ_{v,v'}2^{−(v+v')}e·e'·μ̂_{n−1}(ξ2^{−(v+v')})` — the **renewal/path-sum** (G3),
   not a product. **RE-RUN** of the renewal.
2. **Plancherel square breaks any product.** Even granting `μ̂=Π f_k`, `S_n=Σ_ξ|μ̂_n|²`, and
   `|Σ_v a_v|² = Σ_v|a_v|² + Σ_{v≠v'}a_v\bar a_{v'}` — the off-diagonal is irreducible. Demo (`ξ=0.37`):
   `|Σa_v|² = 0.839259 = diag 0.333333 + off-diag 0.505926` — off-diagonal is **60% of the total**. That
   off-diagonal *is* the non-free/monotone part carrying the value beyond the diagonal
   (free-prob ledger, `|φ(X̃ⱼX̃ₖX̃ⱼ)|=0.1078`). `Σ→Π` does not survive `|·|²`. **EXCLUDED.**

## MAHLER backstop (any survivor)
`ζ(2)=π²/6` closes because product **+ a closing functional equation**. MAHLER (proven) = the closing FE
does **not** exist at finite order (infinite Mahler depth = branch cut = continuous spectrum, G5). Any
product whose value-extraction needs a closing FE ⟺ finite Mahler order ⟺ **contradicts MAHLER**. Even a
formal product is non-closing (infinitely many independent factors) = the wall.

## Verdict
**The space where a STRUCTURED product for `S_∞` could hide is CLOSED.** A EXCLUDED (no orbit product;
leading residue = dead `7/15`; subdominant → continuous spectrum). B VACUOUS (denominators trivial;
numerators carry foreign, unstructured primes). C RE-RUN⊕EXCLUDED (FE-iteration = renewal sum; Plancherel
`|·|²` breaks the product). **No SURVIVED-A-KILL.** Consistent with every framework in the corpus: the
value lives in the tail, and the tail admits **no** independent product representation — the Euler move
that gives `π²/6` requires a closing functional equation, which MAHLER forbids here.

**Not at stake:** the tail SUM `Σ Λ_i` (exact, certified) and the leading-mode identities are untouched;
T4 only excludes a *product* representation of the value.
