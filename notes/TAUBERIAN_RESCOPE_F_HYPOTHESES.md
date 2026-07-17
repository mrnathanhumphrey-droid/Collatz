# TAUBERIAN_RESCOPE_F_HYPOTHESES (Singha Roy 2511.15928 — LSD for L-functions)

**Source PDF:** `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/arxiv_2511.15928_Landau_Selberg_Delange_L_functions.pdf`
**Extracted text:** `C:/Collatz/tauberian_extract/F.txt`
**Mode:** E — verbatim from PDF.

---

## Theorem 1.1 — verbatim (PDF lines 174-243)

> **Theorem 1.1.** Assume that {a_n}_n has property P(ν, {α_χ}_χ; c_0, Ω), and that for all x > 1, we have
>
>   Σ_{x < n ≤ 2x} |a_n| ≤ κ x^{1/ν},  where κ ≥ 2 is independent of x.    (1.5)
>
> Fix any K_0 > 0. The following hold uniformly in all x ≥ q ≥ e^{4 + 5/3ν}, h ∈ (0, x/2], N ∈ ℤ_{≥0}, and in {α_χ}_χ ⊂ ℂ with max{|α_{χ_0}|, |α_{χ_e}|} ≤ K_0; the implied constants depend only on c_0, ν and K_0.
>
> (1) If the exceptional zero η_e exists and satisfies 1 − c_0/(10 λ_q log q) < η_e < 1 − 3ν/log x, then
>   | Σ_{n ≤ x} a_n − (x^{1/ν} / (log x)^{1 − α_{χ_0}}) Σ_{j=0}^N (μ_j (log x)^{−j} / Γ(α_{χ_0} − j)) |
>   ≪ Σ_{x < n ≤ x+h} |a_n| + κ · x^{1 + 1/ν} (log x) / (T h)
>     + (2 λ_q log q)^{λ_q + 2 K_0} x^{1/ν} { Ω(T) (log(eνT))^{1 + λ_q} / T + Ω(1) N! [71 (1 + ν) (1 − η_e)^{−1}]^{N + K_0 + 1} (log x)^{1 − |Re(α_{χ_0})|} · min{x/h, (log x)^{N+1}} }
>   where T := (qν)^{−1/2} exp( 0.5 √(q log²(qν)) + c_0 log x / (ν λ_q) ).
>
> (2) If η_e does not exist or satisfies η_e ≤ 1 − c_0/(10 λ_q log q), then for q < x^{c_0/(80 ν λ_q)}, we have [similar formula (1.7)] …

Property P(ν, {α_χ}_χ; c_0, Ω) is defined earlier (verbatim from PDF lines 144-170):
> a sequence {a_n}_n ⊂ ℂ has the property if its Dirichlet series F(s) := Σ a_n n^{−s} has, in some explicit Dirichlet-character decomposition over (ℤ/qℤ)^* with parameters (α_χ)_{χ mod q}, a log-derivative behavior matching ζ(sν)^{α_{χ_0}} · ∏_{χ ≠ χ_0} L(sν, χ)^{α_χ}, with Ω(t) controlling vertical-strip growth, and c_0 controlling the zero-free region.

---

## Hypotheses extracted (load-bearing list)

| # | Hypothesis | Source |
|---|---|---|
| h_1 | **{a_n}_n ⊂ ℂ has property P(ν, {α_χ}_χ; c_0, Ω)** — Dirichlet-character decomposition with α_χ parameters and prescribed zero-free / vertical-strip data. | line 174 |
| h_2 | **Growth bound** Σ_{x < n ≤ 2x} |a_n| ≤ κ x^{1/ν} for all x > 1, κ ≥ 2 independent of x. | line 178 |
| h_3 | **q (Dirichlet modulus)** satisfies q ≥ e^{4 + 5/3ν}, and the conclusion holds uniformly for q in a specific range depending on whether the exceptional zero η_e exists. | line 179 |
| h_4 | **N ∈ ℤ_{≥0}** (asymptotic-expansion truncation order) is free; constants depend on c_0, ν, K_0. | line 179 |

Conclusion: asymptotic expansion of Σ_{n ≤ x} a_n in powers of (log x)^{-j}.

---

## Notational mapping for our use case

Theorem F.1.1 is built for Dirichlet series in *arithmetic progressions* with explicit Dirichlet-character data. The c=7/45 closure is not in this category — there's no Dirichlet-modulus q structure on the ε_k sequence. The 3-adic / 2-adic structure of Syracuse is **not a Dirichlet-character structure** in the sense of property P(ν, {α_χ}_χ; c_0, Ω) — that would require α_χ assigned to Dirichlet characters χ mod q for some integer q with explicit L-function decomposition.

This is anticipated to fail h_1 — see HYPOTHESIS_CHECK.

---

## End of F HYPOTHESES extraction.
