# L3 DEFINITIONS — the frozen object (Phase 0 contract)

**All later phases reference THIS page. Any change = new STATE entry.**

> ⚠️ **NOT FROZEN — G0 CAUGHT A DRIFT (2026-07-16, see `result_G0_phase0.md`).** §1–§6 stand. **§7/§9 are REFUTED as written:** the uniform-lift refinement operator L_k does NOT carry the spectral gap — its gain γ_k = 1/√3 for EVERY q (gapless, first-moment), and σ_max(L_k) = 0.61–0.67 > 1/√3 for all q. The gap r_q is second-moment (build_M / R42), where R46's `cB_k` rate actually measured it. §7/§9 must be re-authored on the second-moment object (Nathan to adjudicate). Kept below only as the first-moment scaffold + the record of the caught drift.

Fix an odd prime `q ≥ 3`. Let `d := ord_q(2)`. Indices `k ≥ 1`. Weight truncation `V = 64`.

---

### §1 State space
`U_k := { r ∈ Z/q^k : q ∤ r }`, the units-adjacent (coprime-to-q) residues. `|U_k| = (q−1)q^{k−1}`.
Real Hilbert space `ℓ²(U_k)` with `⟨f,g⟩ = Σ_{r∈U_k} f(r)g(r)`, `‖f‖² = Σ_r f(r)²`.

### §2 Transfer operator `K_k`
Branches `v ≥ 1`: `r ↦ (q r + 1)·2^{−v} (mod q^k)` with weight `p_v = 2^{−v}/Z_V`, `Z_V = Σ_{v=1}^{V} 2^{−v} = 1 − 2^{−V}`.
`(K_k)_{r,r'} = Σ_{ v=1..V : (q r+1)2^{−v} ≡ r' (q^k) } p_v`.  Row-stochastic. **Measures act on the LEFT:** `(μ K_k)(r') = Σ_r μ(r)(K_k)_{r,r'}`. As a column operator this is `K_k^{⊤}`.

### §3 Stationary measure
`π_k ∈ ℓ²(U_k)`: `π_k K_k = π_k`, `Σ π_k = 1`, `π_k > 0` (Perron–Frobenius; K_k irreducible).

### §4 Reduction and uniform lift
Reduction `ρ_k : U_{k+1} ↠ U_k`, `ρ_k(r') = r' mod q^k` (well-defined since `q ∤ r' ⇒ q ∤ (r' mod q^k)`; exactly `q` preimages per point).
Lift `lift_k : ℓ²(U_k) → ℓ²(U_{k+1})`, `(lift_k u)(r') = u(ρ_k(r'))/q` (uniform spread over the fiber).

### §5 Deviation subspace and projection
`W_{k+1} := { w ∈ ℓ²(U_{k+1}) : Σ_{r' ∈ ρ_k^{−1}(r)} w(r') = 0  ∀ r ∈ U_k }` (fiberwise mean-zero).
`dim W_{k+1} = (q−1)·|U_k| = (q−1)² q^{k−1}`.
Orthogonal projection `P_{W_{k+1}} : ℓ²(U_{k+1}) → W_{k+1}`,  `(P_W f)(r') = f(r') − (1/q) Σ_{s' ∈ ρ_k^{−1}(ρ_k(r'))} f(s')`.
(Convention: `W_1 :=` fiberwise-mean-zero part of `ℓ²(U_1)` over the single fiber `U_1`.)

### §6 The deviation vector
`d_k := π_k − lift_{k−1}(π_{k−1}) ∈ W_k`,  i.e. `d_k(r) = π_k(r) − π_{k−1}(r mod q^{k−1})/q`.
Norm identity (R74, to reprove Phase 1): `‖d_k‖² = ‖π_k‖² − (1/q)‖π_{k−1}‖²`.

### §7 THE OPERATOR `L_k` (the object of L3)
`L_k : W_k → W_{k+1}`,   `L_k := P_{W_{k+1}} ∘ ( · K_{k+1}) ∘ lift_k`,
that is  `L_k w = P_{W_{k+1}}[ lift_k(w) · K_{k+1} ]  = P_{W_{k+1}} K_{k+1}^{⊤} lift_k(w)`.
**`σ_max(L_k) := ‖L_k‖ =` largest singular value** of `L_k` w.r.t. the `ℓ²` norms of §1 (a max over ALL `w ∈ W_k`, `‖w‖=1` — NOT the gain on one vector).
Also define the **stationary gain** `γ_k := ‖L_k d_k‖ / ‖d_k‖` (the gain on the deviation `d_k` specifically). A priori `γ_k ≤ σ_max(L_k)`; whether they coincide is a Phase-0/2a question, not an assumption.

### §8 Frozen facts (from R45/R46 — to be upgraded to theorems in Phase 1–2)
- **(REFINE)** `L_k d_k = d_{k+1}` exactly (verified 1e-16, R46). ⇒ `γ_k² = ‖d_{k+1}‖²/‖d_k‖²`.
- **(LEM-FORGET)** `(q r+1)2^{−v} mod q^k` is independent of the `q^{k−1}`-digit of `r` (since `q·q^{k−1} = q^k ≡ 0`) ⇒ rows of `K_k` are constant on fibers ⇒ `w K_k = 0` for `w ∈ W_k` (self-block `P_W K P_W = 0`).
- **(ONE-STEP)** `π_{k+1} = lift_k(π_k)·K_{k+1}` exactly (R46).
- **r_q relation:** `r_q = 3·γ_k²` in the `k→∞` limit (`γ_∞² = r_q/3`; measured `r_3=1, r_5≈0.62, r_7≈0.39`). Whether `σ_max = γ` (dominant mode = deviation) is the open Phase-0 fork.

### §9 THE CLAIM (L3)
`σ_max(L_k) < 1/√3` for `d ≥ 3`, uniformly in `k` (and in `q` at fixed `d≥3`); with **equality `σ_max = 1/√3` at `d = 2`** (i.e. `q = 3`) — the isometry/marginal boundary. Then `r_q = 3σ_∞² < 1` ⇒ spectral gap ⇒ `‖π_k‖² ~ C_q 3^{−k}` with a genuine gap for `q ≥ 5`, marginal (linear) at `q = 3`.

---
**GATE G0 (Phase 0):** code §1–§7 literally (fresh, importing nothing from prior probes) and compute BOTH `σ_max(L_k)` (full SVD) and `γ_k` at `q = 3,5,7`, per `k`. Pass bar: (i) `γ_k` reproduces R46's √(r_q/3) column to 1e-10; (ii) `q=3` value = `1/√3` to 1e-10; (iii) REFINE `L_k d_k = d_{k+1}` to 1e-12. **Report whether `σ_max = γ` or `σ_max > γ`** — this decides whether §9 bounds the operator norm or the dominant mode (Phase 2a fork previewed).
