# P54 common scalar spectrum and interaction-defined channels

Date: 2026-09-10. Action and background unchanged. This is a tree scalar-exchange audit, not a fitted particle spectrum or loop-complete scattering prediction.

Validation: **153/153**. No cached Hessian was recomputed.

## Physical vertices, not arbitrary observation projectors

The HdaggerH vertex is the singlet trace of the cubic tensor on the exact tree light doublet. Seven named fermion bilinears use actual Clifford/Yukawa tensors; h_raw and f_raw are factored out, not fitted. Massive vector-pair vertices are derivatives of the actual kinetic-term mass matrix, averaged over complete degenerate vector eigenspaces.

There are 33 real current slots (including selection-rule zeros), source rank 19, and a 21-dimensional channel-generated invariant subspace inside all 295 physical real scalar coordinates. The 33 gauge directions are excluded; all five physical tree zero modes remain.

For physical scalar H and interaction vertices J, $$F_E(t)=J^T(H+tI)^{-1}J=\sum_\lambda R_\lambda/(t+\lambda),\quad R_\lambda=J^T\Pi_\lambda J\succeq0.$$

Residues are aggregated over exact degenerate eigenspaces, not assigned to arbitrary eigenvectors. The JSON contains every pole cluster and its full cross-channel residue matrix. Diagonal entries are squared vertex strengths, not branching fractions.

## Singlet-channel poles and squared couplings

HH weights are in omega^2; NN weights have the family spurion removed and are dimensionless. Other channels can see additional poles listed in JSON.

| m^2 / omega^2 | HH weight | NN weight per unit f_raw |
|---:|---:|---:|
| 0.01352047 | 2.35348682e-07 | 15.9472386 |
| 0.0997156902 | 0.000209220862 | 0.0492598746 |
| 2.82916586 | 0.0754058603 | 0.00350147563 |

The physical PQ mode lies in the S phase, which has no direct Yukawa tensor. Its residues vanish in ALL declared linear channels, including NN and HdaggerH. An initial nonzero-NN-PQ hypothesis was rejected by the actual vertices, not repaired by fitting. This is a physical mode invisible to this current set, not a deleted state, a gauge artifact, a globally decoupled axion, or a claim of zero QCD axion mass. The massive inverse spectral moment independently reproduces the previous nonzero CHN coefficient, including its common complex phase.

The symmetry-enforced zero NN axion residue has numerical residual `2.08e-32`; this is not a predicted tiny nonzero coupling. The massive moment gives `omega CHN / f_raw = -0.0848692011507 i`, with real part consistent with zero.

## Exact nested elimination must also transport vertices and contact terms

For retained r and eliminated h, $$K'=K_{rr}-K_{rh}K_{hh}^{-1}K_{hr},\quad J'=J_r-K_{rh}K_{hh}^{-1}J_h,\quad C'=C+J_h^TK_{hh}^{-1}J_h.$$

The full response is $$F=C+J^TK^{-1}J=C'+J'^TK'^{-1}J'.$$

The actual nested dimensions are `[9, 231, 55]`. The final chart is a factorization device, not an asserted SM-only EFT. Direct, staged, and opposite-order elimination agree at spacelike and complex near-pole momenta, and under noncanonical coordinate rescaling.

| p_E^2 / omega^2 | full vs staged relative residual | error if induced contact omitted |
|---:|---:|---:|
| 0.001+0j | 1.25e-15 | 0.366 |
| 0.03+0j | 9.05e-16 | 0.965 |
| 0.2+0j | 3.79e-16 | 0.984 |
| -0.02+0.003j | 2.81e-16 | 0.992 |
| -0.265+0.001j | 2.65e-15 | 1 |

The induced contact is nonlocal before a derivative expansion and can carry poles. Omitting it loses 100% of the HH-to-HH scalar-exchange response at every tested point. The matrix-norm errors above use the declared unit-spurion current normalization, not fitted cross sections. Keeping K alone is not an observable-preserving projection. Feshbach-Schur transitivity is existing mathematics, not evidence of a fourth spatial direction.

## Scope and next gates

- Preserve the action and symbolic UV family matrices. The selected currents do not exhaust all physical channels; darkness is always relative to the declared current set.
- Use these source/contact identities in the common finite matching workflow. Full vector/Goldstone/ghost/scalar terms, running, widths, and a physical flavor fit remain open.
- Extra spatial dimensions are comparison-only. A fixed interval spectrum would obey m_n^2=M_5^2+(n*pi/L)^2 and therefore (m_2^2-m_0^2)/(m_1^2-m_0^2)=4. No extra-dimensional parameters or particle assignments are fitted here.

## Reproduction

Use the pinned `code/requirements_p54_matching.txt` environment and the existing full-doublet/kinetic Hessian cache. This consumer is read-only on the cache. The exact keys and file hashes are recorded in JSON; an absent or incompatible cache fails explicitly.

```sh
python3 route_f/code/verify_p54_spectral_channels.py --cache-dir /path/to/tmp/p54_full_doublet_cw
```

Derivation: [TeX](../tex/p54_spectral_channels_nested_elimination.tex), [PDF](pdf/p54_spectral_channels_nested_elimination.pdf).

References: [Feshbach-Schur map](https://arxiv.org/abs/2105.02058), [nonlocal covariant EFT matching](https://arxiv.org/abs/1604.01019).
