# Audit 4a Source-Sector / Heavy-Spectrum Schema

Audit 4a has been opened as a schema and decision card, not as a completed UV source-sector model.

## Digest

- card sha256: `be735ff9f65a62f53a7cd1cbf504b4113d0f65dcaec553a4f6bb55e410471fbc`
- audit0 sha256: `a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc`
- source artifacts: 7/7 present

## Near-Term Track

- Default: constrained 54/210 field-theory source-sector interface.
- Fallback: CMSGUT-like 210+126+10 spectrum import if the constrained branch cannot export a complete heavy spectrum.
- Optional parallel: Route-D global/string geometry, but not as a replacement for a heavy-spectrum artifact.

## Required Heavy-Spectrum Fields

state_id, sector, spin10_representation, ps_representation, sm_representation, multiplicity, mass_GeV, mass_expression, threshold_beta_vector_b1_b2_b3, d5_role, source_artifact, status

## Blockers

- No complete heavy_spectrum.json exists yet.
- M_star/v_R/source-scale convention remains unset for publication.
- Doublet-triplet splitting and colored-triplet inverse propagator are not exported in final Audit-2-ready form.
- The constrained 54/210 branch has local conormal evidence but not a complete microscopic UV completion.
- Route-D placement is local and cannot replace a global compactification spectrum.

## Boundary

This file does not claim threshold closure, proton safety, a complete UV action, or a global compactification. It only defines the source/heavy-spectrum output contract required before Audit 3 and Audit 2 can be run as final audits.
