# Route-C P5 Anomaly and Chiral-Consistency Ledger

This ledger compares candidate GUT branches against anomaly and
chiral-spectrum requirements.  It uses the P1 one-family
`SM + nu^c` ledger as the low-energy baseline.

## Baseline

The P1 baseline spectrum is

```text
Q + L + u^c + d^c + nu^c + e^c
```

and the exact checks pass:

```text
Tr Y^2 / Tr T3L^2 = 5/3
all SM-family anomalies pass = True
```

## Candidate Summary

| candidate | anomaly status | single object | exotics/lifting status | Route-C status |
| --- | --- | --- | --- | --- |
| SU(5) | True | False | none for one SM family plus singlet | anomaly-consistent but fails the single irreducible family object criterion |
| Pati-Salam | True | False | none for the standard family pair | anomaly-consistent and physically natural, but fails the simple-group and single-object filters |
| Spin(10) | True | True | none in the standard 16 | standard minimal anomaly-safe single-object survivor under the stated Route-C pre-filter |
| E6 | conditional | True | extra 10 + 1 must be lifted or otherwise accounted for | anomaly-consistent, but not minimal and requires an extra-state lifting audit |
| Superstring-derived completions | conditional | not a universal requirement | compactification-dependent; exotics must be absent, vectorlike, or lifted | UV-completion class with stronger high-energy softness/tower structure, but requires compactification-specific spectrum and exotics audits |

## Detailed Checks

### SU(5)

- category: four-dimensional field-theory GUT candidate
- family realization: `10 + bar5 + 1`
- contains `nu^c`: yes, but as an SU(5) singlet
- minimal single-object filter: False

| check | formula | passes |
| --- | --- | --- |
| SU(5)^3 perturbative anomaly | `A(10)+A(bar5)+A(1)=1-1+0=0` | True |
| global SU(5) anomaly | `pi_4(SU(5))=0` | True |

Missing audits:

- embedding comparison if the single-object criterion is relaxed
- proton and threshold audits for an SU(5)-specific branch

### Pati-Salam

- category: four-dimensional field-theory GUT-like candidate
- family realization: `(4,2,1) + (bar4,1,2)`
- contains `nu^c`: yes, inside the SU(2)_R doublet
- minimal single-object filter: False

| check | formula | passes |
| --- | --- | --- |
| SU(4)_C^3 perturbative anomaly | `2 A(4)+2 A(bar4)=2-2=0` | True |
| SU(2)_L Witten anomaly | `four SU(2)_L doublets from dim(4)=4, even` | True |
| SU(2)_R Witten anomaly | `four SU(2)_R doublets from dim(bar4)=4, even` | True |

Missing audits:

- explicit product-group bootstrap branch if simplicity is relaxed
- broken-sector Ward and proton audits for Pati-Salam leptoquark currents

### Spin(10)

- category: four-dimensional field-theory GUT candidate
- family realization: `16 half-spinor`
- contains `nu^c`: yes, inside the half-spinor
- minimal single-object filter: True

| check | formula | passes |
| --- | --- | --- |
| Spin(10)^3 perturbative anomaly | `D5 has no independent cubic gauge-anomaly invariant; 16 is anomaly-safe` | True |
| Spin(10) global anomaly | `pi_4(Spin(10))=0 in the standard four-dimensional check` | True |

Missing audits:

- explicit broken-generator matrices for Ward checks
- Higgs/source sector for broken-vector completion
- full low-energy proton and threshold matching

### E6

- category: four-dimensional field-theory GUT candidate
- family realization: `27 -> 16 + 10 + 1 under Spin(10)`
- contains `nu^c`: yes, through the embedded Spin(10) 16
- minimal single-object filter: False

| check | formula | passes |
| --- | --- | --- |
| E6^3 perturbative anomaly | `E6 has no independent cubic gauge-anomaly invariant for the 27` | True |
| extra-state anomaly pairing after restriction | `27 contains 16 plus vectorlike 10 and singlet under Spin(10)` | conditional on a lifting/vectorlike audit |

Missing audits:

- mass/lifting mechanism for 10 + 1 extra states
- check that exotics do not reintroduce low-energy chiral matter
- E6-specific threshold and proton audits

### Superstring-derived completions

- category: UV-completion class, not a single four-dimensional group
- family realization: `model-dependent; may realize SU(5), Pati-Salam, Spin(10), E6, or other quiver/brane/heterotic branches`
- contains `nu^c`: compactification-dependent
- minimal single-object filter: not directly comparable

| check | formula | passes |
| --- | --- | --- |
| higher-dimensional anomaly cancellation | `Green-Schwarz/tadpole/modular consistency as appropriate to the string construction` | construction-dependent |
| four-dimensional anomaly cancellation | `massless spectrum plus any 4D Green-Schwarz terms must cancel anomalies` | construction-dependent |
| worldsheet or compactification consistency | `modular invariance, tadpole cancellation, K-theory/Freed-Witten-type constraints where applicable` | construction-dependent |

Missing audits:

- explicit compactification and massless spectrum
- modular/tadpole/Green-Schwarz consistency data
- chiral index and exotics-lifting audit
- matching to the Route-A representation/index skeleton

## P5 Boundary

Passing anomaly checks does not by itself select a GUT.  `SU(5)`,
Pati-Salam, `Spin(10)`, and `E6` all have anomaly-consistent forms,
but they satisfy different structural filters.  Superstring-derived
completions are a UV-completion class rather than one fixed
four-dimensional candidate; they require compactification-specific
modular/tadpole/Green-Schwarz, spectrum, and exotics-lifting audits.
