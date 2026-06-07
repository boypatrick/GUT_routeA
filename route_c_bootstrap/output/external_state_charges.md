# Route-C P1 External State and Charge Ledger

This is a pre-amplitude ledger for one left-handed Standard Model
family plus `nu^c`.  It fixes the charge conventions used by later
pole-factorization, Ward-identity, anomaly, and proton-matching
checks.

## State Table

| state | SU(3)_C | SU(2)_L | Y | B-L | multiplicity | convention |
| --- | --- | --- | --- | --- | ---: | --- |
| Q | `3` | `2` | 1/6 | 1/3 | 6 | left-handed Weyl |
| L | `1` | `2` | -1/2 | -1 | 2 | left-handed Weyl |
| u^c | `bar3` | `1` | -2/3 | -1/3 | 3 | left-handed conjugate Weyl |
| d^c | `bar3` | `1` | 1/3 | -1/3 | 3 | left-handed conjugate Weyl |
| nu^c | `1` | `1` | 0 | 1 | 1 | left-handed conjugate Weyl |
| e^c | `1` | `1` | 1 | 1 | 1 | left-handed conjugate Weyl |

All entries use left-handed Weyl conventions.  The fields
`u^c`, `d^c`, `nu^c`, and `e^c` are left-handed conjugate fields,
so their listed charges are the charges of those conjugate Weyl
fields.

## Hypercharge Normalization

For one full family plus `nu^c`,

```text
Tr Y^2       = 10/3
Tr T3L^2     = 2
Tr Y^2 / Tr T3L^2 = 5/3
```

The normalization check passes exactly:

```text
True
```

## Anomaly Checks

| check | exact value | passes |
| --- | ---: | --- |
| `SU3_cubed` | 0 | True |
| `SU3_squared_U1Y` | 0 | True |
| `SU2_squared_U1Y` | 0 | True |
| `grav_squared_U1Y` | 0 | True |
| `U1Y_cubed` | 0 | True |
| `SU3_squared_BminusL` | 0 | True |
| `SU2_squared_BminusL` | 0 | True |
| `grav_squared_BminusL` | 0 | True |
| `BminusL_cubed` | 0 | True |
| `U1Y_squared_BminusL` | 0 | True |
| `U1Y_BminusL_squared` | 0 | True |

The Witten `SU(2)` condition also passes:

```text
number of left-handed SU(2) doublets = 4
even doublet count = True
```

## Current Sectors

The following current ledger is schematic.  It separates currents
that preserve the SM face from currents that leave the SM face and
therefore require broken-GUT mediator and Goldstone/Higgs-sector
data in later Route-C stages.

| current sector | representative current | classification | later use |
| --- | --- | --- | --- |
| SU(3)_C current | `bar(psi) gamma_mu T_C psi` | SM-preserving | unbroken gauge Ward identities |
| SU(2)_L current | `bar(psi) gamma_mu T_L psi` | SM-preserving | unbroken gauge Ward identities |
| U(1)_Y current | `sum_i Y_i bar(psi_i) gamma_mu psi_i` | SM-preserving | hypercharge Ward identities and normalization |
| Pati-Salam leptoquark currents | `Q <-> L, u^c <-> nu^c, d^c <-> e^c` | broken-GUT schematic | massive mediator pole and proton-decay matching audit |
| SU(2)_R charged currents | `u^c <-> d^c, nu^c <-> e^c` | broken-GUT schematic after Y selection | broken-sector Ward and Goldstone completion audit |
| Spin(10) off-face currents | `16 weights connected by broken D5 roots` | broken-GUT schematic | candidate mediator ledger before explicit D5 generator matrices |

## P1 Boundary

This file verifies charge and anomaly bookkeeping.  It does not yet
construct four-point amplitudes, residues, or Ward-identity
cancellations.  Those begin in P2 and P4.
