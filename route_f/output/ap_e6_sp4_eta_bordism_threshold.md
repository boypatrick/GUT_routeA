# AP-E6 Sp(4) eta/bordism/threshold audit

- Status: `pass_with_exact_gauge_bordism_and_fail_closed_partial_audits`
- Mechanical checks: `39/39`
- Fermion PV/gamma5 controls verified: `true`
- Euclidean regulator complete: `false`
- Omega5 Spin BSp4 computed: `true` = `Z2`
- Displayed gauge character trivial: `true`
- Target eta pair selected: `false`
- Local pure-gauge theta periodicity check: `true`
- Unbroken-group APS determinant ratio computed: `false`
- Gravitational APS phase computed: `false`
- Full heavy threshold eta matched: `false`
- IR k=+2 orientation derived: `true`
- All-scale k=+2 matched: `false`
- Lane closed: `false`
- Physics promotion allowed: `false`

## Partial Euclidean-regulator result

The finite-mode gamma5-reality test and the three PV moment identities pass.
This is not a complete regulator: regulator-field statistics, full
`K(t)/Lambda_a(t)` interpolation, APS boundary/domain data, a finite local
counterterm scheme, and the gauge-scalar-ghost regulator remain unspecified.

## Exact gauge-bordism result

```text
Omega_5^Spin(BSp(2)) = Z2
generator = unit embedded Sp(1) instanton on S4 x S1_R
nu_4=1, nu_5=0, nu_10=0
two Dirac 5s = four Weyl 5s -> total character 0
```

## Target eta result

Gauge bordism and target torsion have different domains.  The current action
does not provide the required target-to-Fredholm mass-family lift.  The two
mod-two invariants are the one-dimensional spin Dirac mod-two index `nu3` and
the two-dimensional Arf invariant `nu2`.  Abstract invertible
sectors/counterterms enumerate all four character pairs:

- `(epsilon_3,epsilon_2)=(0,0)` -> `(G3,G2)=(+1,+1)`
- `(epsilon_3,epsilon_2)=(0,1)` -> `(G3,G2)=(+1,-1)`
- `(epsilon_3,epsilon_2)=(1,0)` -> `(G3,G2)=(-1,+1)`
- `(epsilon_3,epsilon_2)=(1,1)` -> `(G3,G2)=(-1,-1)`

These rows are not microscopic regulator constructions.  Therefore the
microscopic target pair is underdetermined, not guessed.

## Heavy threshold result

```text
5 -> 2_(+1) + 2_(-1) + 1_0
M(+1),M(-1),M(0) = 0,-2yV,-yV at m=-yV
local pure-gauge theta shifts / pi = 2,4 -> periodicity check passes
conditional uniform-flavour charge trace = +2-2+0 = 0
```

The `+2-2` identity is representation-theory bookkeeping under a uniform
exact-flavour hypothesis, not an evaluated determinant.  The unbroken-group
and gravitational APS determinant ratios remain open.  The light B=+1 branch
has the declared infrared orientation k=+2, but no all-scale coefficient has
been matched.

## Checks

- [PASS] `provenance` - all declared sources exist: TeX, bibliography, and verifier are present and nonempty
- [PASS] `provenance` - source hashes are unique sha256 values: hashed=3/3; unique=3
- [PASS] `provenance` - theorem and fail-closed anchors are present: anchors=5/5
- [PASS] `convention` - anti-Hermitian Cartan gives Hermitian charge generator: anti-Hermitian T_H^3 is converted to Hermitian H=-2i T_H^3
- [PASS] `bordism` - low spin coefficients: Omega_q^Spin={0: 'Z', 1: 'Z2', 2: 'Z2', 3: '0', 4: 'Z', 5: '0'}
- [PASS] `bordism` - BSp(2) homology below degree eight: H_4=Z and H_p=0 for 0<p<8, p!=4
- [PASS] `bordism` - unique total-degree-five AHSS entry: nonzero diagonal entries=[(4, 1)]
- [PASS] `bordism` - all outgoing AHSS differentials vanish: d2 and d3 have zero targets; d4:Z2->Z is zero
- [PASS] `bordism` - all incoming AHSS differentials vanish: incoming d2 starts at H_6(BSp2;Z)=0; higher coefficient degree is negative
- [PASS] `bordism` - Omega5 result and extension: Omega_5^Spin(BSp(2))=Z2 with a single filtration quotient
- [PASS] `representation` - SU2 restriction dimensions: restrictions={'4': [2, 1, 1], '5': [2, 2, 1], '10': [3, 2, 2, 1, 1, 1]}
- [PASS] `representation` - instanton indices: 2T/embedded-index={'4': 1, '5': 2, '10': 6}
- [PASS] `representation` - mod-two Witten characters: nu_R={'4': 1, '5': 0, '10': 0}
- [PASS] `representation` - two Dirac five character: Weyl copies=4; character=0
- [PASS] `representation` - fundamental negative control detects generator: one left-Weyl 4 has exponentiated eta phase -1
- [PASS] `regulator` - PV finite-difference moments: moments={0: 0, 1: 0, 2: 0}
- [PASS] `regulator` - finite-matrix gamma5 hermiticity: residual=0.000e+00
- [PASS] `regulator` - one-flavour determinant is real: det=0.941675000000-1.200e-17i
- [PASS] `regulator` - two identical flavours are nonnegative: det^2=0.886751805625
- [PASS] `regulator` - partial controls do not promote full regulator: fermion PV/gamma5 controls pass; four completion classes remain open
- [PASS] `threshold` - positive-charge light tuning: masses={'+1': 0.0, '-1': -2.0, '0': -1.0}
- [PASS] `threshold` - opposite-charge negative control: masses={'+1': 2.0, '-1': 0.0, '0': 1.0}
- [PASS] `threshold` - declared mass signs: signs={'+1': 1, '-1': -1, '0': -1}
- [PASS] `threshold` - local SU2c theta periodicity: Delta theta_c/pi=2
- [PASS] `threshold` - local U1H theta periodicity: Delta theta_H/pi=4
- [PASS] `threshold` - complete five charge trace: 2*(+1)+2*(-1)+1*0=0
- [PASS] `threshold` - light plus heavy mixed coefficient: light=2; heavy=-2; total=0
- [PASS] `threshold` - conditional uniform-flavour trace is mass-sign independent: algebraic controls={'sigma_minus=-1,sigma_zero=-1': 0, 'sigma_minus=-1,sigma_zero=+1': 0, 'sigma_minus=+1,sigma_zero=-1': 0, 'sigma_minus=+1,sigma_zero=+1': 0}; no determinant evaluated
- [PASS] `orientation` - physical infrared k=+2 orientation: Nc*Xq*B=2*(+1)*(+1)=2
- [PASS] `orientation` - opposite branch orientation negative control: Nc*Xq*B=2*(-1)*(+1)=-2
- [PASS] `threshold` - APS determinant-ratio gates remain open: local theta periodicity is not eta matching; unbroken/gravitational APS ratios are open
- [PASS] `target_eta` - all four abstract target torsion characters are enumerated: invertible character pairs=[(-1, -1), (-1, 1), (1, -1), (1, 1)]; microscopic realization not asserted
- [PASS] `target_eta` - target torsion mod-two indices are defined: nu3 is the 1d spin Dirac mod-two index and nu2 is the 2d Arf invariant
- [PASS] `target_eta` - G3 light reference parity: Nc mod 2=0
- [PASS] `target_eta` - G2 unit-flux light reference parity: Nc*Nf*|c1| mod 2=0
- [PASS] `target_eta` - reference controls do not select the microscopic pair: mass-family lift absent; abstract I3/I2 counterterms independently toggle characters
- [PASS] `claim_boundary` - gauge bordism closes without selecting target eta: gauge character and target character have different bordism domains
- [PASS] `claim_boundary` - local theta periodicity does not imply eta matching: local theta shifts pass, but unbroken/gravitational APS determinant ratios are open
- [PASS] `claim_boundary` - no premature lane or portal promotion: full regulator, target eta, APS threshold, naturalness, and strong B=1 gates remain false

## Remaining blockers

- Specify regulator-field statistics, full K(t)/Lambda_a(t) interpolation, APS boundary data, finite local counterterms, and gauge-scalar-ghost regulators.
- Construct one common regulated mass-family lift (g,s)->(A,Sigma,Phi,M,PV) and evaluate both G3 and G2 APS phases.
- Evaluate the unbroken-group plus gravitational APS determinant ratio across the heavy threshold.
- Test whether emergent flavour or explicit UV topological data evades the conditional uniform-flavour charge-trace obstruction.
- Protect the light-branch tuning m=-yV against radiative thresholds.
- Prove the charged two-colour phase and stable B=1 soliton nonperturbatively before any portal construction.

## Source manifest

- `route_f/tex/ap_e6_sp4_eta_bordism_threshold.tex` - `827b4471ceaca7043947308e1ae9eb7c37093338e2bf881e9af615feba42f936` (35879 bytes)
- `route_f/tex/ap_e6_sp4_eta_bordism_threshold.bib` - `9bf27cf0c4540f2bf5b3a72413513d2ff48094653e53b432bad0cd54f06d4687` (3558 bytes)
- `route_f/code/verify_ap_e6_sp4_eta_bordism_threshold.py` - `8cb4e85cb1ef9ae3226813a39083c9f7df2830131c73152deef0c20217e2967a` (28674 bytes)
