# Scattering Bootstrap Skeleton

## External state space

For the first Route-C pass, use one family object as the external-state
container:

```text
R_family = Q + L + u^c + d^c + nu^c + e^c.
```

The Route-A candidate realizes this as

```text
R_family = 16 of Spin(10).
```

Other candidates are tested against the same visible external data.

## Minimal amplitude ansatz

For a candidate heavy mediator `X`, write the four-point tree ansatz as

```text
A_4(1,2,3,4)
  = sum_X g_12X g_34X/(s - M_X^2)
  + sum_Y g_13Y g_24Y/(t - M_Y^2)
  + sum_Z g_14Z g_23Z/(u - M_Z^2)
  + contact.
```

The first audit does not assume a complete Lagrangian.  It checks whether the
pole residues and contact terms can satisfy consistency constraints.

## Constraints

### Factorization

At each physical pole,

```text
Res_{s=M_X^2} A_4 = g_12X g_34X.
```

The residue must be compatible with a positive-norm intermediate state.

### Ward identities

For each gauge boson external leg,

```text
epsilon_mu -> p_mu
```

must annihilate the amplitude after summing all required diagrams.  Failure is
evidence that the assumed mediator/Higgs/source sector is incomplete.

### High-energy behavior

Longitudinal and broken-sector amplitudes must not grow faster than the
assumed UV completion can tolerate.  A useful first flag is any uncancelled
growth of the schematic form

```text
A(s) ~ s^2/M_X^4
```

in a channel that should remain perturbative up to the GUT scale.

### Low-energy matching

Integrating out a heavy GUT mediator gives operators of the form

```text
L_eff ~ (g_GUT^2/M_X^2) qqql.
```

The surviving branches must be checked against proton-decay constraints before
they can be used as action-level support for Route A.

## Candidate comparison logic

This pass should distinguish:

- `SU(5)`: simple and economical, but a family is split as `10 + bar5 + 1`.
- Pati-Salam: natural quark-lepton organization, but not simple.
- `Spin(10)`: one irreducible `16` contains the six visible faces.
- `E6`: one irreducible `27` contains the `Spin(10)` structure but adds extra
  states unless a lifting mechanism is supplied.

