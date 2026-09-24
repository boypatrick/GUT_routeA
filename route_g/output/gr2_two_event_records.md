# G-R2: one initial record and one final readout

Status: bounded-done. Checks: 139/139.
All quantities are dimensionless; no clock slowdown is fitted.

## Exact two-event test

Rows: initial outcome 0,1. Columns: terminal outcome +,-.

| Initial / terminal | + | - | Initial marginal |
|---|---:|---:|---:|
| 0 | 1/9 | 2/9 | 1/3 |
| 1 | 1/9 | 5/9 | 2/3 |

P(+) = 2/9 with the initial record; P(+) = 0 without it.
A different initial instrument with the SAME POVM effects gives P(+) = 4/9.
The initial effects alone cannot determine a sequential experiment.

## What is stored

Initial memory labels include the fixed stamp (epoch0,k1); terminal
labels include (epoch1,k0). These are predeclared scheduled-event stamps,
not an autonomous recorder of an arbitrary unknown clock reading.
Both outcome records are preserved after q=3 while their pointers are isolated.
Clock labels are 0,1,2,0,1,2; stored-event counts are 0,1,1,2,2,2.
Count refers only to the two declared events, not all phase ticks.

Conditioning on phase k=0 alone mixes the blank and completed histories:
P(terminal valid | k=0) = 1/2. Including epoch1 makes it 1.
A completely periodic blank-to-filled work state would require
W_cycle chi0 = chi0, but their overlap is zero and closure residual is sqrt(2).
This is not a no-go theorem for cyclic clocks with resetting or longer memory.

## Constraint and scope

The 162-dimensional positive open-history constraint has a unique zero
mode. Residual: 1.471e-16.
Its first nonzero eigenvalue is 0.0340742
in arbitrary constraint units, NOT a physical clock energy gap.
K_hist includes ordered recording gates and an initial boundary;
it replaces the noninteracting G-R1 constraint, rather than proving that
the old Hamiltonian automatically records events.

The final gate does not commute with the bare system Hamiltonian
(commutator norm 4.89898). No energy-conserving apparatus,
microscopic gate duration, thermodynamic irreversibility or physical
memory lifetime is certified. Reversing the protocol erases the records.
The pointer probabilities are classical-readable; complete decoherence
of all pointer density-matrix coherences is not assumed or proved.

## Next, not executed

G-R3: explicitly choose an Ei/Ed object and measure, then one interaction
and a two-clock observable. No density law, relative clock rate, metric,
mass formula, detector scan, or Route-F promotion is included here.

Run: python3 route_g/code/verify_gr2_two_event_records.py.
The JSON stores every check and hashes of unchanged prior code/results.
The TeX/PDF give the complete finite construction and its assumptions.
