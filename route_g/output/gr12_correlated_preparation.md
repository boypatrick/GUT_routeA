# G-R12 local correlated hopping and preparation energy ledger

189/189 checks pass. No fitted pulse or physical units.

The current C=i[K,D]/sqrt(2) is local on configuration-graph edges,
joint-translation and swap invariant, but not a density-only interaction.
It supplies the exact desired two-state preparation block without inserting
a Fourier-mode projector. Its coefficient/phase and control remain new inputs.

All fixtures use inherited g=0.1, alpha=0.25, and the same n=25 updates
(duration 6.25; ideal continuous duration 6.15479709). No optimization.

## Original contact retained throughout the direct pulse

| gamma | unconditional Dplus fidelity | outside triple | net controller work |
|---|---:|---:|---:|
| 0.0 | 0.999909367 | 8.26931e-31 | 1.2061134e-16 |
| 0.09 | 0.997140947 | 0.000178685 | -0.0099260795 |
| 0.03 | 0.999600916 | 1.97636e-05 | -0.0033309796 |
| 0.01 | 0.999875084 | 2.19518e-06 | -0.0011108886 |

The pulse is not energy invariant under the old full H. Controller work
equals the full system energy change. Fidelity is preserved after switch-off,
but continued control removes the prepared state; autonomous isolation is absent.
The unchanged G-R11 reader succeeds with this fidelity, and its unread-success
57/82 versus 25/82 backreaction remains. No failed branch is renormalized away.

## Explicit 21-state relay completion

Twelve off-resonant auxiliaries give -B†B/Delta = g*C - sqrt(2)*g*(I+D).
The density-dependent Stark term is leading order, not a small correction.
Compensation +sqrt(2)*g*D is separately declared; it is not silently dropped.

| compensated | Delta | effective F | full unconditional F | auxiliary probability | auxiliary energy | net controller work |
|---|---:|---:|---:|---:|---:|---:|
| False | 20 | 0.95183145 | 0.93561813 | 0.0249881 | 0.599715 | 0.496818 |
| False | 40 | 0.95183145 | 0.94438231 | 0.0125617 | 0.552714 | 0.499559 |
| False | 80 | 0.95183145 | 0.94935599 | 0.0042897 | 0.360335 | 0.340302 |
| True | 20 | 0.99960092 | 0.98003504 | 0.0188301 | 0.451921 | 0.373277 |
| True | 40 | 0.99960092 | 0.99051139 | 0.00886303 | 0.389973 | 0.351192 |
| True | 80 | 0.99960092 | 0.99631112 | 0.00333585 | 0.280212 | 0.263537 |

Small auxiliary population does not make the energy ledger negligible:
its energy is (4+Delta)*P_aux. The displayed positive full-model switching
work cannot be replaced by the direct effective pulse's negative work.
The sudden-switch protocol counts all high-energy exits, not just a normalized
low-energy conditional state. No measured hardware or autonomous drive is claimed.

## Bounded decision

The local-operator existence and declared relay audit are complete. The old
2/3 obstruction is specific to a contact family, not joint translation alone.
A physical controller/mediator choice is still required for an implementation.
If pursuing it, test one energy-accounted smooth switching protocol or a passing
trigger; otherwise close this as a controlled preparation model. Do not add
another parameter fit, claim a time arrow, or promote any Route-F gate.
