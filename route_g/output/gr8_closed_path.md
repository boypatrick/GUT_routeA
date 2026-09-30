# G-R8: specified closed guided protocol

80/80 checks pass. Synthetic control card, not observations.

## One complete sequence and phase convention

Prepare the clock at the common height; split two orthogonal guide modes;
5 ms smooth outward motion, 10 ms hold, 5 ms return; apply the inverse
mode coupler and energy-resolved readout at the common height.
The symmetric 1 mm guide separation encloses A=d*(hold+ramp)=1.5e-5 m s.
Both ramps contribute. The guide is defined with a global-z linear
support/acceleration potential and zero branch scalar offsets.
Its ground-state phase difference vanishes by the complete action,
not by deleting the ground rest energy. Total port phases are
Phi_g=alpha_f-alpha_s and Phi_e=Phi_g-delta_i in this model.

| Clock | Total delta (rad) | Ideal 1-V | Ideal negativity |
|---|---:|---:|---:|
| E2 | 7.0788937e-06 | 6.264e-12 | 1.77e-06 |
| E3 | 6.6033949e-06 | 5.451e-12 | 1.651e-06 |

## Closure is checked, not assumed

The compensated harmonic guide closes ground wave packets exactly.
The same guide with the excited inertial mass is not exactly closed
to all orders. Integer oscillator periods cancel its first-order
residual; finite-mass residuals and Gaussian overlaps are in the JSON.
An independent ODE at an explicitly artificial mass difference checks
the convolution and phase boundary term. This mass-only diagnostic
does not establish all-orders general relativity or real trap accuracy.

## Pure/dephased controls and joint readout

After calibrated local phase rotations, the pure state is
exp(i*delta*Zp*Zc/4)|++>. Its energy-dephased control has the same
path visibility but is separable. The witness
W=(I-Xp*Xc+Zp*Yc+Yp*Zc)/4 is a partial transpose of a positive
projector and is nonnegative on every separable state.
For positive small delta its ideal pure value is -sin(delta/2)/2;
the dephased value is (1-sin(delta/2))/4. XX, ZY and YZ use
incompatible local axes but commute globally. Add ZZ as an actually
noncommuting joint control: [ZY,ZZ]=2i IX. Analysis pulses followed
by local population detection implement the specified observables.

## Physical priority, not another derivation gate

The full guide potential, not just the trajectories or forces, fixes
the common phase. Replacing z by z-q in its linear term leaves forces
unchanged but changes the ground fringe. The declared support field
does not cancel the clock-energy gravitational coupling.
Prioritize calibrated state-resolved phase. Ideal witness shot counts
are already enormous at this geometry, before finite lifetime, controller
records, pulse errors, radiation noise and count-rate limitations.
Do not infer operational feasibility or observed entanglement from the
small mass-closure residual. E2's finite lifetime still matters.
Next select/calibrate one actual guide/coupler or reject this protocol
on a measured noise budget. No new time law or Route-F promotion.
