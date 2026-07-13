# Full Knu width calibration

No web lookup was used.

The audited formula is

```text
Gamma = K_dyn * K_had * (S_T |A_dress|)^2
```

| label | pair | phase | amplitude | future margin | K_dyn [yr^-1] | S_T max at 1e35 |
|---|---:|---:|---:|---:|---:|---:|
| `soft_spectrum_width` | `23` |  | 2.361106e-11 | 3.100993e-01 | 5.377683e-01 | 4.176492e-06 |
| `projected_worst_phase_width` | `13` | 4.188790e+00 | 2.543937e-11 | 2.671278e-01 | 5.377683e-01 | 3.876331e-06 |

## Cross-checks

- K_dyn max relative spread: 0.000000e+00
- projected/soft future-margin ratio: 8.614268e-01
- soft/projected amplitude ratio squared: 8.614268e-01
- ratio mismatch: 1.110223e-16

## Verdict

The preferred soft-spectrum bottleneck and the final projected Knu target row are described by the same calibrated width constant to numerical precision.  The remaining gap is therefore not a hidden normalization mismatch; it is a physical need for an additional Knu amplitude suppression or a stronger triplet filter in the same width formula.

Final future margin is 2.671278e-01; additional amplitude suppression needed is 1.934819e+00.
