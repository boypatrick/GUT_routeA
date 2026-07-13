# Knu target map

No web lookup was used.

The input is the combined e-pi projection profile from the previous replay.

## Preferred conservative channel targets

| channel | class | current margin | future margin | future amp suppression | future S_T max |
|---|---:|---:|---:|---:|---:|
| `LLLL_downdownup_Knu` | `Knu` | 1.339617e+00 | 3.215080e-01 | 1.763616e+00 | 4.252626e-06 |
| `LLLL_upupdown_K0mu` | `K0mu` | 2.075664e+01 | 2.075664e+00 | 1.000000e+00 | 1.080537e-05 |
| `LLLL_upupdown_Knu` | `Knu` | 1.113033e+00 | 2.671278e-01 | 1.934819e+00 | 3.876331e-06 |
| `LLLL_uude_epi` | `e_pi` | 4.208568e+00 | 1.010056e+00 | 1.000000e+00 | 7.537617e-06 |
| `RRRR_uude_epi` | `e_pi` | 1.249765e+02 | 2.999437e+01 | 1.000000e+00 | 4.107534e-05 |
| `RRRR_uusd_anycharged` | `Knu` | 1.228658e+01 | 2.948778e+00 | 1.000000e+00 | 1.287900e-05 |

## Verdict

Current bottleneck: `LLLL_upupdown_Knu` with margin 1.113033e+00.
Future 1e35 bottleneck: `LLLL_upupdown_Knu` with margin 2.671278e-01.
Future-safe S_T target: 3.876331e-06; equivalent amplitude suppression: 1.934819e+00.

After the e-pi projection sensitivity replay, the preferred conservative bottleneck is LLLL_upupdown_Knu.  The current 2.4e34 yr stress passes with margin 1.113, but a uniform 1e35 yr stress requires either S_T <= 3.876e-6 or a common Knu amplitude suppression of about 1.935.
