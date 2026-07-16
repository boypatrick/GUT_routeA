# AP-E6 same-soliton Yukawa--Callias audit

- Checks: **25/25**
- All audit checks pass: **True**
- Same physical CP1 moduli derived: **False**
- Actual Yukawa--Callias operator derived: **False**
- Determinant line O(2) derived: **False**
- Restricted Higgsed stack descent proven: **False**
- Restricted local differential-character ansatz typed: **True**
- Full same-soliton lane closed: **False**

## Canonical Workline-B provenance

- Canonical profile checksum (r,F,s): `81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59`.
- Workline-B script SHA-256: `393a1a1948b65bceb3c04167e452140a643a482a5035ed419d892fdd9cc1f10d`.
- Workline-B card SHA-256: `f5a353011e0877e3659e2efdfde8027884eb2d6487ca657458f16e4d1915f63f`.
- Canonical grid: box=16.0, points=4097, s(0)=0.917128964607058.
- The Wilson operator uses only the canonical F samples (linear off-grid interpolation). The canonical s field is not silently discarded: it is decoupled by the explicitly declared Yukawa term M U^{gamma5}, which contains no s coupling.

## Decisive results

- Relaxed hedgehog orbit: `SU(2)/Z2 = SO(3)` (dimension 3), not CP1.
- Uniform scalar-vacuum rotation norm scales as L^3.000000000000; it is not a localized normalisable mode.
- For the explicitly stated Callias domain and coercivity assumptions, the bounded Hermitian potential tends to beta M; its positive boundary eigenbundle is constant, so the ordinary S2-infinity Callias index is zero and cannot yield the requested rank-one O(2) family.
- The four-dimensional interpolation index/spectral flow can equal B, but it does not imply a static endpoint zero mode.
- Finite Wilson control gaps: [0.010896960821179729, 0.28367497127235386]
- Spectator Veronese control Chern number: 1.999999999999986.

## Gate boundary

The machine card deliberately passes its mathematical and numerical regressions while leaving every physical promotion gate false. A degree-one Route-E portal is not authorized.
