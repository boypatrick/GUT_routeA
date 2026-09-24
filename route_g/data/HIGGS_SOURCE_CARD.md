# G2-H external source card

Retrieved/checked 2026-09-23. A source normalization card, not measured dark flux or a new calibration dataset.

| Quantity | Value | Status and primary basis |
|---|---:|---|
| Collision energy | 13 TeV | Actual Run-2 pp collisions |
| Integrated luminosity at one interaction point | 139 fb^-1 | Published ATLAS dataset |
| Higgs cross section | 55.6 ± 2.5 pb | SM prediction quoted in the production paper; not a model-independent measured dark-source rate |
| Invisible branching allowance | 0.107 | Conditional observed ATLAS 95% CL limit, including earlier-energy combination |
| Xenon active mass | about 2 tonnes | Rounded XENON1T active target |
| Target radius / height | 47.9 / 97 cm | Published cylindrical TPC dimensions |
| Generous stress inputs | 60 pb, fN=.326, column2e24 nuclei/cm² | Declared envelope, not combined measurement/confidence interval |

Primary sources:

- ATLAS, [arXiv:2207.08615v2](https://arxiv.org/html/2207.08615v2), JHEP05(2023)028, sections1,2,5. The luminosity and SM normalization supply N_h=7.7284 million. Visible-channel cross-section extraction assumes SM branching fractions; do not independently reuse that extraction after allowing invisible decays. The SM-production premise is instead stated explicitly. Its Higgs mass convention differs slightly from the retained G2-X mh; no precision interpolation smaller than the quoted theory uncertainty is claimed.
- ATLAS, [arXiv:2301.10731v2](https://arxiv.org/abs/2301.10731v2), Phys.Lett.B842(2023)137963. The invisible allowance requires SM production and escaping invisible particles. It limits a hypothetical source; it does not observe X.
- XENON, [arXiv:1907.11485v2](https://arxiv.org/html/1907.11485v2), introduction and event selections. Target is at Gran Sasso; no LHC-to-XENON installation or transport is claimed. The maximum cylinder chord is only a generous single-pass column bound.
- Actual NR/S2 files remain unchanged in [the pinned response subset](xenon1t_s2only/PROVENANCE.md). Their exposure is not multiplied into an already time-integrated collider particle budget. The full-volume averaged response is not claimed valid for arbitrary narrow-beam positions.

The chosen process is pp→h+anything, h→X_j anti-X_j in the drive-off vacuum. It does not reproduce the old selected Phi-X collision source or construct its GeV-scale pump. Sum every open signed complex mode with the common portal coupling. The source and detector couplings must not be varied independently.

No realistic Higgs boost/solid-angle/transport distribution is fabricated. A general Lorentz-boosted source-to-fluence functional is given in the TeX. The principal verdict uses a stronger, distribution-independent envelope: every produced particle gets the longest xenon chord, and every retained-window elastic recoil gets efficiency one. Additional geometric and response losses can only reduce this single-pass, nonmultiplying signal. Electronic/inelastic and outside-window contributions remain uncalculated.

The continuous count ceiling is not an achievable point or a best-fit radius. Factors in an upper envelope need not attain their individual suprema together. Tests very near the closing threshold are labeled regression-only, not candidate benchmarks.
