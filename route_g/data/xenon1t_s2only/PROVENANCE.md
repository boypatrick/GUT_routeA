# XENON1T S2-only NR response: pinned primary-source subset

Retrieved 2026-09-23. This is the **2019 analysis / 2020 public release**, not the latest xenon result. It was chosen because the collaboration provides a small, machine-readable nuclear-recoil response matrix with explicit normalization and usage instructions. No plot digitization or third-party reconstructed efficiency is used.

## Source and unchanged files

- Collaboration paper: E. Aprile et al., *Light Dark Matter Search with Ionization Signals in XENON1T*, [Phys. Rev. Lett. 123, 251801 (2019)](https://doi.org/10.1103/PhysRevLett.123.251801), [arXiv:1907.11485v2](https://arxiv.org/abs/1907.11485v2).
- Official public-data index: [XENON collaboration](https://xenonexperiment.org/public-data/).
- Official repository: [XENON1T/s2only_data_release](https://github.com/XENON1T/s2only_data_release).
- Pinned commit: [`5a364bc8709f2561e5a013ddea6993a5a7c8e313`](https://github.com/XENON1T/s2only_data_release/commit/5a364bc8709f2561e5a013ddea6993a5a7c8e313), dated 2020-10-09. That commit fixes ER bin edges; the NR response is preserved unchanged.

Raw download prefix:

```text
https://raw.githubusercontent.com/XENON1T/s2only_data_release/5a364bc8709f2561e5a013ddea6993a5a7c8e313/
```

| Local file | Upstream relative path | Bytes | SHA-256 |
|---|---|---:|---|
| s2_response_nr.csv | s2_response_nr.csv | 170120 | 76098a792eaf576ac835fcab7c412b2cc7763f2808c93db93998e256684dc202 |
| s2_binning_info.csv | s2_binning_info.csv | 8032 | 6ee893b8dff0a92848369026ec34ddfbeeac583ed6b9cdb5a9b2e1331cdbd341 |
| UPSTREAM_README.md | README.md | 7750 | e2da8940580aaff299ba89cfeb9e733ec50baae02c7f7eb69d39123d013f3797 |
| UPSTREAM_ZENODO.json | .zenodo.json | 831 | d713acec87130aef77719aa547fd07579dc97257eaaa20b0d648ea7ce2b98f24 |
| UPSTREAM_example_analysis.ipynb | example_analysis.ipynb | 83502 | 1cec3d6f4a6251e7233dddb91b680e69f46522cec08eb204cc42bab0780a1610 |

The five files are byte-identical downloads; only their local `UPSTREAM_` names differ where indicated. No events, background templates, limit curves or fitted model parameters are copied into this subset. The notebook is retained as usage provenance, not executed as a new analysis.

## License and attribution

The preserved upstream `.zenodo.json` identifies the creator as **XENON Collaboration**, release version 1.0 dated 2020-08-13, access `open`, license `other-open`. The upstream description expressly offers the release for researchers' own model calculations and requests citation of the paper. Preserve that metadata and attribution. **Do not relabel this release as CC-BY-4.0 or MIT**: those labels describe other XENON repositories, not this preserved license declaration. This local note does not expand or replace the upstream terms.

## Matrix schema and normalization

The NR CSV contains 41 monoenergetic response rows, from 0.4 to 50.0 keV, and 199 S2-bin columns. The first columns are `energy_kev`, `energy_bin_start_kev`, `energy_bin_end_kev`; the remaining columns are `s2_bin_000` through `s2_bin_198`. Entries are fractions of full-target recoils entering each observed S2 bin after selections; row sums are not one. The binning CSV supplies `s2_bin_number`, `start_pe`, `end_pe`, `linear_center_pe`, `log_center_pe`.

The official raw search exposure is **356770 kg day = 0.97678 tonne year**. The response already includes position/fiducial and event-selection losses: do not multiply a second fiducial fraction or replace it with an effective exposure. Follow [the preserved upstream README](UPSTREAM_README.md) and notebook when folding a rate spectrum.

## Declared bounded usage for Route G

The paper uses a 0.7 keV true-NR cutoff; its default charge-yield extrapolation below this was not promoted to measured response. The distributed matrix deliberately omits that cutoff. Route G must explicitly restrict true recoil integration to **0.7–50 keV** and must not extend its response outside the supplied range. Clip outer energy-bin edges accordingly. This conservative window excludes both the hypothetical low-energy extension and unavailable high-energy response.

The paper excludes S2 below 150 PE. The supplied bin 28 spans [147.4066,150.027] PE, so an exact-bin conservative choice is **bins 29–198**, namely **[150.027,3000] PE**. This drops the small straddling-bin fragment without assuming its sub-bin shape. Sum the selected response columns to obtain an accepted-fraction curve. This is a declared acceptance window, **not** one of the paper's model-specific optimized statistical ROIs. Any alternative interpolation or boundary-bin treatment must be recorded explicitly.

The paper describes approximately 2 tonnes of active xenon, a 47.9 cm TPC radius and approximately 97 cm drift height; its search sample has 180.7 live days. These rounded geometry/time values are useful context, not replacements for the release's exposure normalization. Strong depth and radius selections mean there is no single extra fiducial mass factor to apply to this response. [Paper, data selection and detector response](https://arxiv.org/html/1907.11485v2).

## What a Route-G fold may and may not claim

For a separately derived nuclear differential rate, a valid use is

\[
N_b=\mathcal E\int_{0.7\,\mathrm{keV}}^{50\,\mathrm{keV}}
dE_R\,\frac{dR}{dE_R}\,A_b(E_R),
\qquad \mathcal E=356770\ \mathrm{kg\,day}.
\]

The response is averaged over the original target-volume distribution. Reusing it for boosted particles requires a thin, broadly/uniformly illuminated target and ordinary elastic nuclear recoils with the same charge-yield response. A narrow localized beam, attenuation, multiple scattering, inelastic nuclear excitation, or accompanying electronic activity requires a new spatial/event response. Xenon isotope weighting and nuclear form factors belong in the incident scattering model, not in this efficiency matrix.

Route G has not supplied an incoming physical flux. A response-weighted cross section or required **at-target** flux/fluence is therefore meaningful; an actual predicted event yield, exclusion or discovery sensitivity is not. The original halo-WIMP cross-section curves cannot be applied to Route-G boosted particles merely because both use scalar scattering. Their velocity and abundance assumptions differ. Do not use a zero-background 2.3-event limit for this sample, whose backgrounds are nonzero and incomplete.

For boosted kinematics, record the recoil endpoints and the in-window cross section, and explicitly label below-0.7-keV and above-50-keV contributions as outside this response fold. Do not assign physical fractions relative to a total cross section by extending the low-transfer coherent model to unvalidated high recoil energies. The in-window accepted result is a **coverage-limited contribution**, not an assertion that the detector is physically blind outside that window. No extrapolation of acceptance to keV–MeV recoils is authorized by these files. Detector overburden and source-to-target transport are also outside an at-target fold.

Finally, a xenon NR response does not by itself calibrate the former Gaussian estimator of the incident X momentum or measure its compact recoil sign. Those inverse-observation requirements remain separate from observing a target recoil.
