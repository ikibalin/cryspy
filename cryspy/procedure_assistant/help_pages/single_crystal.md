---
title: Single crystal experiment (diffrn)
keywords: single crystal, diffrn, flip ratio, polarized, unpolarized, intensity, ub matrix, orientation matrix, reflections, монокристалл, флиппинг отношение, поляризованные нейтроны, ub-матрица
categories: diffrn_refln, diffrn_orient_matrix, diffrn_radiation, setup, phase, extinction, jana_twin
---
A single crystal experiment is a `diffrn` data block:

- `setup`: `wavelength` (Å), `field` (T), `radiation` (`neutrons` or `X-rays`);
- `diffrn_orient_matrix`: UB matrix `ub_11 ... ub_33`, `type` gives the
  convention (`CCSL` by default, `BusingLevy`, ...);
- `diffrn_radiation`: beam `polarization` and flipper `efficiency`
  (polarized neutrons);
- `diffrn_refln` loop: `index_h`, `index_k`, `index_l` and the measured data,
  either flip ratios `fr`, `fr_sigma` or intensities `intensity`,
  `intensity_sigma`; `excluded True` removes a reflection from chi-square;
- `phase`: `label` of the crystal and `scale`;
- optional `extinction`, `jana_twin` (twins), `chi2`.

Calculated values are written back as `fr_calc` or `intensity_calc`.
