---
title: Matching mode of powder patterns (profile matching, Le Bail fit, fixed intensities)
keywords: matching mode, profile matching, le bail, pawley, integrated intensities, free intensities, fixed intensities, pd_peak, intensities of reflections, cell without structure, peak shape, jbt 2
categories: pd_matching, pd_peak
---
In the matching mode of a `pd` experiment (1D, constant wavelength,
polarized or unpolarized neutrons, X-rays) the integrated intensities of
the reflections are not calculated from the crystal structure. The item
`pd_matching` switches it on:

```
_pd_matching_mode refine_intensities
```

or `experiment.switch_matching_mode("refine_intensities")`. The pattern is

    y(2θ) = bkg(2θ) + Σ_hkl I_hkl P(2θ − 2θ_hkl) L(2θ) / L(2θ_hkl)

(P the pseudo-Voigt profile of the phase, L the Lorentz factor), where
I_hkl are the values `intensity_plus + intensity_minus` of the loop
`pd_peak` of the phase, the same as written by a Rietveld calculation.
The reflections come from the cell and the space group of the phase.

- `refine_intensities`: at every calculation the intensities of all
  reflections of all phases are fitted to the pattern by linear least
  squares with I_hkl ≥ 0 (Pawley fit with the positivity constraint), for
  the current peak positions and shapes. Refine the zero shift, the cell,
  the resolution `U, V, W, X, Y`, `igsize`, the asymmetry and the
  background as usual: they are refined with the best intensities at
  every step. The fitted intensities are written to `pd_peak`
  (`intensity_plus`, `intensity_minus`). Polarized neutrons (`chi2_sum`
  and `chi2_diff` True): up and down intensities are fitted separately
  (both ≥ 0); with the sum only they are equal halves. Reflections at the
  same 2θ in one phase share the intensity in proportion to their
  multiplicity; reflections without measured points get 0.
- `fixed_intensities`: the intensities of `pd_peak` are kept fixed, the
  peak positions and shapes are refined. The loop must be given for every
  phase (after a calculation, a profile matching or from the user).
  Intensities of a Rietveld calculation give the same pattern.
- `off`: the Rietveld calculation again.

The atoms, the scale factor and the texture of the phases have no effect
in the matching mode: don't refine them. The magnetic satellites of k ≠ 0
structures (`pd_peak_magnetic_k`), `magnetic_broadening` and
`layer_correlation` are not used. Only the `pd` experiments (not `pd2d`,
`tof`) have the matching mode.

Recipes: `profile_matching` (free intensities), `refine_with_fixed_intensities`.
