---
title: Profile matching of powder patterns (Le Bail fit, fixed intensities)
keywords: profile matching, matching mode, le bail, pawley, integrated intensities, free intensities, fixed intensities, pd_peak, intensities of reflections, cell without structure, peak shape, jbt 2
categories: pd_peak
---
Profile matching is a procedure, like the Rietveld refinement: nothing
is written to the experiment except the results. It works for `pd`
experiments (1D, constant wavelength, polarized or unpolarized neutrons,
X-rays):

```
cryspy.rhochi_profile_matching(rcif)                          # free intensities
cryspy.rhochi_profile_matching(rcif, fixed_intensities=True)  # intensities of pd_peak
```

The pattern is

    y(2θ) = bkg(2θ) + Σ_hkl I_hkl P(2θ − 2θ_hkl) L(2θ) / L(2θ_hkl)

(P the pseudo-Voigt profile of the phase, L the Lorentz factor), where
I_hkl are the values `intensity_plus + intensity_minus` of the loop
`pd_peak` of the phase, the same as written by a Rietveld calculation.
The reflections come from the cell and the space group of the phase.

- Free intensities (default): at every step the intensities of all
  reflections of all phases are fitted to the pattern by linear least
  squares with I_hkl ≥ 0 (Pawley fit with the positivity constraint), for
  the current peak positions and shapes. The parameters marked for
  refinement (zero shift, cell, `U, V, W, X, Y`, `igsize`, asymmetry,
  background) are refined with the best intensities; with nothing marked
  only the intensities are fitted. They are written to `pd_peak`
  (`intensity_plus`, `intensity_minus`). Polarized neutrons (`chi2_sum`
  and `chi2_diff` True): up and down intensities are fitted separately
  (both ≥ 0); with the sum only they are equal halves. Reflections at the
  same 2θ in one phase share the intensity in proportion to their
  multiplicity; reflections without measured points get 0.
- `fixed_intensities=True`: the intensities of `pd_peak` (of every phase)
  are kept fixed, the marked parameters of the peak positions and shapes
  are refined. Intensities of a Rietveld calculation give the same pattern.

The atoms, the scale factor and the texture of the phases have no effect:
don't refine them. The magnetic satellites of k ≠ 0 structures,
and `magnetic_broadening` are not used. Other
experiments of the object (single crystal, `pd2d`, `tof`) are refined as
in the Rietveld refinement.

AI assistants: the tool `profile_matching`, the recipes
`profile_matching` (free intensities) and `refine_with_fixed_intensities`.
