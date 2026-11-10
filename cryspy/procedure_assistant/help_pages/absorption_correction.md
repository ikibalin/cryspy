---
title: Absorption correction
keywords: absorption, absorption correction, mu, transmission, crystal size, ellipsoid, sphere, gaussian, T_min, T_max, абсорбция, поглощение, поправка на поглощение, размер кристалла
categories: exptl_absorpt, exptl_crystal
---
The absorption of neutrons in a single crystal is described by two items of
a single crystal experiment (`diffrn`), named as in the core CIF dictionary:

```
_exptl_absorpt_correction_type gaussian

_exptl_crystal_size_max 3.0
_exptl_crystal_size_mid 1.0
_exptl_crystal_size_min 0.5()
_exptl_crystal_ellipsoid_euler_alpha 0.0
_exptl_crystal_ellipsoid_euler_beta 40.0
_exptl_crystal_ellipsoid_euler_gamma 0.0
```

- `correction_type`: `none`, `sphere` (radius `_exptl_crystal_size_rad`) or
  `gaussian` (ellipsoid, numerical Gaussian integration over its volume).
- mu, the linear absorption coefficient (1/mm, absorption and scattering),
  is calculated from the measured crystal, the phase of the experiment
  (`_phase_label`): its `atom_site` (multiplicity x occupancy) and cell, the
  wavelength of `setup` (1/v law), and at lambda/2 for the lambda/2
  contamination. `calc_absorption_coefficient` reports it. Give
  `_exptl_absorpt_coefficient_mu` only to replace it: the 1/v law is wrong
  for Gd, Sm, Eu, Cd, Dy, or a measured mu.
- `size_max`, `size_mid`, `size_min`: full lengths of the axes of the
  ellipsoid in mm; the Euler angles (degrees, z-y-z) orient them in the
  Cartesian system x || a*, z || c. `ExptlCrystal.set_axes_by_directions`
  sets them from the directions of the longest and middle axes.
- `correction_T_min`, `correction_T_max`: calculated transmission range.

The transmission of each reflection multiplies the calculated intensities
in the refinement (RhoChi) and in MEM. The paths of the beams: normal beam
geometry (rotation about the vertical axis of U, as for polarized neutrons;
`_exptl_absorpt_diffracted_beam left|right`) or, when the reflections have
`_diffrn_refln_angle_omega`, `_diffrn_refln_angle_chi`, `_diffrn_refln_angle_phi`, a four-circle
measurement (rotation matrices of Busing & Levy, omega is the angle of the
sample, as the motors give it). Data of a four-circle diffractometer without
these columns are treated as normal beam data: give the angles for an
ellipsoid. At `ellipsoid_euler_beta` = 0 only alpha + gamma matters: do not
refine both. Refine the shape only when the test of the equivalents is
significant; mu and the sizes correlate (only mu x size matters).

Flip ratios do not change by absorption (it cancels in I+/I-) except for the
lambda/2 contamination. A sphere changes intensities with 2theta only and
correlates with B_iso; a flat or long crystal changes equivalent reflections
differently.

`estimate_absorption_correction` (assistant tool `estimate_absorption`) tells
whether a correction is needed; `apply_absorption_correction` switches it on;
the recipe `absorption_correction` does both and refines the shape.
