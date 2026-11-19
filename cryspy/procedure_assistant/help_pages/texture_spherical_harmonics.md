---
title: Texture by spherical harmonics
keywords: texture, preferred orientation, spherical harmonics, harmonic order, texture index, fibre texture, cylindrical symmetry, sample symmetry, pole figure, ODF, GSAS-II, goniometer
categories: texture_sh, texture_sh_coefficient, sample_goniometer
---
The texture of a phase in a powder experiment (`pd`, `pd2d`, `tof`) can be
described by the generalized spherical harmonics, as in GSAS-II (Bunge,
1982; Von Dreele, J. Appl. Cryst. 30 (1997) 517). The integrated intensity
of the reflection h is multiplied by

    A(h, y) = 1 + sum_l 4 pi/(2l+1) sum_m sum_n C_l^mn k_l^m(y) k_l^n(h),

l = 2, 4, ..., order_max. k_l^n(h) are the harmonics of the crystal
direction h symmetrized by the Laue class of the phase, k_l^m(y) the
harmonics of the direction y of the scattering vector in the sample,
symmetrized by the sample symmetry. A random powder has all C = 0 and A = 1.

```
loop_
_texture_sh_phase_label
_texture_sh_order_max
_texture_sh_sample_symmetry
  b2  6  cylindrical

loop_
_texture_sh_coefficient_phase_label
_texture_sh_coefficient_l
_texture_sh_coefficient_m
_texture_sh_coefficient_n
_texture_sh_coefficient_value
  b2  4  0  1  -2.42(5)
  b2  6  0  1  -7.2(1)
```

- `sample_symmetry`: symmetry of the sample, not of the crystal.
  `cylindrical` (fibre texture: wire, rod, powder pressed in a capillary or
  a flat holder, powder aligned by a field; only m = 0) is what one pattern
  can determine. `mmm` (rolled sheet), `2/m` (shear) and `-1` (none) need
  several orientations of the sample or a 2D detector.
- `order_max`: 4 to 8 for a weak texture, up to 12-16 for a strong one.
  The number of coefficients grows quickly; refine the low orders first.
- The allowed (l, m, n) follow from the Laue class of the phase, found from
  its symmetry operations (4-, 3- and 6-fold axes along c, the
  rhombohedral 3-fold along [111]). `add_texture_spherical_harmonics`
  writes them with the value 0 and the refinement flag.
- Texture index J = 1 + sum C^2/(2l+1): 1 for a random powder.
  `report_texture_spherical_harmonics` gives J and the range of the
  correction (inverse pole figure of the cylinder axis, multiples of the
  random distribution, MRD). A negative minimum is unphysical: lower
  `order_max`.
- The same coefficients mean the same texture in GSAS-II (crystal frame
  z || c*, x || a*; x || a for the hexagonal axes; monoclinic unique b:
  z || b*, x || c*). For monoclinic and triclinic cells with an oblique
  angle far from 90 degrees GSAS-II uses an approximate azimuth of the
  reflections, and its coefficients differ a little. For the Laue class
  -31m (P-31m, P-31c) GSAS-II takes the harmonics of -3m1, which do not
  have the symmetry -31m for n = 3, 9, ...; CrysPy uses the harmonics of
  -31m, so these coefficients are not the same.
- Do not use the March-Dollase `texture` and `texture_sh` for the same phase.
- The coefficients are written in every experiment of the phase. When
  `add_texture_spherical_harmonics` adds the texture of a phase to several
  experiments (banks of a time-of-flight instrument, sectors of a 2D
  detector, orientations of one sample), it ties their coefficients with
  equal constraint marks, `-2.4()b2_c4_0_1`, so that the experiments refine
  one texture as the phase texture of GSAS-II. Remove the marks for
  different samples.

Orientation of the sample (`sample_goniometer`, optional)
---------------------------------------------------------

```
_sample_goniometer_geometry debye_scherrer
_sample_goniometer_omega 0.0
_sample_goniometer_chi 90.0
_sample_goniometer_phi 0.0
_sample_goniometer_azimuth 0.0
```

Without this item the cylinder axis of the texture is vertical (vertical
capillary or wire, field of a 2D experiment): a capillary or a wire needs
nothing else. A flat plate in Bragg-Brentano (reflection) geometry measures
the scattering vectors normal to the plate, and the texture axis is the
normal: give it with

```
_sample_goniometer_geometry bragg_brentano
_sample_goniometer_chi 0.0
```

With the vertical default the fit of a flat plate is the same, but the
coefficients, J and the inverse pole figure belong to the vertical axis,
not to the normal of the plate, and the pole figure may look unphysical. The angles (degrees) are those of an Eulerian cradle in the
convention of Busing & Levy (Acta Cryst. 22 (1967) 457), the same as for
the four-circle single crystal experiments of CrysPy:

- laboratory system: y along the incident beam, z vertical (up), x
  horizontal, x = y x z;
- the sample frame (phi frame) is turned from the laboratory by
  R = Omega Chi Phi; omega and phi turn about z, chi about y (at omega = 0);
  a vector q of the laboratory is q_phi = R^T q in the sample frame;
  omega is the angle of the sample (as the motor gives it), not omega - theta;
- the axis of the cylindrical texture (psi = 0) is the x axis of the
  sample frame: with omega = chi = phi = 0 it is horizontal and normal to
  the beam, with chi = 90 it is vertical; the azimuth of the sample
  harmonics (m != 0) is counted from -z towards +y of the sample frame;
- `azimuth`: angle of the detector around the beam, from +x towards +z
  (0: scattering in the horizontal plane). For `pd2d` the angle phi of
  every point of the detector is added to it;
- `geometry`: `debye_scherrer` (transmission, default) or `bragg_brentano`
  (omega follows theta: omega + theta);
- `_texture_sh_sample_omega`, `_texture_sh_sample_chi`,
  `_texture_sh_sample_phi`: zeros of the
  orientation of the texture frame in the sample (refinable; GSAS-II
  "sample orientation zeros").

These are the sample parameters Omega, Chi, Phi, Azimuth of GSAS-II, which
uses the same convention: with the same angles both programs give the same
correction. The defaults differ: chi = 90 here, chi = 0 in GSAS-II (set
`_sample_goniometer_chi 0.0` to take the coefficients of a GSAS-II project
with its default angles).
