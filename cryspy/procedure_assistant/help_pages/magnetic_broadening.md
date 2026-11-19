---
title: Extra broadening of magnetic reflections (magnetic domain size) in powder
keywords: magnetic broadening, magnetic domain size, magnetic correlation length, magnetic size, broad magnetic peaks, narrow magnetic peaks, short-range magnetic order, lorentzian, scherrer, anisotropic size, interference, magnetic phase size fullprof
categories: magnetic_broadening
---
The magnetic order of a powder grain is often correlated over a shorter
distance than its crystal structure (magnetic domains, short-range order):
the magnetic reflections are broader than the nuclear ones. In a `pd`
experiment (1D, constant wavelength) a phase can get an extra magnetic
width with the loop `magnetic_broadening`:

```
loop_
_magnetic_broadening_phase_label
_magnetic_broadening_width
  Ho2BaNiO5  0.0125(4)
```

The nuclear intensity |F_N|² keeps the profile of the phase (instrument,
`igsize`, the size and strain terms `X`, `Y` of the phase). The magnetic
intensity (|F_M⊥|², the paramagnetic terms) and the nuclear-magnetic
interference (polarized neutrons) are convolved in addition with a
Lorentzian, whose half width in Q = 4π sinθ/λ is `width` (1/Å):

    width = 1/ξ_mag − 1/ξ_nucl

(ξ are the magnetic and nuclear correlation lengths, correlations decaying
as exp(-r/ξ)).

- `width = 0` (default): the magnetic reflections have the nuclear profile,
  the same result as without the loop;
- `width > 0`: the magnetic correlations are shorter (smaller magnetic
  domains); the extra Scherrer size is 2/width;
- `width < 0`: they are longer; the magnetic peaks get narrower, at most by
  the Lorentzian part of the nuclear width.

Effect on the diffraction signal: the integrated intensity of every
reflection stays the same (the Lorentzian is normalized), only its shape
changes. With `width > 0` the magnetic peaks (and the magnetic part of
the peaks where nuclear and magnetic reflections coincide) get lower and
broader with long Lorentzian tails, so a broad magnetic peak sits on a
sharp nuclear one; the nuclear-only peaks do not change. The extra width
in 2θ grows as 1/cos θ, so high-angle magnetic peaks broaden more in 2θ,
while in Q the extra width is the same for all of them (or depends only
on the direction of hkl for `width_ij`). With polarized neutrons the
flipping difference (the interference term) is broadened as well. The
satellites of a structure with k ≠ 0 are purely magnetic and are
broadened entirely.

In 2θ this is the Lorentzian term `Y_mag / cos θ`, `Y_mag = 180 λ width / π²`
degrees, the size term `Y` of FullProf. So it is the FullProf trick of a
separate "magnetic phase" with its own size, but nuclear and magnetic
scattering stay in one phase (one set of moments, the interference term is
kept). The magnetic satellites of a structure with k ≠ 0 are broadened in
the same way.

When the width depends on the direction of hkl, it is a quadratic form of
the direction: the components `width_11 ... width_23` (1/Å, the convention
of the `U_ij` of atoms) give `width(h) = Σ h_i h_j a*_i a*_j W_ij / |H|²`.
`width_11`, `width_22`, `width_33` are the widths along a*, b*, c*. The
tensor follows the symmetry of the crystal class of the phase (of the
magnetic structure for k ≠ 0): the dependent components are calculated
from the independent ones. When the `width_ij` are given, `width` is not used.

```
loop_
_magnetic_broadening_phase_label
_magnetic_broadening_width_11
_magnetic_broadening_width_22
_magnetic_broadening_width_33
_magnetic_broadening_width_12
_magnetic_broadening_width_13
_magnetic_broadening_width_23
  TbTe3  0.004(1) 0.004(1) 0.031(2) 0. 0. 0.
```

Here the magnetic correlations are shorter along c*.

Procedures: `add_magnetic_broadening` (isotropic or depending on hkl, sets
the refinement flags of the independent parameters, start value 0) and
`report_magnetic_broadening` (width with errors, Y_mag, the extra Scherrer
size, the limit of a negative width, the widths along a*, b*, c* and the
principal axes). Recipe: `magnetic_broadening`. Not implemented for `pd2d`
and `tof`.
