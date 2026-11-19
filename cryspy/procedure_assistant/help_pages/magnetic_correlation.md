---
title: Magnetic correlation length (magnetic domain size) in powder
keywords: magnetic correlation length, magnetic domain size, magnetic size, broad magnetic peaks, short-range magnetic order, lorentzian, scherrer, anisotropic size, ellipsoid, interference, magnetic phase size fullprof, корреляционная длина, размер магнитного домена, уширение магнитных пиков, ближний магнитный порядок
categories: magnetic_correlation
---
The magnetic order of a powder grain is often correlated over a shorter
distance than its crystal structure (magnetic domains, short-range order):
the magnetic reflections are broader than the nuclear ones. In a `pd`
experiment (1D, constant wavelength) a phase can get its own magnetic width
with the loop `magnetic_correlation`:

```
loop_
_magnetic_correlation_phase_label
_magnetic_correlation_length
  Ho2BaNiO5  45.0(20)
```

The nuclear intensity |F_N|² keeps the profile of the phase (instrument,
`igsize`, the size and strain terms `X`, `Y` of the phase). The magnetic
intensity (|F_M⊥|², the paramagnetic terms) and the nuclear-magnetic
interference (polarized neutrons) are convolved in addition with a
Lorentzian, whose half width in Q = 4π sinθ/λ is 1/ξ, ξ being the
`length` in Å (spin correlations decaying as exp(-r/ξ)). In 2θ this is the
Lorentzian term `Y_mag / cos θ`, `Y_mag = 180 λ / (π² ξ)` degrees, the same
as the size term `Y` of FullProf for the apparent (Scherrer) size `2 ξ`.
So it is the FullProf trick of a separate "magnetic phase" with its own
size, but nuclear and magnetic scattering stay in one phase (one set of
moments, the interference term is kept). A length 0 means infinite.
The magnetic satellites of a structure with k ≠ 0 are broadened in the
same way.

When the width depends on the direction of hkl, the correlation lengths
form an ellipsoid: the components `aniso_11 ... aniso_23` (Å², the
convention of the `U_ij` of atoms) give the length along the scattering
vector `xi(h)² = Σ h_i h_j a*_i a*_j Xi_ij / |H|²`. `aniso_11`, `aniso_22`,
`aniso_33` are the squared lengths along a*, b*, c*. The tensor follows the
symmetry of the crystal class of the phase (of the magnetic structure for
k ≠ 0): the dependent components are calculated from the independent ones.
When the `aniso` components are given, `length` is not used.

```
loop_
_magnetic_correlation_phase_label
_magnetic_correlation_aniso_11
_magnetic_correlation_aniso_22
_magnetic_correlation_aniso_33
_magnetic_correlation_aniso_12
_magnetic_correlation_aniso_13
_magnetic_correlation_aniso_23
  TbTe3  3600(150) 3600(150) 400(30) 0. 0. 0.
```

Procedures: `add_magnetic_correlation_length` (isotropic or ellipsoid, sets
the refinement flags of the independent parameters) and
`report_magnetic_correlation` (ξ with errors, the Scherrer sizes, Y_mag,
the lengths along a*, b*, c* and the principal axes of the ellipsoid; for
example the correlations shorter along c). Recipe: `magnetic_correlation`.
Not implemented for `pd2d` and `tof`.
