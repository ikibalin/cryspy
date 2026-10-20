---
title: Paramagnetic susceptibility (local susceptibility approach)
keywords: susceptibility, paramagnetic, chi, chi_11, local susceptibility, field, flip ratio, magnetization ellipsoid, восприимчивость, парамагнетик, локальная восприимчивость, эллипсоид намагниченности
categories: atom_site_susceptibility, atom_site_exchange, setup
---
In the paramagnetic state the induced moment of an ion is M = χ·H with a
local susceptibility tensor χ per site. It is set by the
`atom_site_susceptibility` loop of the crystal:

- `label` as in `atom_site`;
- `chi_type` is `Ciso` (only `chi_11` is used) or `Cani` (full tensor,
  symmetry constraints of the site are applied);
- `chi_11 ... chi_23` in μB/T, in the Cartesian system x ∥ a*, z ∥ c.

The magnetic form factor is taken from `atom_site_scat`. The field in tesla is
`_setup_field` of the experiment; for a single crystal the field direction is
along the vertical axis given by the UB matrix (`diffrn_orient_matrix`).

Flip ratios are measured in `diffrn_refln` (`fr`, `fr_sigma`); polarization and
flipping efficiency are in `diffrn_radiation`.

After refinement `report_main_axes_of_magnetization_ellipsoids` gives the main
axes of the magnetization ellipsoids.
