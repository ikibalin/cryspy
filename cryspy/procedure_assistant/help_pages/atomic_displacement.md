---
title: Atomic displacement parameters (thermal vibrations)
keywords: adp, thermal, vibration, vibrations, displacement, debye-waller, temperature factor, b_iso, u_iso, biso, uiso, bani, uani, anisotropic, isotropic, тепловые колебания, тепловые параметры, температурный фактор, дебай-валлер, анизотропные, изотропные
categories: atom_site, atom_site_aniso
---
Thermal vibrations of atoms are described in the crystal data block by the
`atom_site` loop and, for anisotropic vibrations, the `atom_site_aniso` loop.

Isotropic vibrations: set `_atom_site_adp_type` to `Biso` or `Uiso` and give
`_atom_site_B_iso_or_equiv` (Å²) or `_atom_site_U_iso_or_equiv` (Å²),
B = 8π²U.

```
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
  Fe  Fe3+  0.0  0.0  0.0  Biso  0.35()
```

Anisotropic vibrations: set `_atom_site_adp_type` to `Uani` or `Bani` and add
the `atom_site_aniso` loop with the same labels. Give `U_11 ... U_23` for
`Uani` or `B_11 ... B_23` for `Bani`. Keep the isotropic value in `atom_site`
equal to zero, it is not used in this case.

```
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
  Fe  0.004()  0.004  0.006()  0.0  0.0  0.0
```

Symmetry constraints of the site are applied automatically: dependent
components follow the free ones and cannot be refined separately.

The types `Uovl`, `Umpe`, `Bovl` are accepted but treated as isotropic for now.
