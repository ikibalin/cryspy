---
title: Magnetic form factor in the dipole approximation
keywords: form factor, magnetic form factor, lande, kappa, j0, j2, dipole approximation, atom_site_scat, атомный форм-фактор, магнитный форм-фактор, фактор ланде
categories: atom_site_scat, atom_type_scat
---
The spherical magnetic form factor is f = <j0> + (2/g − 1) <j2>, where g is
the Landé factor. It is defined by the `atom_site_scat` loop of the crystal:

```
loop_
_atom_site_scat_label
_atom_site_scat_type_symbol
_atom_site_scat_Lande
_atom_site_scat_kappa
  Ho1  Ho3+  1.25  1.0
```

- `label` is the label of the magnetic ion (as in `atom_site_susceptibility`
  or `atom_site_moment`);
- `type_symbol` selects the <j0>, <j2> coefficients from the CrysPy table;
- `Lande` (2 by default) and `kappa` (1 by default) can be refined.

Ions missing from the table can be given by the `atom_type_scat` loop with the
coefficients `neutron_magnetic_j0_a1 ... neutron_magnetic_j6_e`.
For an aspherical form factor see the multipole model.
