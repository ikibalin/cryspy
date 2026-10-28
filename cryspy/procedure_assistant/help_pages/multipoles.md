---
title: Multipole model of the magnetization (form factor)
keywords: multipole, multipoles, aspherical, anisotropic form factor, pseudoatom, hansen-coppens, kappa, plm, orbital, local axes, dummy atom, мультиполи, мультипольная модель, анизотропный форм-фактор, локальные оси
categories: atom_rho_multipole, atom_local_axes, atom_site, atom_electron_configuration
---
The anisotropic magnetic form factor of an atom can be described by the
multipole (Hansen-Coppens) model. Two loops of the crystal data block are
needed: `atom_rho_multipole` with the coefficients and `atom_local_axes`
with the local coordinate system of each atom.

`atom_rho_multipole` (one row per atom):

- `atom_label` is the label from `atom_site`;
- `atom_type_symbol` is the ion (for example `Mn2+`), with `orbital`
  (for example `3d`) it selects the Slater radial function from the CrysPy
  database;
- `kappa` scales the radial function (1 by default);
- `coeff_P00`, `coeff_P11`, `coeff_P1-1`, `coeff_P10`, ... `coeff_P4-4` are the
  population coefficients of the real spherical harmonics up to l = 4.
  Give only the coefficients allowed by the site symmetry: CrysPy does not
  check symmetry constraints of multipole parameters.

`atom_local_axes` (one row per atom): the first local axis `ax1` (X, Y or Z)
goes from `atom_label` to `atom0`, the second axis `ax2` is defined by the
vector from `atom1` to `atom2`. Dummy atoms (`_atom_site_type_symbol Dummy`)
can be added to `atom_site` to define directions that do not point to a
real atom; they do not scatter.

```
loop_
_atom_local_axes_atom_label
_atom_local_axes_atom0
_atom_local_axes_ax1
_atom_local_axes_atom1
_atom_local_axes_atom2
_atom_local_axes_ax2
  Mn  D1  X  Mn  D2  Z

loop_
_atom_rho_multipole_atom_label
_atom_rho_multipole_atom_type_symbol
_atom_rho_multipole_orbital
_atom_rho_multipole_kappa
_atom_rho_multipole_coeff_P00
_atom_rho_multipole_coeff_P20
_atom_rho_multipole_coeff_P22
  Mn  Mn2+  3d  1.0()  1.0  0.0()  0.0()
```

`kappa` and the `coeff_P...` values can be refined against flip ratios.
The model is also used as a prior density for MEM (`_mem_parameters_prior_density core`).
For a spherical form factor in the dipole approximation use `atom_site_scat` instead.

To add the multipoles of an atom to a loaded object, add a row to each of
the two loops (the AI assistant does it with the action `add_item`; a
loop is created when the crystal has none):

```
add_item {"category": "atom_local_axes", "values": {"atom_label": "Fe", "atom0": "D1", "ax1": "Z", "atom1": "Fe", "atom2": "D2", "ax2": "X"}}
add_item {"category": "atom_rho_multipole", "values": {"atom_label": "Fe", "atom_type_symbol": "Fe3+", "orbital": "3d"}}
```

`remove_item {"path": "crystal_fe.atom_rho_multipole[Fe]"}` removes the row again.
