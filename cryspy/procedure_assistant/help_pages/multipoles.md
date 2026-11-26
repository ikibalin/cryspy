---
title: Multipole model of the magnetization (form factor)
keywords: multipole, multipoles, add multipoles, check requirements, propagation vector, superspace, aspherical, anisotropic form factor, pseudoatom, hansen-coppens, kappa, plm, orbital, local axes, dummy atom, core, valence, electron configuration, open shell, what is needed
categories: atom_rho_multipole, atom_local_axes, atom_site, atom_electron_configuration
---
The anisotropic magnetic form factor of an atom can be described by the
multipole (Hansen-Coppens) model. Two loops of the crystal data block are
needed: `atom_rho_multipole` with the coefficients and `atom_local_axes`
with the local coordinate system of each atom.

What a magnetic atom needs for its form factor through multipoles:

1. a row of `atom_site` (the atom itself) and of `atom_site_scat` (its
   magnetic ion);
2. local axes: a row of `atom_local_axes`. The coefficients `coeff_Plm`
   are given in the local coordinate system of the atom, so without the
   axes they have no meaning. The axes are best aligned with the bonds
   or the symmetry axes of the site; dummy atoms give directions without
   an atom;
3. the radial part: the ion (`atom_type_symbol`) and its open shell
   (`orbital`: 3d for Mn2+, 4f for rare earths, 2p for a spin transferred
   to O2- or F1-), with `kappa`;
4. the populations `coeff_P00` (1 for a spherical start) and the other
   `coeff_Plm` allowed by the site symmetry.

The core electrons (the closed shells, loop `atom_electron_configuration`
with `core` and `valence`) are not needed for the magnetic form factor:
their spins are paired, only the open shell given by `orbital` carries
the magnetization, and the polarisation of the core is neglected. The
core part is used for the electron (charge) density and for the core
prior density of MEM (`_mem_parameters_prior_density core`); give it
only there.

`add_item {"category": "atom_rho_multipole"}` adds the missing rows of
items 1 and 2 together with the multipoles and asks the values it cannot
find (the local axes are proposed from the nearest atoms).

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

The same rows describe the form factor of magnetic space groups
(`space_group_symop_magn_operation`, `atom_site_moment`) and of magnetic
superspace groups with a propagation vector
(`space_group_symop_magn_ssg_operation`, `atom_site_moment_Fourier`): the
density of every generated atom is rotated by its symmetry operation.

The method `Crystal.add_multipole_form_factor` adds both rows for one or
several atoms in one step (AI assistant: `call_method` on the crystal):

```
crystal_cula2o4.add_multipole_form_factor(atom_label="Cu", coefficients="P20 0.1")
crystal_mnf2.add_multipole_form_factor(atom_label="Mn, F", magnetic_ion="F: O2-", orbital="F: 2p")
```

It needs the atoms in `atom_site` with their magnetic ion in
`atom_site_scat`; `check_requirements {"name": "add_multipole_form_factor"}`
shows which atoms can get multipoles, which lack `atom_site_scat` (dummy
atoms) and which have them already. It takes the ion from `atom_site_scat`
or `atom_site` (Cu2+), the orbital from the open shell of the ion (3d),
sets P00 = 1, and adds the local axes z || c and x || a* with new dummy
atoms; a dummy atom put near the atom before gives one of the axes. When
the ion of `atom_site_scat` is of another element than the atom (O2- for
F) or is not a magnetic ion (spin transferred to F1-, O2-: orbital 2p),
the ion and the orbital are asked with the possible values. Atoms with
multipoles already are not changed. Change `atom_local_axes` to align the
axes with the bonds.

Rows can also be added one by one: to add the multipoles of an atom to a loaded object, add a row to each of
the two loops (the AI assistant does it with the action `add_item`; a
loop is created when the crystal has none):

```
add_item {"category": "atom_local_axes", "values": {"atom_label": "Fe", "atom0": "D1", "ax1": "Z", "atom1": "Fe", "atom2": "D2", "ax2": "X"}}
add_item {"category": "atom_rho_multipole", "values": {"atom_label": "Fe", "atom_type_symbol": "Fe3+", "orbital": "3d"}}
```

`remove_item {"path": "crystal_fe.atom_rho_multipole[Fe]"}` removes the row again.
