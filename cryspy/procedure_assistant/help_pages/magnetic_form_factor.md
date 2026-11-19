---
title: Magnetic form factor in the dipole approximation
keywords: form factor, magnetic form factor, magnetic ion, valence, usual ion, knowledge base, free-ion moment, lande, kappa, j0, j2, dipole approximation, atom_site_scat, атомный форм-фактор, магнитный форм-фактор, фактор ланде
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

## Knowledge base of the magnetic ions

The usual ions of the magnetic elements, their configurations, Lande
factors and free-ion moments are in the knowledge base of CrysPy
(`cryspy/knowledge/magnetic_ions_4f.md`, `magnetic_ions_3d.md`,
`magnetic_ions_4d_5d_5f.md`; the usual charges of the other ions for the
charge balance in `ion_charges.md`). CrysPy reads the Lande factor of an ion
from these tables (`lande` 0 in `add_propagation_vector`), and an AI
assistant proposes from them the magnetic atoms and the magnetic ion of a
crystal (Tb3+ in TbTe3; Fe3+ in Fe2O3 from the charge balance) for the user
to confirm. The tool `search_knowledge` searches these files (the same
search as for the documents of the user).
