---
title: Interatomic distances, bond angles, torsion angles, hydrogen bonds
keywords: distance, distances, interatomic distance, bond, bonds, bond length, bond lengths, angle, angles, bond angle, torsion, torsion angle, dihedral, neighbours, neighbors, coordination, environment, hydrogen bond, hydrogen bonds, hbond, geom_bond, geom_angle, geom_torsion, geom_hbond, publication cif, symmetry code
categories: atom_site, cell, space_group, geom_bond, geom_angle, geom_torsion, geom_hbond
---
The geometry of a crystal is calculated from the cell, the atoms of
`atom_site` and the symmetry elements (space group, magnetic or superspace
group). Every atom is an atom of `atom_site` moved by a symmetry operation and
a lattice translation; the operation is written as '-x+1,y+1/2,-z'.

Assistant tool `report_distances_and_angles`:

- `atoms: "N1 C1"` – the distance N1–C1 (C1 is its equivalent position
  closest to N1) and all other N1–C1 distances up to `max_distance`;
- `atoms: "O1 Fe1 O2"` – the angle at the middle atom Fe1;
- `atoms: "N1 C3 C4 H4A"` – the torsion angle (IUPAC sign);
- `atoms: "Fe1"` – all neighbours of Fe1 up to `max_distance` (3 Å by
  default) and the angles between these bonds;
- no atoms – all bonds and angles of the structure;
- `cif: true` – the same as CIF loops geom_bond, geom_angle,
  geom_torsion with the loop space_group_symop that their symmetry
  codes 'n_klm' refer to;
- `hydrogen_bonds: true` – the hydrogen bonds D–H···A.

Standard uncertainties come from the sigmas of the refined coordinates and
of the cell parameters, taken as independent (no correlations).

In Python:

```python
crystal.calc_distance_or_angle("N1 C1")       # dict: value, sigma, symmetry
print(crystal.report_distances_and_angles("Fe1", max_distance=2.5))
print(crystal.report_hydrogen_bonds(acceptors="N O", max_distance_ha=2.6))
crystal.add_distances_and_angles(max_distance=2.0, torsions="N1 C3 C4 H4A")
```

`add_distances_and_angles` stores the loops `geom_bond`, `geom_angle`,
`geom_torsion` and `geom_hbond` in the crystal: they are saved in the rcif
file and can be copied into a publication CIF. They are results, not
parameters: call the method again after a refinement, it replaces them.

Hydrogen bonds: for every H (or D) atom the donor D is the closest other atom
within `max_distance_dh` (1.2 Å), the acceptors are atoms of the elements
`acceptors` ('N O F Cl S') within `max_distance_ha` (2.6 Å) of H with the angle
D–H···A of at least `min_angle_dha` (120°).

To see the atoms themselves, `list_atoms_in_unit_cell` takes `atoms` too
(e.g. 'N1 C1') and lists only their equivalent positions in the unit cell.
