---
title: 3D view of the crystal structure
keywords: 3d, view, ellipsoid, ellipsoids, adp, negative biso, negative uiso, negative adp, non-positive definite, npd, cubes, thermal ellipsoids, displacement, susceptibility ellipsoid, magnetization ellipsoid, p1, unit cell atoms, dummy, visualization, visualisation, structure, picture, draw, show, py3dmol, 3dmol, moments, arrows, unit cells, bonds, jupyter
categories: atom_site, atom_site_aniso, atom_site_moment, atom_site_susceptibility, cell, space_group
---
A crystal data block can be shown as an interactive 3D picture: atoms of the
chosen number of unit cells (generated from the asymmetric unit with the space
group or the magnetic space group of mcif), bonds, the cell edges, the a, b, c
axes and the magnetic moments of `atom_site_moment` as arrows.

The view needs the optional extra (py3Dmol and ipywidgets):

```
pip install cryspy[view3d]
```

In a Jupyter notebook:

```python
import cryspy
rcif = cryspy.load_file("main.rcif")
crystal = rcif.items[0]                 # a Crystal data block
crystal.view_structure()                # one unit cell
crystal.view_structure(n_cells_a=2, n_cells_b=2, show_labels=True)
crystal.view_structure(interactive=True)  # sliders and check boxes
```

All settings are fields of `cryspy.StructureViewSettings`: number of cells
along a, b, c, how atoms are drawn, sphere size, bonds, labels, cell edges,
axes, atoms on the far faces, dummy atoms (hidden by default), moments, the
moment scale, colours and the size of the view.

`atom_display` chooses how atoms are drawn:

- `spheres`: Jmol colours, `atom_scale` times the van der Waals radius;
- `ADP ellipsoids`: displacement ellipsoids of `Uiso`/`Biso` (spheres) and
  `Uani`/`Bani` atoms holding the atom with probability `adp_probability`
  (0.5 by default, as ORTEP);
- `susceptibility ellipsoids`: magnetization ellipsoids of the atoms of
  `atom_site_susceptibility`, the moments M = χB induced by a field of
  `susceptibility_field` tesla in every direction (semi-axes |χ_i|·B).

Atoms with non-positive definite displacement parameters (a negative
`Biso`/`Uiso`, or a `Uani` with a negative eigenvalue U_i) are not physical
but are still drawn, with |U|, and marked as `npd_display` chooses:

- `cubes` (default): boxes along the principal axes with half edges of the
  |U| ellipsoid, the usual convention of ORTEP, SHELX and Olex2;
- `ellipsoids with bars`: half transparent ellipsoids of |U| crossed by
  magenta bars along the axes with a negative U_i, so an anisotropic atom
  shows in which direction its displacement is wrong (a negative `Biso`
  gives three bars).

A mark is at least 0.4 Å in size, so a small negative value is still seen.
A warning names these atoms, and `report_ellipsoids` lists their U_i with
the usual causes (wrong atom type or occupancy, missing absorption or
extinction correction, scale, background, correlations).

Atoms without such parameters stay spheres. `moment_scale` (Å per μB) gives
the length of the moment arrows and the size of the magnetization
ellipsoids, so both are on one scale.

Rotate with the left mouse button, zoom with the wheel, move with the right
button.

`crystal.get_structure_view_html()` gives an HTML page of the view, used by
CrysPy Editor. Pass `js=` with the path of a local `3Dmol-min.js` to work
without internet.

Related assistant tools: `list_atoms_in_unit_cell` (all atoms of the cell
with their moments), `convert_to_p1` (a copy of the crystal in P1) and
`report_ellipsoids` (principal axes of the displacement or magnetization
ellipsoids, also for every atom of the cell). In Python:
`crystal.calc_atoms_in_unit_cell()`, `crystal.to_p1()`,
`crystal.report_ellipsoids("ADP")`.
