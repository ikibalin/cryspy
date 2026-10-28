---
title: 3D view of the crystal structure
keywords: 3d, view, visualization, visualisation, structure, picture, draw, show, py3dmol, 3dmol, moments, arrows, unit cells, bonds, jupyter, визуализация, 3d вид, структура, картинка, нарисовать, показать, моменты, стрелки, ячейки, связи
categories: atom_site, atom_site_moment, cell, space_group
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
along a, b, c, atom size, bonds, labels, cell edges, axes, atoms on the far
faces, moments, arrow length of the largest moment (Å), colours and the size
of the view. Rotate with the left mouse button, zoom with the wheel, move with
the right button.

`crystal.get_structure_view_html()` gives an HTML page of the view, used by
CrysPy Editor. Pass `js=` with the path of a local `3Dmol-min.js` to work
without internet.
