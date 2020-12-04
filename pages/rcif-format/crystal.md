---
layout: page
permalink: /rcif-format/crystal/
---
[Back to RCIF format][rcif-format]

Data items in the CRYSTAL category record details about
crystal structure.

# Mandatory Items/Loops:
- [Space Group][space-group]
- [Cell][cell]
- [Atom Site][atom-site]

# Optional Items/Loops:
- [Atom Type][atom-type]
- [Atom Site Aniso][atom-site-aniso]
- [Atom Site Susceptibility][atom-site-susceptibility]
- [Atom Site Scat][atom-site-scat]
- [Atom Type Scat][atom-type-scat]
- [Atom Local Axes][atom-local-axes]
- [Atom Electron Confiduration][atom-electron-confiduration]

Minimal description:

      data_Fe3O4

      _cell_angle_alpha 90.0
      _cell_angle_beta 90.0
      _cell_angle_gamma 90.0
      _cell_length_a 8.56212()
      _cell_length_b 8.56212
      _cell_length_c 8.56212

      _space_group_IT_number    227
      _space_group_it_coordinate_system_code 2


      loop_
      _atom_site_adp_type
      _atom_site_B_iso_or_equiv
      _atom_site_fract_x
      _atom_site_fract_y
      _atom_site_fract_z
      _atom_site_label
      _atom_site_occupancy
      _atom_site_type_symbol
      Uiso 0.0 0.12500 0.12500 0.12500 Fe3A 1.0 Fe3+
      Uiso 0.0 0.50000 0.50000 0.50000 Fe3B 1.0 Fe3+
      Uiso 0.0 0.25521 0.25521 0.25521   O1 1.0  O2-


[rcif-format]: /cryspy/rcif-format

[space-group]: /cryspy/rcif-format/space-group
[cell]: /cryspy/rcif-format/cell
[atom-site]: /cryspy/rcif-format/atom-site
[atom-type]: /cryspy/rcif-format/atom-type
[atom-site-aniso]: /cryspy/rcif-format/atom-site-aniso
[atom-site-susceptibility]: /cryspy/rcif-format/atom-site-susceptibility
[atom-site-scat]: /cryspy/rcif-format/atom-site-scat
[atom-type-scat]: /cryspy/rcif-format/atom-type-scat
[atom-local-axes]: /cryspy/rcif-format/atom-local-axes
[atom-electron-confiduration]: /cryspy/rcif-format/atom-electron-confiduration
