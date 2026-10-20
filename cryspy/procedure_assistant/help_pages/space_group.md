---
title: Space group and unit cell
keywords: space group, symmetry, cell, unit cell, lattice, it number, hermann-mauguin, setting, wyckoff, пространственная группа, симметрия, элементарная ячейка, параметры ячейки
categories: space_group, cell, space_group_symop, space_group_wyckoff
---
The crystal data block needs the `cell` item: `length_a`, `length_b`,
`length_c` in Å and `angle_alpha`, `angle_beta`, `angle_gamma` in degrees.
When the space group is given, the cell parameters fixed by the crystal
system follow the free ones and are not refined separately.

The space group is the `space_group` item. Give one of `name_H-M_alt`
(for example `"P n m a"`), `IT_number`, `name_Hall` and, for a non-standard
setting, `IT_coordinate_system_code` (for example `abc` for the axes setting
or `2` for the second origin choice).

```
_space_group_name_H-M_alt "F d -3 m"
_space_group_IT_coordinate_system_code 2
```

Atomic positions on special sites are constrained by the site symmetry.
For magnetic structures the magnetic space group is given instead (see the
magnetic structure page).
