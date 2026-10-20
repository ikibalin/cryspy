---
title: Ordered magnetic structure (magnetic space group, mcif)
keywords: magnetic structure, ordered moment, moment, mcif, magnetic space group, shubnikov, antiferromagnet, ferromagnet, symmform, магнитная структура, магнитный момент, магнитная группа
categories: atom_site_moment, space_group_symop_magn_operation, space_group_symop_magn_centering
---
An ordered (commensurate) magnetic structure is described in the crystal data
block in the mcif way:

- instead of `space_group`, give the magnetic space group by the loops
  `space_group_symop_magn_operation` (`xyz` such as `-y,x-y+1/2,z,-1`) and
  `space_group_symop_magn_centering`;
- `atom_site_moment` gives the moment of each magnetic ion: `crystalaxis_x`,
  `crystalaxis_y`, `crystalaxis_z` in μB along the crystal axes (refinable) and
  `symmform` (for example `0,0,mz` or `mx,mx,0`) with the symmetry constraints;
- `atom_site_scat` gives the magnetic form factor.

Magnetic structures with a propagation vector k (incommensurate, or
commensurate described by k) are not supported yet.
