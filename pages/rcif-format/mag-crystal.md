---
layout: page
permalink: /rcif-format/mag-crystal/
---
[Back to RCIF format][rcif-format]

Data items in the MAG_CRYSTAL category record details about
magnetic crystal structure.

# Mandatory Items/Loops:
- [Cell][cell]
- [Atom Site][atom-site]
- [Atom Site Moment][atom-site-moment]
- [Space Group Symop Magn Operation][space-group-symop-magn-operation]

# Optinal Items/Loops:
- [Atom Type][atom-type]
- [Atom Site Aniso][atom-site-aniso]
- [Space Group Symop Magn Centering][space-group-symop-magn-centering]
- [Atom Site Scat][atom-site-scat]
- [Atom Type Scat][atom-type-scat]
- [Atom Local Axes][atom-local-axes]
- [Atom Electron Confiduration][atom-electron-confiduration]


Minimal description:

      data_5yOhtAoR
      
      _cell_length_a 12.42600
      _cell_length_b 4.14200
      _cell_length_c 9.62400
      _cell_angle_alpha 90.00
      _cell_angle_beta 90.00
      _cell_angle_gamma 90.00
      
      loop_
      _space_group_symop_magn_centering_xyz
      _space_group_symop_magn_centering_id
      x,y,z,+1               1  
      x+1/2,y+1/2,z+1/2,+1   2  
      
      loop_
      _atom_site_scat_label
      _atom_site_scat_Lande
      _atom_site_scat_kappa
      U1_1   2.0   1.0  
      U1_2   2.0   1.0  
      
      loop_
      _atom_site_moment_label
      _atom_site_moment_symmform
      _atom_site_moment_crystalaxis_x
      _atom_site_moment_crystalaxis_y
      _atom_site_moment_crystalaxis_z
      U1_1   0,0,mz   0.0   0.0   -0.630(54)  
      U1_2   0,0,mz   0.0   0.0   0.505(53)   
      
      loop_
      _space_group_symop_magn_operation_xyz
      _space_group_symop_magn_operation_id
      x,y,z,+1      1  
      -x,-y,z,+1    2  
      -x,-y,-z,+1   3  
      x,y,-z,+1     4  
      x,-y,-z,-1    5  
      -x,y,-z,-1    6  
      -x,y,z,-1     7  
      x,-y,z,-1     8  
      
      loop_
      _atom_site_label
      _atom_site_type_symbol
      _atom_site_fract_x
      _atom_site_fract_y
      _atom_site_fract_z
      _atom_site_occupancy
      U1_1    U3+   0.000000   0.000000   0.000000   1.0   
      U1_2    U3+   0.333330   0.000000   0.000000   1.0   
      Rh1_1   Rh    0.500000   0.000000   0.250000   0.04  
      Rh1_2   Rh    0.166670   0.000000   0.250000   0.04  
      Ru1_1   Ru    0.500000   0.000000   0.250000   0.96  
      Ru1_2   Ru    0.166670   0.000000   0.250000   0.96  
      Si1_1   Si    0.000000   0.000000   0.629000   1.0   
      Si1_2   Si    0.333330   0.000000   0.629000   1.0   



[rcif-format]: /cryspy/rcif-format

[cell]: /cryspy/rcif-format/cell
[atom-site]: /cryspy/rcif-format/atom-site
[atom-type]: /cryspy/rcif-format/atom-type
[atom-site-aniso]: /cryspy/rcif-format/atom-site-aniso
[atom-site-scat]: /cryspy/rcif-format/atom-site-scat
[atom-type-scat]: /cryspy/rcif-format/atom-type-scat
[atom-local-axes]: /cryspy/rcif-format/atom-local-axes
[atom-electron-confiduration]: /cryspy/rcif-format/atom-electron-confiduration
[atom-site-moment]: /cryspy/rcif-format/atom-site-moment
[space-group-symop-magn-operation]: /cryspy/rcif-format/space-group-symop-magn-operation
[space-group-symop-magn-centering]: /cryspy/rcif-format/space-group-symop-magn-centering
