---
layout: page
permalink: /rcif-format/space-group-symop/
---
[Back to RCIF format][rcif-format]

# `space_group_symop` item/loop:

     loop_
     _space_group_symop.id
     _space_group_symop.operation_xyz
     _space_group_symop.operation_description     
     1    x,y,z              'identity mapping'
     2    -x,-y,-z           'inversion'
     3    -x,1/2+y,1/2-z '2-fold screw rotation with axis in (0,y,1/4)'
     4    x,1/2-y,1/2+z  'c glide reflection through the plane (x,1/4,y)'

**Mandatory parameters:** 
`operation_xyz`.

**Optional parameters:** 
`id`, `operation_description`, `generator_xyz`, `sg_id`.

# Item object

**Methods:** 
`.define_by_el_symm(...)`, `.get_symop_inversed(...)`, `.get_symops_by_centring_type(...)`,
`.get_symops_by_generator_xyz(...)`, `.get_coords_xyz_by_coord_xyz(...)`.

**Internal parameters:** 
`r`, `b`, `r_11`, `r_12`, `r_13`, `r_21`, `r_22`,
`r_23`, `r_31`, `r_32`, `r_33`, `b_1`, `b_2`, `b_3`.

# Loop object
**Methods:** 
`.create_by_generators_xyz(...)`, `.get_coords_xyz_by_coord_xyz(...)`,
`.get_symop_for_x1_x2(...)`.

[rcif-format]: /cryspy/rcif-format
