---
layout: page
permalink: /rcif-format/space-group/
---
[Back to RCIF format][rcif-format]

# `space_group` item/loop:

     _space_group_IT_number    227
     _space_group_it_coordinate_system_code 2

or 

     _space_group_name_H-M_alt "F d -3 m"


**Mandatory parameters:** 
No

**Optional parameters:** 
`id`, `name_H-M_alt`, `name_H-M_alt_description`, `IT_number`,
`name_H-M_ref`, `IT_coordinate_system_code`, `name_H-M_full`,
`name_Hall`, `name_Schoenflies`, `point_group_H-M`, `Laue_class`,
`Patterson_name_H-M`, `centring_type`, `Bravais_type`,
`crystal_system`, `reference_setting`, `transform_pp_abc`,
`transform_qq_xyz`.

# Item object

**Methods:** 
`.calc_hkl_equiv(...)`, `.calc_xyz_mult(...)`, `.calc_symop_for_xyz(...)`,
`.calc_el_symm_for_xyz(...)`, `.calc_f_hkl_by_f_hkl_as(...)`,
`.calc_asymmetric_cell(...)`, `.calc_rotated_matrix_for_position(...)`,
`.report_space_group(...)`.

**Internal parameters:** 
`centrosymmetry`, `pcentr`, `reduced_space_group_symop`,
`full_space_group_symop`, `space_group_wyckoff`, `shift`.


# Loop object
**Methods:** 
No

[rcif-format]: /cryspy/rcif-format