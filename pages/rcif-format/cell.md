---
layout: page
permalink: /rcif-format/cell/
---
[Back to RCIF format][rcif-format]

# `cell` item/loop:

     _cell_length_a 8.0
     _cell_length_b 3.76(5)
     _cell_length_c 4.0
     _cell_angle_alpha 35.0
     _cell_angle_beta 35.0
     _cell_angle_gamma 40.0

Record details about the crystallographic cell parameters and their measurement.

**Mandatory parameters:** 
`length_a`: angstroms, `length_b`: angstroms, `length_c`: angstroms,
`angle_alpha`: degrees, `angle_beta`: degrees, `angle_gamma`: degrees.

# Item object

**Methods:** 
`.calc_sthovl(...)`,
`.calc_k_loc(...)`,
`.calc_m_t(...)`,
`.calc_hkl(...)`,
`.calc_hkl_in_range(...)`,
`.calc_position_by_coordinate(...)`,
`.calc_coordinate_by_position(...)`,
`.calc_length_sq(...)`,
`.closest_distance_between_fractions(...)`,
`.calc_reciprocal_length_sq(...)`,
`.ortogonalize_matrix(...)`,
`.report(...)`

**Internal parameters:** 
`m_g`, `m_g_reciprocal`, `m_b`, `m_b_norm`, `m_m`, `m_m_norm`,
`cos_a`, `cos_b`, `cos_g`, `sin_a`, `sin_b`, `sin_g`, `cos_a_sq`,
`cos_b_sq`, `cos_g_sq`, `sin_a_sq`, `sin_b_sq`, `sin_g_sq`, `cos_ia`,
`cos_ib`, `cos_ig`, `sin_ia`, `sin_ib`, `sin_ig`, `cos_ia_sq`,
`cos_ib_sq`, `cos_ig_sq`, `sin_ia_sq`, `sin_ib_sq`, `sin_ig_sq`,
`reciprocal_length_a`, `reciprocal_length_b`, `reciprocal_length_c`,
`reciprocal_angle_alpha`, `reciprocal_angle_beta`,
`reciprocal_angle_gamma`, `volume`

**Protected parameters:** 
`type_cell`, `it_coordinate_system_code`,
`formula_units_z`

# Loop object
**Methods:** 
No


[rcif-format]: /cryspy/rcif-format