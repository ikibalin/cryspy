---
layout: page
permalink: /rcif-format/diffrn-orient-matrix/
---
[Back to RCIF format][rcif-format]

# `difrrn_orient_matrix` item/loop:

     _diffrn_orient_matrix_UB_11           -0.04170
     _diffrn_orient_matrix_UB_12           -0.01429
     _diffrn_orient_matrix_UB_13           -0.02226
     _diffrn_orient_matrix_UB_21           -0.00380
     _diffrn_orient_matrix_UB_22           -0.05578
     _diffrn_orient_matrix_UB_23           -0.05048
     _diffrn_orient_matrix_UB_31            0.00587
     _diffrn_orient_matrix_UB_32           -0.13766
     _diffrn_orient_matrix_UB_33            0.02277

**Mandatory parameters:** 
`ub_11`, `ub_12`, `ub_13`, `ub_21`, `ub_22`,
`ub_23`, `ub_31`, `ub_32`, `ub_33`.

**Optional parameters:** 
`id`, `type`.

Constraint of `type` parameters is "CCSL", "6T2@LLB", "5C1@LLB", "BusingLevy",
"-YXZ", "X-YZ", "XYZ"

# Item object

**Methods:** 
`.calc_ub_ccsl(...)`, `.calc_matrix_from_ccsl(...)`, `.calc_ub(...)`,
`.calc_e_up(...)`, `.calc_angle(...)`, `.calc_q2(...)`

**Internal parameters:** 
`u_11`, `u_12`, `u_13`, `u_21`, `u_22`, `u_23`,
`u_31`, `u_32`, `u_33`, `ub`, `u` 

**Protected parameters:** 
`cell`

# Loop object
**Methods:** 
No

[rcif-format]: /cryspy/rcif-format