---
layout: page
permalink: /rcif-format/density-point/
---
[Back to RCIF format][rcif-format]

# `density_point` item/loop:

      loop_Ho1
      _density_point_index_x
      _density_point_index_y
      _density_point_index_z
      _density_point_density
      _density_point_multiplicity
      0  0  0  0.32  2
      0  0  1  0.31  1
      0  0  2  0.29  1

**Mandatory parameters:** 
`index_x`, `index_y`, `index_z`.

**Optional parameters:** 
`density`, `density_ferro`, `density_antiferro`,
`multiplicity`, `fract_x`, `fract_y`, `fract_z`,
`basin_atom_label`, `basin_atom_symop`,
`basin_atom_distance`.

# Item object

**Methods:** 
`.calc_multiplicity(...)`.

# Loop object
**Methods:** 
`.form_asymmetric_unit_cell(...)`, `.separate_basins(...)`,
`.set_flat_density(...)`, `.set_core_density(...)`,
`.create_flat_density(...)`, `.create_core_density(...)`,
`.calc_density_chi(...)`, `.save_to_file_den(...)`,
`.calc_fm(...)`, `.calc_phase_3d(....)`, `.calc_rbs_i(...)`,
`.calc_susc_i(...)`, `.calc_moment_2d(...)`,
`.calc_factor_in_front_of_density_for_fm(...)`,
`.calc_factor_in_front_of_density_for_fm_perp(...)`,
`.calc_electrons_per_unit_cell(...)`, `.renormalize_numpy_densities(...)`.

[rcif-format]: /cryspy/rcif-format
