---
layout: page
permalink: /rcif-format/diffrn/
---
[Back to RCIF format][rcif-format]

Data items in the DIFFRN category record details about
single diffraction measurements.

# Mandatory Items/Loops:
- [Setup][setup]
- [Diffrn Radiation][diffrn-radiation]
- [Diffrn Orient Matrix][diffrn-orient-matrix]
- [Diffrn Refln][diffrn-refln]

# Optinal Items/Loops:
- [Extinction][extinction]
- [Phase][phase]
- [Refln][refln]
- [Refln Susceptibility][refln-susceptibility]
- [RefineLs][refine-ls]

Example:

      data_mono

      _setup_wavelength     0.840
      _setup_field          1.000

      _diffrn_radiation_polarization 1.0
      _diffrn_radiation_efficiency   1.0

      _extinction_mosaicity 100.0
      _extinction_radius    50.0
      _extinction_model     gauss

      _diffrn_orient_matrix_UB_11 6.59783
      _diffrn_orient_matrix_UB_12 -6.99807
      _diffrn_orient_matrix_UB_13 3.3663
      _diffrn_orient_matrix_UB_21 2.18396
      _diffrn_orient_matrix_UB_22 -2.60871
      _diffrn_orient_matrix_UB_23 -9.5302
      _diffrn_orient_matrix_UB_31 7.4657
      _diffrn_orient_matrix_UB_32 6.94702
      _diffrn_orient_matrix_UB_33 -0.18685

      _phase_label  Fe3O4

      loop_
      _diffrn_refln_index_h
      _diffrn_refln_index_k
      _diffrn_refln_index_l
      _diffrn_refln_fr
      _diffrn_refln_fr_sigma
      0 0 8 0.64545 0.01329
      2 0 6 1.75682 0.04540
      0 2 6 1.67974 0.03711


[rcif-format]: /cryspy/rcif-format

[setup]: /cryspy/rcif-format/setup
[diffrn-radiation]: /cryspy/rcif-format/diffrn-radiation
[diffrn-orient-matrix]: /cryspy/rcif-format/diffrn-orient-matrix
[diffrn-refln]: /cryspy/rcif-format/diffrn-refln
[extinction]: /cryspy/rcif-format/extinction
[phase]: /cryspy/rcif-format/phase
[refln]: /cryspy/rcif-format/refln
[refln-susceptibility]: /cryspy/rcif-format/refln-susceptibility
[refine-ls]: /cryspy/rcif-format/refine-ls
