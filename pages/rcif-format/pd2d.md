---
layout: page
permalink: /rcif-format/pd2d/
---
[Back to RCIF format][rcif-format]

Data items in the PD2D category record details about
powder polarized neutron diffraction
measurements (2d case).

# Mandatory Items/Loops:
- [Setup][setup]
- [Pd2d Instr Resolution][pd2d-instr-resolution]
- [Phase][phase]
- [Pd2d Background][pd2d-background]
- [Pd2d Meas][pd2d-meas]
- [Diffrn Radiation][diffrn-radiation]
- [Range][range]
- [Chi2][chi2]

# Optional Items/Loops:
- [Extinction][extinction]
- [Pd2d Instr Reflex Asymmetry][pd2d-instr-reflex-asymmetry]
- [Texture][texture]
- [Exclude][exclude]
- [Pd2dProc][pd2d-proc]
- [Pd2dPeak][pd2d-peak]
- [RefineLs][refine-ls]
- [Refln][refln]
- [ReflnSusceptibility][refln-susceptibility]

Example:

      data_powder2d

      _setup_wavelength     0.840
      _setup_field          1.000
      _setup_offset_2theta -0.385

      _chi2_sum  True
      _chi2_diff False
      _chi2_up   False
      _chi2_down False

      _range_2theta_min     4.000
      _range_2theta_max    80.000
      _range_phi_min       -2.000
      _range_phi_max       40.000

      _pd2d_background_2theta_phi_intensity
      ;
      2 5 80
      5 2. 2.(1)
      40 2. 2.
      ;

      loop_
      _exclude_2theta_min
      _exclude_2theta_max
      _exclude_phi_min
      _exclude_phi_max
      0.0 1.0 0. 1.0

      _pd2d_instr_reflex_asymmetry_p1 0.0
      _pd2d_instr_reflex_asymmetry_p2 0.0
      _pd2d_instr_reflex_asymmetry_p3 0.0
      _pd2d_instr_reflex_asymmetry_p4 0.0

      _diffrn_radiation_polarization 0.0
      _diffrn_radiation_efficiency   1.0

      _pd2d_instr_resolution_u 16.9776
      _pd2d_instr_resolution_v -2.8357(1)
      _pd2d_instr_resolution_w  0.5763
      _pd2d_instr_resolution_x  0.0
      _pd2d_instr_resolution_y  0.0

      loop_
      _phase_label
      _phase_scale
      _phase_igsize
      Fe3O4 0.02381() 0.0

      _pd2d_meas_2theta_phi_intensity_up
      ;
      2 5 80
      5 2. 2.
      40 2. 2.
      ;

      _pd2d_meas_2theta_phi_intensity_up_sigma
      ;
      2 5 80
      5 2. 2.
      40 2. 2.
      ;

      _pd2d_meas_2theta_phi_intensity_down
      ;
      2 5 80
      5 2. 2.
      40 2. 2.
      ;

      _pd2d_meas_2theta_phi_intensity_down_sigma
      ;
      2 5 80
      5 2. 2.
      40 2. 2.
      ;


[rcif-format]: /cryspy/rcif-format

[setup]: /cryspy/rcif-format/setup
[pd2d-instr-resolution]: /cryspy/rcif-format/pd2d-instr-resolution
[phase]: /cryspy/rcif-format/phase
[pd2d-background]: /cryspy/rcif-format/pd2d-background
[pd2d-meas]: /cryspy/rcif-format/pd2d-meas
[diffrn-radiation]: /cryspy/rcif-format/diffrn-radiation
[chi2]: /cryspy/rcif-format/chi2
[range]: /cryspy/rcif-format/range
[extinction]: /cryspy/rcif-format/extinction
[pd2d-instr-reflex-asymmetry]: /cryspy/rcif-format/pd2d-instr-reflex-asymmetry
[texture]: /cryspy/rcif-format/texture
[exclude]: /cryspy/rcif-format/exclude
[pd2d-proc]: /cryspy/rcif-format/pd2d-proc
[pd2d-peak]: /cryspy/rcif-format/pd2d-peak
[refine-ls]: /cryspy/rcif-format/refine-ls
[refln]: /cryspy/rcif-format/refln
[refln-susceptibility]: /cryspy/rcif-format/refln-susceptibility
