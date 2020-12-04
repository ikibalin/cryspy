---
layout: page
permalink: /rcif-format/pd/
---
[Back to RCIF format][rcif-format]

Data items in the PD category record details about
powder polarized and unpolarized neutron diffraction
measurements (1d case).

# Mandatory Items/Loops:
- [Setup][setup]
- [Pd Instr Resolution][pd-instr-resolution]
- [Phase][phase]
- [Pd Background][pd-background]
- [Pd Meas][pd-meas]

# Optinal Items/Loops:
- [Diffrn Radiation][diffrn-radiation]
- [Chi2][chi2]
- [Range][range]
- [Extinction][extinction]
- [Pd Instr Reflex Asymmetry][pd-instr-reflex-asymmetry]
- [Texture][texture]
- [Exclude][exclude]
- [PdProc][pd-proc]
- [PdPeak][pd-peak]
- [RefineLs][refine-ls]
- [Refln][refln]
- [ReflnSusceptibility][refln-susceptibility]

Example:

      data_pnd

      _setup_wavelength     0.840
      _setup_field          1.000
      _setup_offset_2theta -0.385(12)

      _chi2_sum  True
      _chi2_diff False
      _chi2_up   False
      _chi2_down False

      _range_2theta_min     4.000
      _range_2theta_max    80.000

      loop_
      _pd_background_2theta
      _pd_background_intensity
       4.5 256.0
      40.0 158.0
      80.0  65.0

      loop_
      _exclude_2theta_min
      _exclude_2theta_max
      0.0 1.0

      _pd_instr_reflex_asymmetry_p1 0.0
      _pd_instr_reflex_asymmetry_p2 0.0
      _pd_instr_reflex_asymmetry_p3 0.0
      _pd_instr_reflex_asymmetry_p4 0.0

      _diffrn_radiation_polarization 0.0
      _diffrn_radiation_efficiency   1.0

      _pd_instr_resolution_u 16.9776
      _pd_instr_resolution_v -2.8357
      _pd_instr_resolution_w  0.5763
      _pd_instr_resolution_x  0.0
      _pd_instr_resolution_y  0.0

      loop_
      _phase_label
      _phase_scale
      _phase_igsize
      Fe3O4 0.02381 0.0

      loop_
      _pd_meas_2theta
      _pd_meas_intensity_up
      _pd_meas_intensity_up_sigma
      _pd_meas_intensity_down
      _pd_meas_intensity_down_sigma
      4.0 465.80 128.97 301.88 129.30
      4.2 323.78 118.22 206.06 120.00
      4.4 307.14 115.90 230.47 116.53


[rcif-format]: /cryspy/rcif-format

[setup]: /cryspy/rcif-format/setup
[pd-instr-resolution]: /cryspy/rcif-format/pd-instr-resolution
[phase]: /cryspy/rcif-format/phase
[pd-background]: /cryspy/rcif-format/pd-background
[pd-meas]: /cryspy/rcif-format/pd-meas
[diffrn-radiation]: /cryspy/rcif-format/diffrn-radiation
[chi2]: /cryspy/rcif-format/chi2
[range]: /cryspy/rcif-format/range
[extinction]: /cryspy/rcif-format/extinction
[pd-instr-reflex-asymmetry]: /cryspy/rcif-format/pd-instr-reflex-asymmetry
[texture]: /cryspy/rcif-format/texture
[exclude]: /cryspy/rcif-format/exclude
[pd-proc]: /cryspy/rcif-format/pd-proc
[pd-peak]: /cryspy/rcif-format/pd-peak
[refine-ls]: /cryspy/rcif-format/refine-ls
[refln]: /cryspy/rcif-format/refln
[refln-susceptibility]: /cryspy/rcif-format/refln-susceptibility
