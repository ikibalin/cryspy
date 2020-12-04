---
layout: page
permalink: /rcif-format/tof/
---
[Back to RCIF format][rcif-format]

Data items in the TOF category record details about
powder diffraction measurements by Time of FLight
(1f case).

# Mandatory Items/Loops:
- [TOF Parameters][tof-parameters]
- [TOF Profile][tof-profile]
- [Phase][phase]
- [TOF Background][tof-background]
- [TOF Meas][tof-meas]

# Optinal Items/Loops:
- [Diffrn Radiation][diffrn-radiation]
- [Chi2][chi2]
- [Range][range]
- [Extinction][extinction]
- [Texture][texture]
- [Exclude][exclude]
- [TOF Proc][tof-proc]
- [TOF Peak][tof-peak]
- [RefineLs][refine-ls]
- [Refln][refln]
- [ReflnSusceptibility][refln-susceptibility]

Example:

      data_tof

      _range_time_min 3001.589
      _range_time_max 19000.0

      _tof_parameters_Zero 2.92100
      _tof_parameters_Dtt1 6167.24700
      _tof_parameters_2theta_bank 145.00
      _tof_parameters_neutrons thermal
      _tof_parameters_dtt2 -2.28000

      _tof_profile_alpha0 0.00000
      _tof_profile_alpha1 0.29710
      _tof_profile_beta0 0.04182
      _tof_profile_beta1 0.00224
      _tof_profile_sigma0 0.409(35)
      _tof_profile_sigma1 8.118(30)
      _tof_profile_peak_shape Gauss

      _tof_background_time_max 19000.00
      _tof_background_coeff1 24832.85000
      _tof_background_coeff2 6139.24400
      _tof_background_coeff3 8063.47200
      _tof_background_coeff4 3125.05000
      _tof_background_coeff5 2566.95600
      _tof_background_coeff6 311.07700
      _tof_background_coeff7 837.34800
      _tof_background_coeff8 -103.74200
      _tof_background_coeff9 -11.80600

      loop_
      _phase_label
      _phase_scale
      _phase_igsize
      cecual   1359.60(56)   0.0  

      loop_
      _tof_meas_time
      _tof_meas_intensity
      _tof_meas_intensity_sigma
      3001.589    40409.00    462.00    
      3003.090    40171.00    460.00    
      3004.591    39733.00    458.00    


[rcif-format]: /cryspy/rcif-format

[tof-parameters]: /cryspy/rcif-format/tof-parameters
[tof-profile]: /cryspy/rcif-format/tof-profile
[phase]: /cryspy/rcif-format/phase
[tof-background]: /cryspy/rcif-format/tof-background
[tof-meas]: /cryspy/rcif-format/tof-meas
[diffrn-radiation]: /cryspy/rcif-format/diffrn-radiation
[chi2]: /cryspy/rcif-format/chi2
[range]: /cryspy/rcif-format/range
[extinction]: /cryspy/rcif-format/extinction
[texture]: /cryspy/rcif-format/texture
[exclude]: /cryspy/rcif-format/exclude
[tof-proc]: /cryspy/rcif-format/tof-proc
[tof-peak]: /cryspy/rcif-format/tof-peak
[refine-ls]: /cryspy/rcif-format/refine-ls
[refln]: /cryspy/rcif-format/refln
[refln-susceptibility]: /cryspy/rcif-format/refln-susceptibility
