---
layout: page
permalink: /rcif-format/tof-profile/
---
[Back to RCIF format][rcif-format]

# `tof_profile` item/loop:

     _tof_profile_alpha0     0.000
     _tof_profile_alpha1     0.29710
     _tof_profile_beta0      0.04182
     _tof_profile_beta1      0.00224
     _tof_profile_sigma0     0.409(35)
     _tof_profile_sigma1     8.118(30)
     _tof_profile_peak_shape  Gauss


**Mandatory parameters:** 
`alpha0`, `alpha1`, `beta0`, `beta1`, 
`sigma0`,
`sigma1`.

**Optional parameters:** 
`sigma2`, `peak_shape`, `size_g`, `size_l`, `strain_g`, `strain_l`,
`gamma2`, `gamma1`, `gamma0`.

**Constraints on `peak_shape`**: "pseudo-Voigt", "Gauss".

For details see [the expression page](https://ikibalin.github.io/cryspy/theory/tof-profile/).

# Item object

**Methods:** 
`.calc_hpveta(...)`,
`.calc_agbg(...)`,
`.calc_albl(...)`,
`.calc_peak_shape_function(...)`.

# Loop object
**Methods:** 
No

[rcif-format]: /cryspy/rcif-format
