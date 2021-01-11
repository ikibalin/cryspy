---
layout: page
permalink: /rcif-format/tof-parameters/
---
[Back to RCIF format][rcif-format]

# `tof_parameters` item/loop:

     _tof_parameters_Zero    2.92100
     _tof_parameters_Dtt1 6167.24700
     _tof_parameters_neutrons thermal
     _tof_parameters_dtt2 -2.28000
     _tof_parameters_2theta_bank 145.00
     _tof_parameters_extinction 0.


**Mandatory parameters:** 
`zero`: microseconds, `dtt1`: in microseconds/Angstrem, `ttheta_bank`: degrees.

**Optional parameters:** 
`neutrons`, `dtt2`: in microseconds/Angstrem^2, `zerot`: in microseconds, `dtt1t`: in microseconds/Angstrem,
`dtt2t`: in microseconds/Angstrem^2,
`width`, `x_cross`, `field`: in Tesla, `extinction`.

**Constraint on `neutrons`**: "thermal", "epithermal"

For details see [the expression page](https://ikibalin.github.io/cryspy/theory/tof-parameters/).

# Item object

**Methods:** 
No

**Internal parameters:** 
No 

**Protected parameters:** 
No

# Loop object
**Methods:** 
No

[rcif-format]: /cryspy/rcif-format
