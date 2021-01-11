---
layout: page
permalink: /rcif-format/tof-background/
---
[Back to RCIF format][rcif-format]

# `tof_background` item/loop:

     _tof_background_time_max  30000.00
     _tof_background_coeff1    24832.850
     _tof_background_coeff2     6139.244
     _tof_background_coeff3     8063.472
     _tof_background_coeff4     3125.050
     _tof_background_coeff5     2566.956
     _tof_background_coeff6      311.077
     _tof_background_coeff7      837.348
     _tof_background_coeff8     -103.742
     _tof_background_coeff9      -11.806


**Mandatory parameters:** 
`time_max` in microseconds, `coeff1`.

**Optional parameters:** 
`coeff2`, `coeff3`, `coeff4`, `coeff5`, `coeff6`, `coeff7`, `coeff8`,
`coeff9`, `coeff10`, `coeff11`, `coeff12`, `coeff13`, `coeff14`,
`coeff15`, `coeff16`, `coeff17`, `coeff18`, `id`.

For details see [the expression page](https://ikibalin.github.io/cryspy/theory/tof-background/).


# Item object

**Methods:** 
`.calc_background(...)`.

# Loop object
**Methods:** 
No

[rcif-format]: /cryspy/rcif-format
