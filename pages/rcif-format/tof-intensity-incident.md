---
layout: page
permalink: /rcif-format/tof-intensity-incident/
---
[Back to RCIF format][rcif-format]

# `tof_intensity_incident` item/loop:

     _tof_intensity_incident_spectrum Maxwell
     _tof_intensity_incident_a1 0
     _tof_intensity_incident_a2 0
     _tof_intensity_incident_a3 0
     _tof_intensity_incident_a4 0
     _tof_intensity_incident_a5 0
     _tof_intensity_incident_a6 0
     _tof_intensity_incident_a7 0
     _tof_intensity_incident_a8 0


**Mandatory parameters:** 
`a0`, `a1`, `a2`, `a3`, `a4`, `a5`, `a6`, `a7`,
`a8`.

**Optional parameters:** 
`spectrum`.

**Constraints on `spectrum`**: "Maxwell" or "Empirical-Exponents".

For details see the [expression page](https://ikibalin.github.io/cryspy/theory/tof-incident-intensity/).

# Item object

**Methods:** 
`.calc_spectrum(...)`.

# Loop object
**Methods:** 
No

[rcif-format]: /cryspy/rcif-format
