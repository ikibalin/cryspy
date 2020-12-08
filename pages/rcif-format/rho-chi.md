---
layout: page
permalink: /rcif-format/rho-chi/
---
[Back to RCIF format][rcif-format]

# `RhoChi` global object:

It describes data blocks of crystals and experiments for data refinement

**Crystals:**
- [Crystal][crystal]
- [Mag Crystal][mag-crystal]

**Experiments:**
- [Diffrn][diffrn]
- [Pd][pd]
- [Pd2d][pd2d]
- [TOF][tof]

**Optional Items/Loops:**
- [Inversed Hessian][inversed-hessian]

**Methods:**
`.apply_constraints(...)`, 
`.calc_chi_sq(...)`, 
`.crystals(...)`, 
`.estimate_f_mag_for_diffrn()`, 
`.estimate_inversed_hessian(...)`, 
`.experiments(...)`, 
`.refine(...)`, 
`.save_to_file(...)`, 
`.save_to_files(...)`.

[rcif-format]: /cryspy/rcif-format

[crystal]: /cryspy/rcif-format/crystal
[mag-crystal]: /cryspy/rcif-format/mag-crystal
[diffrn]: /cryspy/rcif-format/diffrn
[pd]: /cryspy/rcif-format/pd
[pd2d]: /cryspy/rcif-format/pd2d
[tof]: /cryspy/rcif-format/tof
[inversed-hessian]: /cryspy/rcif-format/inversed-hessian
