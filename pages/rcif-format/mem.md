---
layout: page
permalink: /rcif-format/mem/
---
[Back to RCIF format][rcif-format]

# `MEM` object:

It describes data blocks of crystals and experiments for magnetization
density reconstruction by Maximum Entropy Method.

**Mandatory data blocks:**
- [crystal][crystal]
- [Diffrn][diffrn]

**Optional Items/Loops:**
- [Density Point][density-point]
- [MEM Parameters][mem-parameters]
- [Section][section]

**Methods:**
`.apply_constraint(...)`, 
`.calc_fr(...)`, 
`.calc_moments_in_unit_cell(...)`, 
`.calc_sections(...)`, 
`.create_prior_density(...)`, 
`.crystals(...)`, 
`.experiments(...)`, 
`.make_cycle(...)`, 
`.maximize_entropy(...)`, 
`.refine_susceptibility(...)`, 
`.save_to_file_den(...)`.

[rcif-format]: /cryspy/rcif-format

[crystal]: /cryspy/rcif-format/crystal
[diffrn]: /cryspy/rcif-format/diffrn
[density-point]: /cryspy/rcif-format/density-point
[mem-parameters]: /cryspy/rcif-format/mem-parameters
[section]: /cryspy/rcif-format/section
