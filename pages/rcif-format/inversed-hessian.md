---
layout: page
permalink: /rcif-format/inversed-hessian/
---
[Back to RCIF format][rcif-format]

# `inversed_hessian` item/loop:

     _inversed_hessian_with_labels
     ;
     length_a 8 3 0
     length_b 3 9 2
     length_c 2 8 3
     ;

**Mandatory parameters:** 
`with_labels`

# Item object

**Methods:** 
`.set_labels(...)`, 
`.set_inversed_hessian(...)`, 
`.set_correlation_matrix(...)`, 
`.set_sigmas(...)`, 
`.report(...)`.

**Internal parameters:** 
`label`, `matrix`, `correlation_matrix`, `sigma`

[rcif-format]: /cryspy/rcif-format