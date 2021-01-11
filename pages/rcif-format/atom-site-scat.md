---
layout: page
permalink: /rcif-format/atom-site-scat/
---
[Back to RCIF format][rcif-format]

# `atom_site_scat` item/loop:

     loop_
     _atom_site_scat_label
     _atom_site_scat_type_symbol
     _atom_site_scat_Lande 
     _atom_site_scat_kappa 
     Ho Ho3+ 1.250 1.
     Dy Dy3+ 1.333 1.

Description of magnetic structure factor.

**Mandatory parameters:** 
`label`.

**Optional parameters:** 
`type_symbol`, `lande`, `kappa`.

# Item object

**Methods:** 
`.calc_form_factor(...)`,
`.load_atom_type_scat_by_symbol(...)`,
`.plot_form_factor(...)`.

**Protected parameters:** 
`atom_type_scat`.

# Loop object
**Methods:** 
`.calc_form_factor(...)`,
`.load_atom_type_scat_by_atom_site(...)`,
`.plot_form_factor(...)`.

[rcif-format]: /cryspy/rcif-format
