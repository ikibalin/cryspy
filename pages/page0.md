---
layout: page
title: RCIF format
permalink: /rcif-format/
---

The CrysPy library works with objects defined in the RCIF format.
Inspired by [CIF dictionaries][cif-format] RCIF extends it to describe 
various properties needed for data refinement.

# Global containers
- RhoChi
- MEM

# Data containers
- Crystal
- Diffrn
- Pd
- Pd2d
- TOF

# Items and loops
- [Space Group][space-group]
- [Cell][cell]
- [Atom Site][atom-site]
- [Atom Type][atom-type]
- [Atom Susceptibility][atom-susceptibility]
- [Atom Site Aniso][atom-site-aniso]
- [Inversed Hessian][inversed-hessian]
- [Refln][refln]
- [Atom Site Susceptibility][atom-site-susceptibility]
- [Atom Type Scat][atom-type-scat]
- [Refln Susceptibility][refln-susceptibility]
- [Atom Local Axes][atom-local-axes]
- [Atom Electron Configuration][atom-electron-configuration]
- [Chi2][chi2]
- [Exclude][exclude]
- [Diffrn Radiation][diffrn-radiation]
- [Diffrn Refln][diffrn-refln]
- Extinction
- MEM Parameters
- Phase
- Pd Background
- Pd Instr Reflex Asymmetry
- Pd Instr Resolution
- Pd Meas
- Pd Proc
- Pd Peak
- Pd2d Background
- Pd2d Instr Reflex Asymmetry
- Pd2d Instr Resolution
- Pd2dMeas
- Pd2dProc
- Pd2dPeak
- Range
- Refln
- Refln Susceptibility
- Refine Ls
- Setup
- Texture
- TOF Background
- TOF Meas
- TOF Parameters
- TOF Peak
- TOF Proc
- TOF Profile
- AtomSiteScat
- Diffrn Orient Matrix
- Section
- SpaceGroup
- Space Group Symop Magn Operation


[cif-format]: https://www.iucr.org/resources/cif/dictionaries

[space-group]: /cryspy/rcif-format/space-group
[cell]: /cryspy/rcif-format/cell
[atom-site]: /cryspy/rcif-format/atom-site
[atom-type]: /cryspy/rcif-format/atom-type
[atom-susceptibility]: /cryspy/rcif-format/atom-susceptibility
[atom-site-aniso]: /cryspy/rcif-format/atom-site-aniso
[inversed-hessian]: /cryspy/rcif-format/inversed-hessian
[refln]: /cryspy/rcif-format/refln
[atom-site-susceptibility]: /cryspy/rcif-format/atom-site-susceptibility
[atom-type-scat]: /cryspy/rcif-format/atom-type-scat
[refln-susceptibility]: /cryspy/rcif-format/refln-susceptibility
[atom-local-axes]: /cryspy/rcif-format/atom-local-axes
[atom-electron-configuration]: /cryspy/rcif-format/atom-electron-configuration
[chi2]: /cryspy/rcif-format/chi2
[exclude]: /cryspy/rcif-format/exclude
[diffrn-radiation]: /cryspy/rcif-format/diffrn-radiation
[diffrn-refln]: /cryspy/rcif-format/diffrn-refln