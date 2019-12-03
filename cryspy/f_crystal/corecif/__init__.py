"""
Core dictionary (coreCIF)
==============================

This is a not full realization of  the description of 
any crystal structure determination given in 
International Tables for Crystallography (2006) Vol.G, 
Chapter 4.1, pp210-257 "Core dictionary (coreCIF)"
by S.R. Hall, F.H. Allen and I.D. Brown

It includes classes:
------------------------

- Cell() - CELL category
- AtomSite() - ATOM_SITE category
- AtomSiteL() - looped ATOM_SITE category
- AtomSiteAniso() - subcategory ATOM_SITE_ANISO in ATOM_SITE category
- AtomSiteAnisoL() - looped subcategory ATOM_SITE_ANISO in ATOM_SITE category
- AtomSiteAnisoMagnetism() - subcategory ATOM_SITE_ANISO_MAGNETISM in ATOM_SITE category (not defined in IT A)
- AtomSiteAnisoMagnetismL() - looped subcategory ATOM_SITE_ANISO_MAGNETISM in ATOM_SITE category (not defined in IT A)

and additional files:
--------------------------

- ...

Use help function to get more detailed description about used classes:

>>> import cryspy 
>>> help(cryspy.Cell)

"""
__version__ = '1.0.0'
__author__ = 'Iurii Kibalin'
__contact__ = 'iurii.kibalin@cea.fr'
