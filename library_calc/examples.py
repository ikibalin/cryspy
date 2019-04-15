# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:08 2019

@author: Iurii Kibalin
"""

import crystal as *

cell = Cell(a=8.25690, singony="Cubic")

#multiplicity 8
fe_a = AtomType(type_n="Fe", type_m="Fe3", x=0.125, y=0.125, z=0.125, occupation=8./192)
#multiplicity 16
fe_b = AtomType(type_n="Fe", type_m="Fe3", x=0.500, y=0.500, z=0.500, occupation=16./192)
#multiplicity 32
o = AtomType(type_n="O", x=0.25223, y=0.25223, z=0.25223, occupation=32./192)

atom_site = AtomSite()
atom_site.add_atom(fe_a)
atom_site.add_atom(fe_b)
atom_site.add_atom(o)

space_groupe = SpaceGroupe(spgr_given_name = "Fd-3m", spgr_choice="2")


crystal = Crystal(space_groupe, cell, atom_site)


h = numpy.array([1, 2, 2, 1, 0, 0], dtype = int) 
k = numpy.array([1, 0, 2, 0, 1, 0], dtype = int)
l = numpy.array([1, 0, 0, 0, 0, 1], dtype = int)
fn = crystal.calc_fn(h, k, l)




for val in fn:
    print("  {:10.4f}".format(val*192))

