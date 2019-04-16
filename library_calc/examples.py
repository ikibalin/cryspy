# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:08 2019

@author: Iurii Kibalin
"""

from crystal import *
from calculated_data import *
from experiment import *
from setup_1d_pd import *


cell = Cell(a=8.25690, singony="Cubic")

#multiplicity 8
fe_a = AtomType(type_n="Fe", type_m="Fe3", x=0.125, y=0.125, z=0.125)
#multiplicity 16
fe_b = AtomType(type_n="Fe", type_m="Fe3", x=0.500, y=0.500, z=0.500)
#multiplicity 32
o = AtomType(type_n="O", x=0.25223, y=0.25223, z=0.25223)

atom_site = AtomSite()
atom_site.add_atom(fe_a)
atom_site.add_atom(fe_b)
atom_site.add_atom(o)

space_groupe = SpaceGroupe(spgr_given_name = "Fd-3m", spgr_choice="2")



crystal = Crystal(space_groupe, cell, atom_site)

calculated_data_1d_pd = CalculatedData1DPD(1.0, crystal)

h = numpy.array([1, 2, 2, 1, 0, 0], dtype = int) 
k = numpy.array([1, 0, 2, 0, 1, 0], dtype = int)
l = numpy.array([1, 0, 0, 0, 0, 1], dtype = int)
iint_hkl = calculated_data_1d_pd.calc_iint(h,k,l)




resolution = Resolution1DPD()
factor_lorentz_1d_pd = FactorLorentz1DPD()
asymmetry_1d_pd = Asymmetry1DPD()
beam_polarization = BeamPolarization()
Setup1DPD()


for val in fn:
    print("  {:10.4f}".format(val))

