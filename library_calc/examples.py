# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:08 2019

@author: Iurii Kibalin
"""

import matplotlib
import matplotlib.pyplot


from crystal import *
from calculated_data import *
from experiment_pd import *
from observed_data import *


cell = Cell(a=8.25690, singony="Cubic")

fe_a = AtomType(type_n="Fe", type_m="Fe3", x=0.125, y=0.125, z=0.125,
                chi_11=2.3, chi_22=2.3, chi_33=2.3)
fe_b = AtomType(type_n="Fe", type_m="Fe3", x=0.500, y=0.500, z=0.500,
                chi_11=-2.3, chi_22=-2.3, chi_33=-2.3)
o = AtomType(type_n="O", x=0.25223, y=0.25223, z=0.25223)

atom_site = AtomSite()
atom_site.add_atom(fe_a)
atom_site.add_atom(fe_b)
atom_site.add_atom(o)

space_groupe = SpaceGroupe(spgr_given_name = "Fd-3m", spgr_choice="2")


crystal = Crystal(space_groupe, cell, atom_site)
crystal.calc_fn(2,3,0)


calculated_data_1d_pd = CalculatedData1DPD(scale=1.0, crystal=crystal)

experiment_1d_pd = Experiment1DPD()
experiment_1d_pd.add_calculated_data(calculated_data_1d_pd)


setup = Experiment1DPD().get_val("setup")
setup.set_val(wavelength=0.84)
resolution = setup.get_val("resolution")
resolution.set_val(u=0.2, w=0.5)

observed_data = ObservedData1DPD()
f_inp = os.path.join("Fe3O4_0T", "full.dat")
observed_data.read_data(f_inp)

experiment_1d_pd.set_val(observed_data=observed_data)



chi_sq_u, chi_sq_d, chi_sq_sum, chi_sq_diff = experiment_1d_pd.calc_chi_sq()




np_tth=numpy.linspace(1,100,300)
int_u, int_d = experiment_1d_pd.calc_profile(np_tth)


matplotlib.pyplot.plot(observed_data.get_val("tth"), 
   observed_data.get_val("int_u")-observed_data.get_val("int_d"))

matplotlib.pyplot.plot(observed_data.get_val("tth"), 
                       observed_data.get_val("int_u"))



