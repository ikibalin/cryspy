# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:08 2019

@author: Iurii Kibalin
"""

import matplotlib.pyplot
import numpy

from crystal import *
from calculated_data import *
from experiment_powder_1d import *
from observed_data import *
from Variable import *


cell = Cell(a=8.5502, singony="Cubic")

fe_a = AtomType(type_n="Fe", type_m="Fe3", x=0.125, y=0.125, z=0.125,
                chi_11=-3.616, chi_22=-3.616, chi_33=-3.616)
fe_b = AtomType(type_n="Fe", type_m="Fe3", x=0.500, y=0.500, z=0.500,
                chi_11=3.1846, chi_22=3.1846, chi_33=3.1846)
o = AtomType(type_n="O", x=0.25223, y=0.25223, z=0.25223)


atom_site = AtomSite()
atom_site.add_atom(fe_a)
atom_site.add_atom(fe_b)
atom_site.add_atom(o)

space_groupe = SpaceGroupe(spgr_given_name = "Fd-3m", spgr_choice="2")


crystal = Crystal(space_groupe, cell, atom_site)



# 2 dimensional polarized powder diffraction
calculated_data_powder_2d = CalculatedDataPowder2D(scale=0.0232, crystal=crystal)
experiment_powder_2d= ExperimentPowder2D()
experiment_powder_2d.add_calculated_data(calculated_data_powder_2d )
np_tth = numpy.linspace(2,60,200)
np_phi = numpy.linspace(0,90,100)

int_u, int_d = experiment_powder_2d.calc_profile(np_tth, np_phi)


observed_data_powder_2d = ObservedDataPowder2D()
f_inp = os.path.join("Fe3O4_150K_6T_2d", "full.dat")
observed_data_powder_2d.read_data(f_inp)
experiment_powder_2d.set_val(observed_data=observed_data_powder_2d)




chi_sq_u, n_u, chi_sq_d, n_d, chi_sq_sum, n_sum, chi_sq_dif, n_dif = experiment_powder_2d.calc_chi_sq()




# 1 dimensional polarized powder diffraction
np_tth = numpy.linspace(2,60,200)

int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u+int_d, "k-", np_tth, int_u-int_d, "b-")


calculated_data_powder_1d = CalculatedDataPowder1D(scale=0.0232, crystal=crystal)
#ExperimentPowder1D
experiment_powder_1d= ExperimentPowder1D()
experiment_powder_1d.add_calculated_data(calculated_data_powder_1d)


setup = ExperimentPowder1D().get_val("setup")
setup.set_val(wavelength=0.84, zero_shift=-0.412)
resolution = setup.get_val("resolution")
resolution.set_val(u=17.487, v=-2.8357, w=0.576309)




observed_data = ObservedDataPowder1D()
f_inp = os.path.join("Fe3O4_0T", "full.dat")
observed_data.read_data(f_inp)
experiment_powder_1d.set_val(observed_data=observed_data)


chi_sq_u, n_u, chi_sq_d, n_d, chi_sq_sum, n_sum, chi_sq_dif, n_dif = experiment_powder_1d.calc_chi_sq()

np_tth = observed_data.get_val('tth')
int_u_exp = observed_data.get_val("int_u")
int_d_exp = observed_data.get_val("int_d")
matplotlib.pyplot.plot(np_tth, int_u_exp+int_d_exp, "k-", np_tth, int_u+int_d, "b-")










dd=Variable(val=0)


setup.set_val(zero_shift=dd)




dd["val"]=-0.777
int_u2, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u+int_d, "k-", np_tth, int_u-int_d, "b-")
