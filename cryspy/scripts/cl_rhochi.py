__author__ = 'ikibalin'
__version__ = "2019_09_15"
import os
import numpy

from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_data_constr import DataConstr
from cryspy.common.cl_global_constr import GlobalConstr
from cryspy.common.cl_fitable import Fitable

import scipy.optimize
import time


from cryspy.cif_like.cl_crystal import Crystal
from cryspy.cif_like.cl_diffrn import Diffrn
from cryspy.cif_like.cl_pd import Pd
from cryspy.cif_like.cl_pd2d import Pd2d



class RhoChi(GlobalConstr):
    """
Class to describe RhoChi container

Description in cif file::

 global_

 data_Fe3O4             
 _cell_angle_alpha 90.0                    
 _cell_angle_beta 90.0
 _cell_angle_gamma 90.0
 _cell_length_a 8.56212()
 _cell_length_b 8.56212
 _cell_length_c 8.56212
 _space_group_it_coordinate_system_code 2  
 _space_group_IT_number    227
 
 loop_                                     
 _atom_site_adp_type
 _atom_site_B_iso_or_equiv
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_label
 _atom_site_occupancy
 _atom_site_type_symbol
  Uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
  Uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
  Uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-
 
 loop_                                     
 _atom_type_scat_length_neutron
 _atom_type_symbol
   0.945 Fe3+
  0.5803 O2-
 
 loop_
 _atom_site_aniso_U_11
 _atom_site_aniso_U_12
 _atom_site_aniso_U_13
 _atom_site_aniso_U_22
 _atom_site_aniso_U_23
 _atom_site_aniso_U_33
 _atom_site_aniso_label
  0.0 0.0 0.0 0.0 0.0 0.0 Fe3A
  0.0 0.0 0.0 0.0 0.0 0.0 Fe3B
 
 loop_
 _atom_site_scat_label
 _atom_site_scat_lande
 Fe3A 2.0 
 Fe3B 2.0 
 
 loop_     
 _atom_site_susceptibility_label
 _atom_site_susceptibility_chi_type
 _atom_site_susceptibility_chi_11
 _atom_site_susceptibility_chi_12
 _atom_site_susceptibility_chi_13
 _atom_site_susceptibility_chi_22
 _atom_site_susceptibility_chi_23
 _atom_site_susceptibility_chi_33
  Fe3A Cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
  Fe3B Cani 3.041      0.0 0.0  3.041 0.0  3.041

 data_mono
 _setup_wavelength     0.840
 _setup_field          1.000
 
 _diffrn_radiation_polarization 1.0
 _diffrn_radiation_efficiency   1.0

 _extinction_mosaicity 100.0
 _extinction_radius    50.0
 _extinction_model     gauss

 _diffrn_orient_matrix_UB_11 6.59783
 _diffrn_orient_matrix_UB_12 -6.99807
 _diffrn_orient_matrix_UB_13 3.3663
 _diffrn_orient_matrix_UB_21 2.18396
 _diffrn_orient_matrix_UB_22 -2.60871
 _diffrn_orient_matrix_UB_23 -9.5302
 _diffrn_orient_matrix_UB_31 7.4657
 _diffrn_orient_matrix_UB_32 6.94702
 _diffrn_orient_matrix_UB_33 -0.18685

 _phase_label  Fe3O4

 loop_
 _diffrn_refln_index_h
 _diffrn_refln_index_k
 _diffrn_refln_index_l
 _diffrn_refln_fr
 _diffrn_refln_fr_sigma
 0 0 8 0.64545 0.01329 
 2 0 6 1.75682 0.04540 
 0 2 6 1.67974 0.03711 
    """
    MANDATORY_CLASSES = (Crystal, )
    OPTIONAL_CLASSES = (Diffrn, Pd, Pd2d)
    INTERNAL_CLASSES = ()
    def __init__(self, crystals=None, experiments=None,
                 global_name=""):
        super(RhoChi, self).__init__(mandatory_classes=self.MANDATORY_CLASSES,
                                     optional_classes=self.OPTIONAL_CLASSES,
                                     internal_classes=self.INTERNAL_CLASSES)
        self.global_name = global_name
        self.crystals = crystals
        self.experiments = experiments

        if self.is_defined:
            self.form_object

    @property
    def crystals(self):
        """
        """
        l_res = self[Crystal]
        if len(l_res) >= 1:
            return l_res
        else:
            return None
    @crystals.setter
    def crystals(self, l_x):
        self.delete_crystals
        l_x_in = []
        if l_x is None:
            pass
        else:
            l_x_in = [x for x in l_x if isinstance(x, Crystal)]
        self.mandatory_objs.extend(l_x_in)

    @property
    def experiments(self):
        """
        """
        l_res = self[Diffrn]+self[Pd]+self[Pd2d]
        if len(l_res) >= 1:
            return l_res
        else:
            return None
    @experiments.setter
    def experiments(self, l_x):
        self.delete_experiments
        l_x_in = []
        if l_x is None:
            pass
        else:
            l_x_in = [x for x in l_x if (isinstance(x, Diffrn) | isinstance(x, Pd)
                                       | isinstance(x, Pd2d))]
        self.optional_objs.extend(l_x_in)


    @property
    def delete_crystals(self):
        _h = [self.mandatory_objs.remove(obj) for obj in reversed(self.mandatory_objs) if isinstance(obj, Crystal)]
        return True

    @property
    def delete_experiments(self):
        _h = [self.optional_objs.remove(obj) for obj in reversed(self.mandatory_objs) if (
                  isinstance(obj, Diffrn) | isinstance(obj, Pd) | isinstance(obj, Pd2d))]
        return True


    def __repr__(self):
        ls_out = ["RhoChi:"]
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)



    def apply_constraint(self):
        res = all([crystal.apply_constraint() for crystal in self.crystals])
        return res


    def calc_chi_sq(self):
        """
calculate the integral intensity for h, k, l reflections
        """
        self.apply_constraint()
            
        l_crystal = self.crystals

        chi_sq_res, n_res = 0., 0.
        for experiment in self.experiments:
            chi_sq, n = experiment.calc_chi_sq(l_crystal)
            experiment.chi_sq = chi_sq
            experiment.n = n
            chi_sq_res += chi_sq
            n_res += n
        return chi_sq_res, n_res

    @property    
    def remove_internal_objs(self):
        for _obj in self.crystals:
            _obj.remove_internal_objs
        for _obj in self.experiments:
            _obj.remove_internal_objs

    def refine(self):
        """
Minimization procedure
        """
        flag = True
        
        self.remove_internal_objs

        self.apply_constraint()
        l_fitable = self.get_variables()

        if l_fitable == []:
            chi_sq, n = self.calc_chi_sq()
            #self._show_message(f"chi_sq/n {chi_sq/n:.2f} (n = {int(n):}).")
            _dict_out = {"flag": flag, "res":None, "chi_sq":chi_sq, "n":n}
            return _dict_out
        

        val_0 = numpy.array([fitable.value for fitable in l_fitable], dtype=float)
        
        sign = 2*(numpy.array(val_0 >= 0., dtype=int)-0.5)
        param_0 = numpy.log(abs(val_0)*(numpy.e-1.)+1.)*sign
        coeff_norm = numpy.where(val_0 == 0., 1., val_0)/numpy.where(param_0==0., 1., param_0)
        hes_coeff_norm = numpy.matmul(coeff_norm[:, numpy.newaxis], coeff_norm[numpy.newaxis, :])

        chi_sq, n = self.calc_chi_sq()

        def tempfunc(l_param):
            for fitable, param, _1 in zip(l_fitable, l_param, coeff_norm):
                fitable.value = param*_1
            chi_sq, n_points = self.calc_chi_sq()
            if n_points < n:
                res_out = 1.0e+308
            else:
                res_out = (chi_sq*1./float(n_points))
            #print("l_param: ", l_param)
            #print("chi_sq/n_points: ", res_out)
            return res_out


        #if self.ref[0].refin:
        #print("starting chi_sq/n: {:.2f} (n = {:}).".format(chi_sq*1./n, int(n)))
        #print("\nrefinement started for parameters:")
        #ls_out = " ".join(["{:12}".format(fitable.name.rjust(12)) if len(fitable.name)<=12 
        #                   else "{:12}".format(fitable.name[-12:]) for fitable in l_fitable]) + "       chi_sq"
        #print(ls_out)
        aa = time.time()
        """
        res, m_error, infodict, errmsg, ier = \
            scipy.optimize.leastsq(tempfunc, param_0, full_output=1)

        """
        res = scipy.optimize.minimize(tempfunc, param_0, method='BFGS', callback=lambda x : self._f_callback(coeff_norm, x), options = {"disp": True})
        
        _dict_out = {"flag": flag, "res":res}
        #res = scipy.optimize.minimize(tempfunc, l_param_0, method='Nelder-Mead', 
        #                              callback=self._f_callback, options = {"fatol": 0.01*n})

        bb = time.time()
        #print("refinement complete, time {:.2f} sec.\n\nfinal chi_sq/n: {:.2f}\nstarted chi_sq/n: {:.2f}".format(bb-aa, res.fun, chi_sq*1./n))
        
        hess_inv = res["hess_inv"] * hes_coeff_norm
        sigma = (abs(numpy.diag(hess_inv)*1./float(n)))**0.5
        for fitable, _1  in zip(l_fitable, sigma):
            fitable.sigma = _1
        #print("experiment  chi_sq_n")
        #for _1 in self.experiments:
        #    print("{:10}: {:8.2f}".format(_1.label, _1.chi_sq/_1.n))
        return _dict_out

    
    def _f_callback(self, *arg):
        coeff_norm = arg[0]
        res_x = arg[1]
        ls_out = " ".join(["{:12.5f}".format(_1*_2) for _1, _2 in zip(res_x, coeff_norm)])
        if len(arg) > 2:
            res_fun = arg[1]
            ls_out += " {:12.1f}".format(res_fun.fun)
        #print(ls_out)
  


    def read_file(self, f_name):
        self.file_input = f_name
        star_ = to_global(f_name)
        string = str(star_)
        self.from_cif(string)
        self.apply_constraint()

    def save_to_file(self, f_name):
        self.file_input = f_name
        if os.path.basename(f_name) == "main.rcif":
            self.save_to_files()
        else:
            with open(f_name, "w") as fid:
                fid.write(self.to_cif())


    def save_to_files(self):
        if self.file_input is None:
            f_dir = "."
        else:
            f_dir = os.path.dirname(self.file_input)
        f_main = os.path.join(f_dir, "main.rcif")
        ls_main = []
        ls_main.append("global_{:}\n".format(self.global_name))
        for experiment in self.experiments:   
            ls_main.append(f"\n_add_url {experiment.data_name:}_data.rcif\n")
            ls_main.append(f"_add_url {experiment.data_name:}_calc.rcif\n")
        for crystal in self.crystals:
            ls_main.append("\n"+crystal.to_cif())
        for experiment in self.experiments:   
            ls_main.append("\ndata_{:}".format(experiment.data_name))
            ls_main.append(experiment.params_to_cif())

            f_data = os.path.join(f_dir, "{:}_data.rcif".format(experiment.data_name))
            ls_data = []
            ls_data.append("\ndata_{:}".format(experiment.data_name))
            ls_data.append(experiment.data_to_cif())
            with open(f_data, 'w') as fid:
                fid.write("\n".join(ls_data))
            f_calc = os.path.join(f_dir, "{:}_calc.rcif".format(experiment.data_name))
            ls_calc = []
            ls_calc.append("\ndata_{:}".format(experiment.data_name))
            ls_calc.append(experiment.calc_to_cif())
            with open(f_calc, 'w') as fid:
                fid.write("\n".join(ls_calc))
        with open(f_main, 'w') as fid:
            fid.write("\n".join(ls_main))
