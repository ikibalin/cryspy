"""
RhoChi program
"""
__author__ = 'ikibalin'
__version__ = "2019_09_09"
import os
import sys
import numpy
import scipy.optimize
import time

from pystar import CIFglobal

from neupy.f_common.error_simplex import error_estimation_simplex
from neupy.f_crystal.cl_crystal import Crystal
from neupy.f_experiment.f_single.cl_diffrn import Diffrn
from neupy.f_experiment.f_powder_1d.cl_pd import Pd


class RhoChi(dict):
    """
    Class to describe rhochi container
    """
    def __init__(self, label="", file_out=None, file_dir=None, experiments=[], variables=[], crystals=[]):
        super(RhoChi, self).__init__()
        self.__label = None
        self.__experiments = []
        self.__crystals = []
        self.__file_out = None
        self.__file_dir = None

    @property
    def label(self):
        return self.__label
    @label.setter
    def label(self, x: str):
        self.__label = x

    @property
    def experiments(self):
        return self.__experiments
    @experiments.setter
    def experiments(self, l_x):
        self.__experiments = l_x

    @property
    def crystals(self):
        return self.__crystals
    @crystals.setter
    def crystals(self, l_x):
        self.__crystals = l_x

    @property
    def file_out(self):
        return self.__file_out
    @file_out.setter
    def file_out(self, x: str):
        self.__file_out = x

    @property
    def file_dir(self):
        return self.__file_dir
    @file_dir.setter
    def file_dir(self, x: str):
        self.__file_dir = x

    def __repr__(self):
        ls_out = []
        ls_out.append("RhoChi:")
        ls_out.append(" label: ".format(self.label))
        ls_out.append(" file_out: ".format(self.file_out))
        ls_out.append("\n Crystals:")
        ls_out.append("\n".join([str(_) for _ in self.crystals]))
        ls_out.append("\n Experiments:")
        ls_out.append("\n".join([str(_) for _ in self.experiments]))
        return "\n".join(ls_out)






    def add_experiment(self, experiment):
        self.__experiments.append(experiment)

    def del_experiment(self, ind):
        self.__experiments.pop(ind)        

    def replace_experiment(self, ind, experiment):
        self.__experiments.pop(ind)
        self.__experiments.insert(ind, experiment)

    def add_crystal(self, crystal):
        self.__crystals.append(crystal)

    def del_crystal(self, ind):
        self.__crystals.pop(ind)        

    def replace_crystal(self, ind, crystal):
        self.__crystals.pop(ind)
        self.__crystals.insert(ind, crystal)

    
    def calc_chi_sq(self):
        """
        calculate the integral intensity for h, k, l reflections
        """
        self.apply_constraint()
            
        l_crystal = self.crystals

        chi_sq_res, n_res = 0., 0.
        for experiment in self.experiments:
            chi_sq, n = experiment.calc_chi_sq(l_crystal)
            chi_sq_res += chi_sq
            n_res += n
        return chi_sq_res, n_res
    
    def apply_constraint(self):
        for crystal in self.crystals:
            crystal.apply_constraint()

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    def refine(self):
        """
        optimization
        """
        self.apply_constraint()
        l_fitable = self.get_variables()
        if l_fitable == []:
            self._show_message("Variables are not found")
            return None
        
        l_param_0 = [fitable.value for fitable in l_fitable]

        def tempfunc(l_param):
            for fitable, param in zip(l_fitable, l_param):
                fitable.value = param
            chi_sq = self.calc_chi_sq()[0]
            return chi_sq

        chi_sq, n = self.calc_chi_sq()

        #if self.ref[0].refin:
        print("starting chi_sq/n: {:.2f} (n = {:}).".format(chi_sq*1./n, int(n)))
        print("\nrefinement started for parameters:")
        ls_out = " ".join(["{:12}".format(fitable.name.rjust(12)) if len(fitable.name)<=12 
                           else "{:12}".format(fitable.name[-12:]) for fitable in l_fitable]) + "       chi_sq"
        print(ls_out)
        aa = time.time()
        """
        res, m_error, infodict, errmsg, ier = \
            scipy.optimize.leastsq(tempfunc, l_param_0, full_output=1)

        """
        res = scipy.optimize.minimize(tempfunc, l_param_0, method='Nelder-Mead', 
                                      callback=self._f_callback, options = {"fatol": 0.01*n})
        bb = time.time()
        print("refinement complete, time {:.2f} sec.\n\nfinal chi_sq/n: {:.2f}\nstarted chi_sq/n: {:.2f}".format(bb-aa, res.fun*1./n, chi_sq*1./n))
        
        ##it is not checked
        #m_error, dist_hh = error_estimation_simplex(res["final_simplex"][0], res["final_simplex"][1], tempfunc)
        #
        #l_sigma = []
        #for i, val_2 in zip(range(m_error.shape[0]), dist_hh):
        #    #slightly change definition, instead of (n-k) here is n
        #    error = (abs(m_error[i,i])*1./n)**0.5
        #    if m_error[i,i] < 0.:
        #        print(50*"*"+"\nError is incorrect\n(negative diagonal elements of Hessian)\n"+50*"*")
        #    if val_2 > error:
        #        print(50*"*"+"\nErrors is incorrect\n(minimum is not found)\n"+50*"*")
        #        
        #    l_sigma.append(max(error,val_2))
        #    
        #for variable, sigma in zip(l_variable, l_sigma):
        #    variable[4] = sigma
        return res

    
    def _f_callback(self, *arg):
        res_x = arg[0]
        ls_out = " ".join(["{:12.5f}".format(hh) for hh in res_x])
        if len(arg) > 1:
            res_fun = arg[1]
            ls_out += " {:12.1f}".format(res_fun.fun)
        print(ls_out)
    
    @property
    def is_variable(self):
        res = (any([_.is_variable for _ in self.experiments]) |
               any([_.is_variable for _ in self.crystals]))
        return res

    def get_variables(self):
        l_variable = []
        for _ in self.crystals:
            l_variable.extend(_.get_variables())
        for _ in self.experiments:
            l_variable.extend(_.get_variables())
        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        ls_out.append("global_{:}\n".format(self.label))
        if self.crystals is not None:
            ls_out.extend([_.to_cif for _ in self.crystals])
        if self.experiments is not None:
            ls_out.extend([_.to_cif for _ in self.experiments])
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        l_crystal, l_experiment = [], []
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        self.label = cif_global.name
        l_cif_data = cif_global.datas
        for cif_data in l_cif_data:
            flag_crystal, flag_diffrn, flag_pd = False, False, False
            flag_crystal = cif_data.is_prefix("_cell_length_a")
            flag_diffrn = cif_data.is_prefix("_diffrn_refln")
            flag_pd = cif_data.is_prefix("_pd_meas")
            if flag_crystal:
                crystal = Crystal()
                crystal.from_cif(str(cif_data))
                l_crystal.append(crystal)

            if flag_diffrn:
                diffrn = Diffrn()
                diffrn.from_cif(str(cif_data))
                l_experiment.append(diffrn)

            if flag_pd:
                pd = Pd()
                pd.from_cif(str(cif_data))
                l_experiment.append(pd)
        self.experiments = l_experiment
        self.crystals = l_crystal


    
def rhochi_refinement(f_name_in, f_name_out):
    """
    refinement,
    parameters are defined in given .rcif fiel
    """
    print(70*"*"+"\n"+"RhoChi program. Console version.".center(70)+"\n"+
          70*"*")
    rho_chi = RhoChi()
    with open(f_name_in, 'r') as fid:
        string = fid.read()

    rho_chi.from_cif(string)
    print("Before refinement:\n")
    print(rho_chi)
    print("\nRefined parameters -- before:\n")
    for fitable in rho_chi.get_variables():
        print(fitable)

    rho_chi.refine()

    print("After refinement:\n")
    print(rho_chi)
    print("\nRefined parameters -- after:\n")
    for fitable in rho_chi.get_variables():
        print(fitable)
    
    string_out = rho_chi.to_cif
    with open(f_name_out, 'w') as fid:
        fid.write(string_out)

    print(70*"*"+"\n"+70*"*")


def create_temporary(f_name_in, exp_type="2"):
    f_dir = os.path.dirname(f_name_in)
    model = Model()

    atom_type_1 = AtomType(flag_m=True)
    crystal_name = "Phase1"
    crystal = Crystal(label=crystal_name)
    crystal.add_atom(atom_type_1)
    model.add_crystal(crystal)
    
    if "1" in exp_type:
        file_out = "full_sd.out"
        observed_data = ObservedDataSingle(file_dir=f_dir, file_name="full_sd.dat")
        observed_data.create_input_file()
        experiment = ExperimentSingle(label="exp_sd", observed_data=observed_data,
                                      file_out=file_out, file_dir=f_dir)

        calculated_data = CalculatedDataSingle(label=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)

    if "2" in exp_type:
        file_out = "full_pd.out"
        observed_data = ObservedDataPowder1D(file_dir=f_dir, 
                        file_name="full_pd.dat", tth_min=2., tth_max=100)
        observed_data.create_input_file()
        
        background = BackgroundPowder1D(file_dir=f_dir, file_name="full_pd.bkg")
        background.create_input_file()
        
        setup = SetupPowder1D(background=background)
        
        experiment = ExperimentPowder1D(label="exp_pd", setup=setup, file_dir=f_dir, 
                        file_out=file_out, observed_data=observed_data, excl_tth_min=[35], excl_tth_max=[40])

        calculated_data = CalculatedDataPowder1D(label=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)

    if "3" in exp_type:
        file_out = "full_2dpd.out"
        observed_data = ObservedDataPowder2D(file_dir=f_dir, file_name="full_2dpd.dat",
        tth_min=2., tth_max=80, phi_min=0., phi_max=20.)
        observed_data.create_input_file()
        
        background = BackgroundPowder2D(file_dir=f_dir, file_name="full_2dpd.bkg")
        background.create_input_file()
        
        setup = SetupPowder2D(background=background)
        
        experiment = ExperimentPowder2D(label="exp_2dpd", setup=setup, file_out=file_out,
                               file_dir=f_dir, observed_data=observed_data)

        calculated_data = CalculatedDataPowder2D(label=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)
    write_to_rcif(model, f_name_in)

