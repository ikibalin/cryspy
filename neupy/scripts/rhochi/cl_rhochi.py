"""
RhoChi program
"""
__author__ = 'ikibalin'
__version__ = "2019_09_10"
import os
import sys
import numpy
import scipy.optimize
import time

from pycifstar import Global

from neupy.f_common.error_simplex import error_estimation_simplex
from neupy.f_crystal.cl_crystal import Crystal
from neupy.f_experiment.f_single.cl_diffrn import Diffrn
from neupy.f_experiment.f_powder_1d.cl_pd import Pd
from neupy.f_experiment.f_powder_2d.cl_pd2d import Pd2d
from neupy.f_experiment.f_powder_texture_2d.cl_pd2dt import Pd2dt


class RhoChi(dict):
    """
    Class to describe rhochi container
    """
    def __init__(self, label="", file_out=None, file_dir=None, experiments=[], variables=[], crystals=[], file_input=None):
        super(RhoChi, self).__init__()
        self.__label = None
        self.__experiments = []
        self.__crystals = []
        self.__file_input = None
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
    def file_input(self):
        return self.__file_input
    @file_input.setter
    def file_input(self, x: str):
        self.__file_input = x

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
            chi_sq, n_points = self.calc_chi_sq()
            return (chi_sq*1./float(n_points))

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
        res = scipy.optimize.minimize(tempfunc, l_param_0, method='BFGS', callback=self._f_callback, options = {"disp": True})

        #res = scipy.optimize.minimize(tempfunc, l_param_0, method='Nelder-Mead', 
        #                              callback=self._f_callback, options = {"fatol": 0.01*n})

        bb = time.time()
        print("refinement complete, time {:.2f} sec.\n\nfinal chi_sq/n: {:.2f}\nstarted chi_sq/n: {:.2f}".format(bb-aa, res.fun, chi_sq*1./n))
        hess_inv = res["hess_inv"]
        sigma = (abs(numpy.diag(hess_inv)*1./float(n)))**0.5
        for fitable, _1  in zip(l_fitable, sigma):
            fitable.sigma = _1

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
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        self.label = cif_global.name
        l_data = cif_global.datas
        for data in l_data:
            flag_crystal, flag_diffrn, flag_pd, flag_pd2d = False, False, False, False
            flag_crystal = data.is_prefix("_cell_length_a")
            flag_diffrn = data.is_prefix("_diffrn_refln") | data.is_prefix("_diffrn_orient_matrix_UB")
            flag_pd = data.is_prefix("_pd_meas") | data.is_prefix("_pd_instr_resolution")
            flag_pd2d = data.is_prefix("_pd2d_meas") | data.is_prefix("_pd2d_instr_resolution")
            flag_pd2dt = data.is_prefix("_2dpdt_texture")
            if flag_crystal:
                crystal = Crystal()
                crystal.from_cif(str(data))
                l_crystal.append(crystal)

            if flag_diffrn:
                diffrn = Diffrn()
                diffrn.from_cif(str(data))
                l_experiment.append(diffrn)

            if flag_pd:
                pd = Pd()
                pd.from_cif(str(data))
                l_experiment.append(pd)

            if (flag_pd2d & (not(flag_pd2dt))):
                pd2d = Pd2d()
                pd2d.from_cif(str(data))
                l_experiment.append(pd2d)
            elif (flag_pd2d & flag_pd2dt):
                pd2dt = Pd2dt()
                pd2dt.from_cif(str(data))
                l_experiment.append(pd2dt)

        self.experiments = l_experiment
        self.crystals = l_crystal

    def read_file(self, f_name):
        self.file_input = f_name
        with open(f_name, "r") as fid:
            string = fid.read()
            self.from_cif(string)

    def save_to_file(self, f_name):
        self.file_input = f_name
        if os.path.basename(f_name) == "main.rcif":
            self.save_to_files()
        else:
            with open(f_name, "w") as fid:
                fid.write(self.to_cif)


    def save_to_files(self):
        if self.file_input is None:
            f_dir = "."
        else:
            f_dir = os.path.dirname(self.file_input)
        f_main = os.path.join(f_dir, "main.rcif")
        ls_main = []
        for crystal in self.crystals:
            ls_main.append("\n"+crystal.to_cif)
        for experiment in self.experiments:   
            ls_main.append("\ndata_{:}".format(experiment.label))
            ls_main.append(experiment.params_to_cif)

            f_data = os.path.join(f_dir, "{:}_data.rcif".format(experiment.label))
            ls_data = []
            ls_data.append("\ndata_{:}".format(experiment.label))
            ls_data.append(experiment.data_to_cif)
            with open(f_data, 'w') as fid:
                fid.write("\n".join(ls_data))
            f_calc = os.path.join(f_dir, "{:}_calc.rcif".format(experiment.label))
            ls_calc = []
            ls_calc.append("\ndata_{:}".format(experiment.label))
            ls_calc.append(experiment.calc_to_cif)
            with open(f_calc, 'w') as fid:
                fid.write("\n".join(ls_calc))
        with open(f_main, 'w') as fid:
            fid.write("\n".join(ls_main))

    def read_files(self, f_dir="."):
        f_main = os.path.join(f_dir, "main.rcif")
        self.file_input = f_main
        
        if not(os.path.isfile(f_main)):
            self._show_message("File '{:}' is not found.".format(f_main))
            return None
        with open(f_main, "r") as fid:
            string = fid.read()
        self.from_cif(string)
        
        for experiment in self.experiments:
            f_data = os.path.join(f_dir, "{:}_data.rcif".format(experiment.label))
            if not(os.path.isfile(f_main)):
                self._show_message("File '{:}' is not found.".format(f_data))
            else:
                with open(f_data, "r") as fid:
                    string = fid.read()
                experiment.from_cif(string)


def rhochi_read_file(f_name):
    rho_chi = RhoChi(file_input=f_name)
    if os.path.basename(f_name) == "main.rcif":
        rho_chi.read_files(os.path.dirname(f_name))
    else:
        rho_chi.read_file(f_name)
    return rho_chi

def rhochi_read_files():
    rho_chi = RhoChi()
    rho_chi.read_files(".")
    return rho_chi

def rhochi_refinement(f_name_in="", f_name_out=""):
    """
    refinement,
    parameters are defined in given .rcif fiel
    """
    print(70*"*"+"\n"+"RhoChi program. Console version.".center(70)+"\n"+
          70*"*")
    if f_name_in != "":
        rho_chi = rhochi_read_file(f_name_in)
    else:
        rho_chi = rhochi_read_files()
    print("Before refinement:\n")
    #print(rho_chi)
    print("\nRefined parameters -- before:\n")
    for fitable in rho_chi.get_variables():
        print(fitable)

    rho_chi.refine()

    print("After refinement:\n")
    print(rho_chi)
    print("\nRefined parameters -- after:\n")
    for fitable in rho_chi.get_variables():
        print(fitable)

    if f_name_out != "":
        rho_chi.save_to_file(f_name_out)
    else:
        rho_chi.save_to_files()

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

