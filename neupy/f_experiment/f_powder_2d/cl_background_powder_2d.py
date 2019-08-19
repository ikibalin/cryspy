"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.interpolate

#class BeamPolarization
from neupy.f_experiment.cl_beam_polarization import BeamPolarization
from neupy.f_common.cl_variable import Variable
#Description of setup class


        
class BackgroundPowder2D(dict):
    """
    Class to describe characteristics of powder diffractometer
    """
    def __init__(self, tth_bkgd=None, phi_bkgd=None, int_bkgd=None, file_dir=None, 
                 file_name=None):
        super(BackgroundPowder2D, self).__init__()
        self._p_file_dir = None
        self._p_file_name = None

        self._p_tth_bkgd = None
        self._p_phi_bkgd = None
        self._p_int_bkgd = None
        
        self._refresh(tth_bkgd, phi_bkgd, int_bkgd, file_dir, file_name)

    def __repr__(self):
        lsout = """BackgroundPowder2D:\n file_dir: {:}
 file_name: {:}""".format(self._p_file_dir, self._p_file_name)
        if self._p_tth_bkgd is not None:
            lsout += "\n ttheta  IntBKGR"
            lsout += ("\n "+7*" "+
                      " ".join(["{:7.2f}".format(hh_1) for hh_1 in 
                                                           self._p_tth_bkgd]))
            for hh_1, l_hh in zip(self._p_phi_bkgd, self._p_int_bkgd):
                lsout += ("\n {:7.2f} ".format(hh_1) +
                          " ".join(["{:}".format(hh_2) for hh_2 in l_hh]))
        return lsout

    def _refresh(self, tth_bkgd, phi_bkgd, int_bkgd, file_dir, file_name):
        if tth_bkgd is not None:
            self._p_tth_bkgd = tth_bkgd
        if phi_bkgd is not None:
            self._p_phi_bkgd = phi_bkgd
        if int_bkgd is not None:
            self._p_int_bkgd = int_bkgd
        if file_dir is not None:
            self._p_file_dir = file_dir
        if file_name is not None:
            self._p_file_name = file_name
            
    def set_val(self, tth_bkgd=None, phi_bkgd=None, int_bkgd=None, 
                file_dir=None, file_name=None):

        self._refresh(tth_bkgd, phi_bkgd, int_bkgd, file_dir, file_name)
        
    def get_val(self, label):
        lab = "_p_"+label
        
        if lab in self.__dict__.keys():
            val = self.__dict__[lab]
            if isinstance(val, type(None)):
                self.set_val()
                val = self.__dict__[lab]
        else:
            print("The value '{:}' is not found".format(lab))
            val = None
        return val

    def list_vals(self):
        """
        give a list of parameters with small descripition
        """
        lsout = """
Parameters:
tth_bkgd is ttheta in degrees to describe background
int_bkgd is intensity to describe background
file_dir
file_name
        """
        print(lsout)
        
        
    def read_data(self):
        """
        read file from file
        """
        f_name = os.path.join(self._p_file_dir, self._p_file_name)
        
        fid=open(f_name,'r')
        l_cont = fid.readlines()
        fid.close()
        
        l_tth = [float(hh) for hh in l_cont[1].strip().split()[1:]]
        l_phi = []
        ll_int = []
        for line in l_cont[2:]:
            l_int = []
            if line.strip() != "":
                l_help = line.strip().split()
                phi = float(l_help[0])
                l_phi.append(phi)
                for hh, tth in zip(l_help[1:], l_tth):
                    l_help_2 = hh.split("(")
                    if len(l_help_2) > 1:
                        val = Variable(val=float(l_help_2[0]), name="IntBKGR_{:.1f}_{:.1f}".format(tth, phi))
                    else:
                        val = float(l_help_2[0])
                    l_int.append(val)
            ll_int.append(l_int)
        
        #n_tth, n_phi = len(l_tth), len(l_phi)
        #ll_int_b = [[ll_int[i_phi][i_tth] for i_phi in range(n_phi)] 
        #                                  for i_tth in range(n_tth)]
        self.set_val(tth_bkgd=l_tth, phi_bkgd=l_phi, int_bkgd=ll_int)
        

    def interpolate_by_points(self, tth, phi):
        l_tth_b = self._p_tth_bkgd
        l_phi_b = self._p_phi_bkgd
        ll_int_b = self._p_int_bkgd
        
        if l_tth_b is None:
            file_dir = self._p_file_dir 
            file_name = self._p_file_name 
            if file_name is None:
                f_inp = os.path.join(file_dir, file_name)
                self.read_data(f_inp)
                l_tth_b = self._p_tth_bkgd
                l_phi_b = self._p_phi_bkgd
                ll_int_b = self._p_int_bkgd
        if l_tth_b is not None:
            tth_b = numpy.array(l_tth_b, dtype=float)
            phi_b = numpy.array(l_phi_b, dtype=float)
            int_b = numpy.array([[1.*hh2 for hh2 in hh1] for hh1 in ll_int_b], dtype=float)
            
            func = scipy.interpolate.interp2d(tth_b, phi_b, int_b)
            
            tth_2d, phi_2d = numpy.meshgrid(tth, phi, indexing="ij")
            
            int_1d = func(tth, phi)
            int_2d = int_1d.transpose()
        else:
            int_2d = numpy.zeros((tth.size, phi.size), dtype=float)
        return int_2d
    
    def get_variables(self):
        l_int_b = self._p_int_bkgd
        l_variable = []
        for l_hh_1 in l_int_b:
            for hh_1 in l_hh_1:
                if isinstance(hh_1, Variable):
                    l_variable.append(hh_1)
        return l_variable    
    
    def create_input_file(self, f_inp=None):
        if f_inp is not None:
            f_dir= os.path.dirname(f_inp)
            f_name= os.path.bathename(f_inp)
            self.set_val(file_dir=f_dir, file_name=f_name)
        f_dir = self._p_file_dir
        f_name = self._p_file_name
        f_full = os.path.join(f_dir, f_name)
        
        s_out = """#ttheta phi IntBKGR  
   3   4.50   40.00   80.00 
  -3   -350    -350    -400
  41   -350    -350    -400"""
    
        fid = open(f_full, "w")
        fid.write(s_out)
        fid.close()

    def save_data(self):
        l_tth_b = self.get_val("tth_bkgd")
        l_phi_b = self.get_val("phi_bkgd")
        ll_int_b = self.get_val("int_bkgd")
        
        if ((l_tth_b is None)|(l_phi_b is None)|(ll_int_b is None)):
            return

        ls_out = ["#   phi\ttheta     IntBKGR"]
        ls_out.append(" {:12}".format(len(l_phi_b))+" ".join(["{:12.3f}".format(hh) for hh in l_tth_b]))
        for phi_b, l_int_b in zip(l_phi_b, ll_int_b):
            ls_out.append(" {:12.3f} ".format(phi_b)+"   ".join(["{:}".format(hh) for hh in l_int_b]))

        f_name = os.path.join(self._p_file_dir, self._p_file_name)
        fid = open(f_name, "w")
        fid.write("\n".join(ls_out))
        fid.close()
