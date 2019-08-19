"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
    
#Description of setup class

class BackgroundPowder1D(dict):
    """
    Class to describe characteristics of powder diffractometer
    """
    def __init__(self, tth_bkgd=None, int_bkgd=None, file_dir=None, file_name=None):
        super(BackgroundPowder1D, self).__init__()
        self._p_file_dir = None
        self._p_file_name = None

        self._p_tth_bkgd = None
        self._p_int_bkgd = None
        
        self._refresh(tth_bkgd, int_bkgd, file_dir, file_name)

    def __repr__(self):
        lsout = """BackgroundPowder1D:\n file_dir: {:}
 file_name: {:}""".format(self._p_file_dir, self._p_file_name)
        if self._p_tth_bkgd is not None:
            lsout += "\n ttheta  IntBKGR"
            for hh_1, hh_2 in zip(self._p_tth_bkgd, self._p_int_bkgd):
                lsout += "\n {:}    {:}".format(hh_1, hh_2)
        return lsout

    def _refresh(self, tth_bkgd, int_bkgd, file_dir, file_name):
        if tth_bkgd is not None:
            self._p_tth_bkgd = tth_bkgd
        if int_bkgd is not None:
            self._p_int_bkgd = int_bkgd
        if file_dir is not None:
            self._p_file_dir = file_dir
        if file_name is not None:
            self._p_file_name = file_name
            
    def set_val(self, tth_bkgd=None, int_bkgd=None, file_dir=None, file_name=None):
        self._refresh(tth_bkgd, int_bkgd, file_dir, file_name)
        
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
        ddata={}
        fid=open(f_name,'r')
        lcontentH=fid.readlines()
        fid.close()
        lparam = [line[1:].strip() for line in lcontentH if line.startswith('#')]
        if (len(lparam) > 1):
            for line in lparam[:-1]:
                lhelp = line.strip().split()
                if (len(lhelp) > 2):
                    ddata[lhelp[0]] = [float(hh) for hh in lhelp[1:]]
                elif (len(lhelp) == 2):
                    ddata[lhelp[0]] = float(lhelp[1])
                else:
                    print("Mistake in background file '{:}' in line:\n {:}".format(finp, line))
                    print("The program is stopped.")
                    quit()
        lnames = lparam[-1].split()
        for name in lnames:
            ddata[name] = []
        
        lcontent = [line for line in lcontentH if not(line.startswith('#'))]
        for line in lcontent:
            for name, val in zip(lnames, line.strip().split()):
                ddata[name].append(val)
        
        l_tth_b, l_int_b = [], []
        i_numb = 0
        for hh_1, hh_2 in zip(ddata["ttheta"], ddata["IntBKGR"])             :
            i_numb += 1
            val_1 = float(hh_1)
            l_help = hh_2.split("(")
            if len(l_help) > 1:
                val_2 = Variable(val=float(l_help[0]), ref=True, name="IntBKGR_{:}".format(i_numb))
            else:
                val_2 = float(hh_2)
            l_tth_b.append(val_1)
            l_int_b.append(val_2)
        self.set_val(tth_bkgd=l_tth_b, int_bkgd=l_int_b)

    def save_data(self):
        l_tth_b = self.get_val("tth_bkgd")
        l_int_b = self.get_val("int_bkgd")
        
        if ((l_tth_b is None)|(l_int_b is None)):
            return
        ls_out = ["#   ttheta     IntBKGR"]
        for tth_b, int_b in zip(l_tth_b, l_int_b):
            ls_out.append(" {:}  {:}".format(tth_b, int_b))

        f_name = os.path.join(self._p_file_dir, self._p_file_name)
        fid = open(f_name, "w")
        fid.write("\n".join(ls_out))
        fid.close()

    def interpolate_by_points(self, tth):
        l_tth_b = self._p_tth_bkgd
        l_int_b = self._p_int_bkgd
        tth_b = numpy.array([1.*hh for hh in l_tth_b], dtype=float)
        int_b = numpy.array([1.*hh for hh in l_int_b], dtype=float)
        
        if l_tth_b is None:
            file_dir = self._p_file_dir 
            file_name = self._p_file_name 
            if file_name is None:
                f_inp = os.path.join(file_dir, file_name)
                self.read_data(f_inp)
                l_tth_b = self._p_tth_bkgd
                l_int_b = self._p_int_bkgd
                tth_b = numpy.array([1.*hh for hh in l_tth_b], dtype=float)
                int_b = numpy.array([1.*hh for hh in l_int_b], dtype=float)
        if l_tth_b is not None:
            int_1d = numpy.interp(tth, tth_b, int_b)
        else:
            int_1d = numpy.zeros(tth.size, dtype=float)
        return int_1d
    
    def get_variables(self):
        l_tth_b = self._p_tth_bkgd
        l_int_b = self._p_int_bkgd
        l_variable = []
        for hh_1, hh_2 in zip(l_tth_b, l_int_b):
            if isinstance(hh_1, Variable):
                l_variable.append(hh_1)
            if isinstance(hh_2, Variable):
                l_variable.append(hh_2)
        return l_variable
    
    def create_input_file(self, f_inp=None):
        if f_inp is not None:
            f_dir= os.path.dirname(f_inp)
            f_name= os.path.bathename(f_inp)
            self.set_val(file_dir=f_dir, file_name=f_name)
        f_dir = self._p_file_dir
        f_name = self._p_file_name
        f_full = os.path.join(f_dir, f_name)
        
        s_out = """#   ttheta     IntBKGR
  4.50   150
 40.00   150
 80.00   100"""
    
        fid = open(f_full, "w")
        fid.write(s_out)
        fid.close()
