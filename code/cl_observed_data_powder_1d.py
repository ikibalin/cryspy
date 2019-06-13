"""
define classes to describe observed data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy


class ObservedDataPowder1D(dict):
    """
    Containt the experimental data
    """
    def __init__(self, tth_exp=None, int_u_exp=None, sint_u_exp=None, 
                 int_d_exp=None, sint_d_exp=None, tth_min=None, tth_max=None,
                 field=None, wave_length=None, file_dir=None, file_name=None):
        super(ObservedDataPowder1D, self).__init__()
        self._p_tth_exp = None
        self._p_int_u_exp = None
        self._p_sint_u_exp = None
        self._p_int_d_exp = None
        self._p_sint_d_exp = None

        self._p_tth_min = None
        self._p_tth_max = None


        self._p_tth = None
        self._p_int_u = None
        self._p_sint_u = None
        self._p_int_d = None
        self._p_sint_d = None

        self._p_file_dir = None
        self._p_file_name = None
        
        self._p_field = None
        self._p_wave_length = None
        
        self._refresh(tth_exp, int_u_exp, sint_u_exp, int_d_exp, sint_d_exp, 
                      tth_min, tth_max, field, wave_length, file_dir, file_name)

    def __repr__(self):
        ls_out = """ObservedDataPowder1D:\n file_dir: {:}
 file_name: {:}""".format(self._p_file_dir, self._p_file_name)
        if self._p_tth is not None:
            ls_out += "\n tth range: {:} --- {:} ({:} points)".format(
                    self._p_tth.min(), self._p_tth.max(), self._p_tth.size)
            ls_out += "\n tth_min: {:}\n tth_max: {:}".format(self._p_tth_min, 
                                   self._p_tth_max)
        ls_out += "\n field: {:}".format(self._p_field)
        ls_out += "\n wave_length: {:}".format(self._p_wave_length)
        return ls_out

    def _refresh(self, tth_exp, int_u_exp, sint_u_exp, int_d_exp, sint_d_exp, 
                 tth_min, tth_max, field, wave_length, file_dir, file_name):
        flag = any([(hh is not None) for hh in [tth_exp, int_u_exp, sint_u_exp, 
                                                int_d_exp, sint_d_exp]])
        if tth_exp is not None:
            self._p_tth_exp = tth_exp
        if int_u_exp is not None:
            self._p_int_u_exp = int_u_exp
        if sint_u_exp is not None:
            self._p_sint_u_exp = sint_u_exp
        if int_d_exp is not None:
            self._p_int_d_exp = int_d_exp
        if sint_d_exp is not None:
            self._p_sint_d_exp = sint_d_exp
        if tth_min is not None:
            self._p_tth_min = tth_min
        if tth_max is not None:
            self._p_tth_max = tth_max
        if field is not None:
            self._p_field = field
        if wave_length is not None:
            self._p_wave_length = wave_length
        if file_dir is not None:
            self._p_file_dir = file_dir
        if file_name is not None:
            self._p_file_name = file_name
        if flag:
            self.exclude_data()
            
    def set_val(self, tth_exp=None, int_u_exp=None, sint_u_exp=None, 
                int_d_exp=None, sint_d_exp=None, tth_min=None, tth_max=None, 
                field=None, wave_length=None, file_dir=None, file_name=None):
        self._refresh(tth_exp, int_u_exp, sint_u_exp, int_d_exp, sint_d_exp, 
                      tth_min, tth_max, field, wave_length, file_dir, file_name)
        
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
tth is 1D array of ttheta in degrees 
int_u, sint_u are 1D array of intensity with errorbars at flipper postion 'up'
int_d, sint_d are 1D array of intensity with errorbars at flipper postion 'down'

field is the magnetic field along z axis 
wave_length is the neutron wave_length

tth_min, tth_max -- range of diffraction angle for calculations

file_dir
file_name
        """
        print(lsout)
        
    def exclude_data(self):
        tth_exp = self._p_tth_exp
        int_u_exp = self._p_int_u_exp
        sint_u_exp = self._p_sint_u_exp
        int_d_exp = self._p_int_d_exp
        sint_d_exp = self._p_sint_d_exp
        flag_tth = numpy.ones(tth_exp.size, dtype=bool)

        tth_min = self._p_tth_min 
        tth_max = self._p_tth_max 
        
        cond_tth = (tth_min is not None)&(tth_max is not None)

        if cond_tth:
            flag_1 = tth_exp >= tth_min
            flag_2 = tth_max >= tth_exp
            flag_tth = numpy.logical_and(flag_1, flag_2)
        
        tth = tth_exp[flag_tth]
        int_u = int_u_exp[flag_tth]
        sint_u = sint_u_exp[flag_tth]
        int_d = int_d_exp[flag_tth]
        sint_d = sint_d_exp[flag_tth]
    
        self._p_tth = tth
        self._p_int_u = int_u 
        self._p_sint_u = sint_u 
        self._p_int_d = int_d 
        self._p_sint_d = sint_d 
    

    def read_data(self):
        finp = os.path.join(self._p_file_dir, self._p_file_name)
        ddata = {}
        fid = open(finp,'r')
        lcontentH = fid.readlines()
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
                    print("Mistake in experimental file '{:}' in line:\n {:}".format(finp, line))
                    print("The program is stopped.")
                    quit()
        lnames = lparam[-1].split()
        for name in lnames:
            ddata[name] = []
        lcontent = [line for line in lcontentH if not(line.startswith('#'))]
        for line in lcontent:
            for name, val in zip(lnames, line.strip().split()):
                ddata[name].append(val)
        field = ddata["field"]
        wave_length = ddata["wave_length"]
        tth = numpy.array(ddata["ttheta"], dtype=float)
        int_u = numpy.array(ddata["IntUP"], dtype=float)
        sint_u = numpy.array(ddata["sIntUP"], dtype=float)
        int_d = numpy.array(ddata["IntDOWN"], dtype=float)
        sint_d = numpy.array(ddata["sIntDOWN"], dtype=float)
        self.set_val(tth_exp=tth, int_u_exp=int_u, sint_u_exp=sint_u, 
                     int_d_exp=int_d, sint_d_exp=sint_d, field=field, 
                     wave_length=wave_length)
        
        if self._p_tth_min is None:
            self._p_tth_min = tth.min()
        if self._p_tth_max is None:
            self._p_tth_max = tth.max()

    def create_input_file(self, f_inp=None):
        if f_inp is not None:
            f_dir= os.path.dirname(f_inp)
            f_name= os.path.bathename(f_inp)
            self.set_val(file_dir=f_dir, file_name=f_name)
        f_dir = self._p_file_dir
        f_name = self._p_file_name
        f_full = os.path.join(f_dir, f_name)
        
        s_out = """#wave_length 0.84
#field 1.0
#   ttheta     IntUP    sIntUP   IntDOWN  sIntDOWN
 4.00   465.80000   128.97000   301.88000   129.30000
 4.20   323.78000   118.22000   206.06000   120.00000 """
    
        fid = open(f_full, "w")
        fid.write(s_out)
        fid.close()
        

if (__name__ == "__main__"):
  pass

