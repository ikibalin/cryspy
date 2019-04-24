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
    def __init__(self, tth=None, int_u=None, sint_u=None, int_d=None, 
                 sint_d=None, field=None, wavelength=None):
        super(ObservedDataPowder1D, self).__init__()
        self._p_tth = None
        self._p_int_u = None
        self._p_sint_u = None
        self._p_int_d = None
        self._p_sint_d = None
        
        self._p_field = None
        self._p_wavelength = None
        
        self._refresh(tth, int_u, sint_u, int_d, sint_d, field, wavelength)

    def __repr__(self):
        lsout = """Observed data:"""
        return lsout

    def _refresh(self, tth, int_u, sint_u, int_d, sint_d, field, wavelength):
        if not(isinstance(tth, type(None))):
            self._p_tth = tth
        if not(isinstance(int_u, type(None))):
            self._p_int_u = int_u
        if not(isinstance(sint_u, type(None))):
            self._p_sint_u = sint_u
        if not(isinstance(int_d, type(None))):
            self._p_int_d = int_d
        if not(isinstance(sint_d, type(None))):
            self._p_sint_d = sint_d
        if not(isinstance(field, type(None))):
            self._p_field = field
        if not(isinstance(wavelength, type(None))):
            self._p_wavelength = wavelength
            
    def set_val(self, tth=None, int_u=None, sint_u=None, int_d=None, 
                 sint_d=None, field=None, wavelength=None):
        self._refresh(tth, int_u, sint_u, int_d, sint_d, field, wavelength)
        
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
wavelength is the neutron wavelength
        """
        print(lsout)
    
    def read_data(self, finp):
        """
        read file from file
        """
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
        lcontent = [line for line in lcontentH if line[0]!='#']
        for line in lcontent:
            for name, val in zip(lnames, line.strip().split()):
                ddata[name].append(val)
        field = ddata["field"]
        wavelength = ddata["wavelength"]
        tth = numpy.array(ddata["ttheta"], dtype=float)
        int_u = numpy.array(ddata["IntUP"], dtype=float)
        sint_u = numpy.array(ddata["sIntUP"], dtype=float)
        int_d = numpy.array(ddata["IntDOWN"], dtype=float)
        sint_d = numpy.array(ddata["sIntDOWN"], dtype=float)
        self.set_val(tth=tth, int_u=int_u, sint_u=sint_u, int_d=int_d, 
                     sint_d=sint_d, field=field, wavelength=wavelength)

        
if (__name__ == "__main__"):
  pass

