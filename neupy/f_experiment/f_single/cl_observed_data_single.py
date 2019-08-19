"""
define classes to describe observed data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy


class ObservedDataSingle(dict):
    """
    Containt the experimental data
    """
    def __init__(self, h=None, k=None, l=None, flip_ratio=None, 
                 sflip_ratio=None, wave_length=None, field=None, 
                 orientation=None, file_dir=None, file_name=None):
        super(ObservedDataSingle, self).__init__()
        self._p_h = None
        self._p_k = None
        self._p_l = None
        self._p_flip_ratio = None
        self._p_sflip_ratio = None
        
        self._p_file_dir = None
        self._p_file_name = None
        
        self._p_wave_length = None
        self._p_field = None
        self._p_orientation = None
        
        self._refresh(h, k, l, flip_ratio, sflip_ratio, wave_length, field, 
                 orientation, file_dir, file_name)

    def __repr__(self):
        ls_out = """ObservedDataSingle:\n file_dir: {:}
 file_name: {:}""".format(self._p_file_dir, self._p_file_name)
        if self._p_h is not None:
            ls_out += "\n h range: {:} --- {:}".format(
                    self._p_h.min(), self._p_h.max())
        if self._p_k is not None:
            ls_out += "\n k range: {:} --- {:}".format(
                    self._p_k.min(), self._p_k.max())
        if self._p_l is not None:
            ls_out += "\n l range: {:} --- {:}".format(
                    self._p_l.min(), self._p_l.max())
            ls_out += "\n number of reflections: {:}".format(
                    self._p_l.size)
        ls_out += "\n field: {:}".format(self._p_field)
        ls_out += "\n wave_length: {:}".format(self._p_wave_length)
        ls_out += "\n orientation: {:}".format(self._p_orientation)
        return ls_out

    def _refresh(self, h, k, l, flip_ratio, sflip_ratio, wave_length, field, 
                 orientation, file_dir, file_name):
        if h is not None:
            self._p_h = h
        if k is not None:
            self._p_k = k
        if l is not None:
            self._p_l = l
        if flip_ratio is not None:
            self._p_flip_ratio = flip_ratio
        if sflip_ratio is not None:
            self._p_sflip_ratio = sflip_ratio
        if wave_length is not None:
            self._p_wave_length = wave_length
        if field is not None:
            self._p_field = field
        if orientation is not None:
            self._p_orientation = orientation
        if file_dir is not None:
            self._p_file_dir = file_dir
        if file_name is not None:
            self._p_file_name = file_name
            
    def set_val(self,  h=None, k=None, l=None, flip_ratio=None, 
                 sflip_ratio=None, wave_length=None, field=None, 
                 orientation=None, file_dir=None, file_name=None):
        self._refresh(h, k, l, flip_ratio, sflip_ratio, wave_length, field, 
                 orientation, file_dir, file_name)
        
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
 h, k, l are the Miller indices
 flip_ratio is the flipping ratio
 sflip_ratio is the errorbar of flipping ratio
 wave_length is the neutron wave length in the Angtrems
 field is the magnetic field along z axis 
 orientation is the tranbsformation matrix from local coordinate system to global (matrix U)
 file_dir is directory of file with measured data
 file_name is basename of file with measured data
 
        """
        print(lsout)


    def read_data(self):
        """
        read file from file
        """
        finp = os.path.join(self._p_file_dir, self._p_file_name)
        ddata = {}
        if not(os.path.isfile(finp)):
            return
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
        wave_length = ddata["wave_length"]
        l_o = ddata["orientation"]
        orientation = numpy.array([[l_o[0], l_o[1], l_o[2]], 
                                   [l_o[3], l_o[4], l_o[5]], 
                                   [l_o[6], l_o[7], l_o[8]]], dtype=float)

        h = numpy.array(ddata["h"], dtype=int)
        k = numpy.array(ddata["k"], dtype=int)
        l = numpy.array(ddata["l"], dtype=int)
        flip_ratio = numpy.array(ddata["FR"], dtype=float)
        sflip_ratio = numpy.array(ddata["sFR"], dtype=float)
        self.set_val(h=h, k=k, l=l, flip_ratio=flip_ratio, 
                 sflip_ratio=sflip_ratio, wave_length=wave_length, field=field, 
                 orientation=orientation)

    def create_input_file(self, f_inp=None):
        if f_inp is not None:
            f_dir= os.path.dirname(f_inp)
            f_name= os.path.bathename(f_inp)
            self.set_val(file_dir=f_dir, file_name=f_name)
        f_dir = self._p_file_dir
        f_name = self._p_file_name
        f_full = os.path.join(f_dir, f_name)
        
        s_out = """#wave_length 1.40 
#field  1.000
#orientation 0.6468462   -0.6860854   0.3300297 0.2141139   -0.2557555   -0.9343334 0.7319312   0.6810804    -0.0183183 
#   h    k    l        FR       sFR
    0    0    8   0.64545   0.01329 """
    
        fid = open(f_full, "w")
        fid.write(s_out)
        fid.close()
        
        

if (__name__ == "__main__"):
  pass

