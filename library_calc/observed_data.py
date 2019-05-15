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
                 orientation=None):
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
                 orientation)

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
                 orientation):
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
            
    def set_val(self,  h=None, k=None, l=None, flip_ratio=None, 
                 sflip_ratio=None, wave_length=None, field=None, 
                 orientation=None):
        self._refresh(h, k, l, flip_ratio, sflip_ratio, wave_length, field, 
                 orientation)
        
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
        """
        print(lsout)
    
    def read_data(self, finp):
        """
        read file from file
        """
        self._p_file_dir = os.path.dirname(finp)
        self._p_file_name = os.path.basename(finp)
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



class ObservedDataPowder1D(dict):
    """
    Containt the experimental data
    """
    def __init__(self, tth_exp=None, int_u_exp=None, sint_u_exp=None, 
                 int_d_exp=None, sint_d_exp=None, field=None, wave_length=None):
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
                      field, wave_length)

    def __repr__(self):
        ls_out = """ObservedDataPowder1D:\n file_dir: {:}
 file_name: {:}""".format(self._p_file_dir, self._p_file_name)
        if self._p_tth is not None:
            ls_out += "\n tth range: {:} --- {:} ({:} points)".format(
                    self._p_tth.min(), self._p_tth.max(), self._p_tth.size)
        ls_out += "\n field: {:}".format(self._p_field)
        ls_out += "\n wave_length: {:}".format(self._p_wave_length)
        return ls_out

    def _refresh(self, tth_exp, int_u_exp, sint_u_exp, int_d_exp, sint_d_exp, 
                 field, wave_length):
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
        if field is not None:
            self._p_field = field
        if wave_length is not None:
            self._p_wave_length = wave_length
        if flag:
            self.exclude_data()
            
    def set_val(self, tth_exp=None, int_u_exp=None, sint_u_exp=None, 
                int_d_exp=None, sint_d_exp=None, field=None, wave_length=None):
        self._refresh(tth_exp, int_u_exp, sint_u_exp, int_d_exp, sint_d_exp, 
                      field, wave_length)
        
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
    
    def read_data(self, finp):
        """
        read file from file
        """
        self._p_file_dir = os.path.dirname(finp)
        self._p_file_name = os.path.basename(finp)
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


class ObservedDataPowder2D(dict):
    """
    Containt the experimental data
    """
    def __init__(self, tth_exp=None, phi_exp=None, int_u_exp=None, 
                 sint_u_exp=None, int_d_exp=None, 
                 sint_d_exp=None, tth_min=None, tth_max=None, phi_min=None, 
                 phi_max=None, field=None, wave_length=None):
        super(ObservedDataPowder2D, self).__init__()
        self._p_tth_exp = None
        self._p_phi_exp = None
        self._p_int_u_exp = None
        self._p_sint_u_exp = None
        self._p_int_d_exp = None
        self._p_sint_d_exp = None
        
        self._p_tth_min = None
        self._p_tth_max = None
        self._p_phi_min = None
        self._p_phi_max = None

        self._p_tth = None
        self._p_phi = None
        self._p_int_u = None
        self._p_sint_u = None
        self._p_int_d = None
        self._p_sint_d = None

        self._p_file_dir = None
        self._p_file_name = None        
        
        self._p_field = None
        self._p_wave_length = None
        
        self._refresh(tth_exp, phi_exp, int_u_exp, sint_u_exp, int_d_exp, 
                      sint_d_exp, tth_min, tth_max, 
                      phi_min, phi_max, field, wave_length)

    def __repr__(self):
        ls_out = """ObservedDataPowder2D:\n file_dir: {:}
 file_name: {:}""".format(self._p_file_dir, self._p_file_name)
        if self._p_tth is not None:
            ls_out += "\n tth range: {:} --- {:} ({:} points)".format(
                    self._p_tth.min(), self._p_tth.max(), self._p_tth.size)
        if self._p_phi is not None:
            ls_out += "\n phi range: {:} --- {:} ({:} points)".format(
                    self._p_phi.min(), self._p_phi.max(), self._p_phi.size)
        ls_out += "\n field: {:}".format(self._p_field)
        ls_out += "\n wave_length: {:}".format(self._p_wave_length)
        
        return ls_out

    def _refresh(self, tth_exp, phi_exp, int_u_exp, sint_u_exp, int_d_exp, 
                 sint_d_exp, tth_min, 
                 tth_max, phi_min, phi_max, field, wave_length):
        flag = any([(hh is not None) for hh in [tth_exp, phi_exp, int_u_exp, 
                                                sint_u_exp, int_d_exp, 
                                                sint_d_exp, tth_min, tth_max, 
                                                phi_min, phi_max]])
        if tth_exp is not None:
            self._p_tth_exp = tth_exp
        if phi_exp is not None:
            self._p_phi_exp = phi_exp
        if int_u_exp is not None:
            self._p_int_u_exp = int_u_exp
        if sint_u_exp is not None:
            self._p_sint_u_exp = sint_u_exp
        if int_d_exp is not None:
            self._p_int_d_exp = int_d_exp
        if sint_d_exp is not None:
            self._p_sint_d_exp = sint_d_exp
        if field is not None:
            self._p_field = field
        if wave_length is not None:
            self._p_wave_length = wave_length
        if tth_min is not None:
            self._p_tth_min = tth_min
        if tth_max is not None:
            self._p_tth_max = tth_max
        if phi_min is not None:
            self._p_phi_min = phi_min
        if phi_max is not None:
            self._p_phi_max = phi_max
        if flag:
            self.exclude_data()
            
    def set_val(self, tth_exp=None, phi_exp=None, int_u_exp=None, 
                sint_u_exp=None, int_d_exp=None, sint_d_exp=None, 
                tth_min=None, tth_max=None, phi_min=None, phi_max=None, 
                field=None, wave_length=None):
        self._refresh(tth_exp, phi_exp, int_u_exp, sint_u_exp, int_d_exp, 
                      sint_d_exp, tth_min, tth_max, 
                      phi_min, phi_max, field, wave_length)
        
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
phi is 1D array of phi in degrees

int_u, sint_u are 2D array of intensity with errorbars at flipper postion 'up'
int_d, sint_d are 2D array of intensity with errorbars at flipper postion 'down'

field is the magnetic field along z axis 
wave_length is the neutron wave_length in Angstrem
        """
        print(lsout)

    def exclude_data(self):
        tth_exp = self._p_tth_exp
        phi_exp = self._p_phi_exp
        int_u_exp = self._p_int_u_exp
        sint_u_exp = self._p_sint_u_exp
        int_d_exp = self._p_int_d_exp
        sint_d_exp = self._p_sint_d_exp
        flag_tth = numpy.ones(tth_exp.size, dtype=bool)
        flag_phi = numpy.ones(phi_exp.size, dtype=bool)

        tth_min = self._p_tth_min 
        tth_max = self._p_tth_max 
        phi_min = self._p_phi_min 
        phi_max = self._p_phi_max 
        
        cond_tth = (tth_min is not None)&(tth_max is not None)
        cond_phi = (phi_min is not None)&(phi_max is not None)

        if cond_tth:
            flag_1 = tth_exp >= tth_min
            flag_2 = tth_max >= tth_exp
            flag_tth = numpy.logical_and(flag_1, flag_2)

        if cond_phi:
            flag_1 = phi_exp >= phi_min
            flag_2 = phi_max >= phi_exp
            flag_phi = numpy.logical_and(flag_1, flag_2)
        
        tth = tth_exp[flag_tth]
        phi = phi_exp[flag_phi]
        int_u = int_u_exp[flag_tth][:, flag_phi]
        sint_u = sint_u_exp[flag_tth][:, flag_phi]
        int_d = int_d_exp[flag_tth][:, flag_phi]
        sint_d = sint_d_exp[flag_tth][:, flag_phi]
    
        self._p_tth = tth
        self._p_phi = phi
        self._p_int_u = int_u 
        self._p_sint_u = sint_u 
        self._p_int_d = int_d 
        self._p_sint_d = sint_d 

    def read_data(self, finp):
        """
        read file from file
        """
        self._p_file_dir = os.path.dirname(finp)
        self._p_file_name = os.path.basename(finp)
        ddata = {}
        fid = open(finp,'r')
        lcontentH = fid.readlines()
        fid.close()
        lparam = [line[1:].strip() for line in lcontentH if line.startswith('#')]
        if (len(lparam) > 0):
            for line in lparam:
                lhelp = line.strip().split()
                if (len(lhelp) > 2):
                    ddata[lhelp[0]] = [float(hh) for hh in lhelp[1:]]
                elif (len(lhelp) == 2):
                    ddata[lhelp[0]] = float(lhelp[1])
                else:
                    print("Mistake in experimental file '{:}' in line:\n {:}".format(finp, line))
                    print("The program is stopped.")
                    quit()

        lcontent = [line for line in lcontentH if not line.startswith('#')]

        lmat, mat = [], []
        for hh in lcontent:
            if hh.strip() == "":
                if mat != []:
                    lmat.append(mat)
                    mat = []
            else:
                mat.append(hh)
        if mat != []:
            lmat.append(mat)
            mat = []
        
        ll_int, l_ang1, l_ang2 = [], [], []
        for mat in lmat:
            ang1 = [float(hh) for hh in mat[0].split()[1:]]
            l_int, ang2 = [], []
            for hh in mat[1:]:
                lhelp = hh.split()
                ang2.append(float(lhelp[0]))
                l_int.append([float(hh2) if hh2 != "None" else None for hh2 in lhelp[1:]])
            ll_int.append(l_int)
            l_ang1.append(ang1)
            l_ang2.append(ang2)

        field = ddata["field"]
        wave_length = ddata["wave_length"]
        tth = numpy.array(l_ang1[0], dtype=float)
        phi = numpy.array(l_ang2[0], dtype=float)
        int_u = numpy.array(ll_int[0], dtype=float).transpose()
        sint_u = numpy.array(ll_int[1], dtype=float).transpose()
        int_d = numpy.array(ll_int[2], dtype=float).transpose()
        sint_d = numpy.array(ll_int[3], dtype=float).transpose()
        
        self.set_val(tth_exp=tth, phi_exp=phi, int_u_exp=int_u, 
                     sint_u_exp=sint_u, int_d_exp=int_d, 
                     sint_d_exp=sint_d, field=field, wave_length=wave_length)
        
        if self._p_tth_min is None:
            self._p_tth_min = tth.min()
        if self._p_tth_max is None:
            self._p_tth_max = tth.max()

        if self._p_phi_min is None:
            self._p_phi_min = phi.min()
        if self._p_phi_max is None:
            self._p_phi_max = phi.max()
            

if (__name__ == "__main__"):
  pass

