"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.interpolate

#class BeamPolarization
from f_experiment.f_powder_1d.cl_setup_powder_1d import BeamPolarization
from f_common.cl_variable import *
#Description of setup class

class ResolutionPowder2D(dict):
    """
    Resoulution of the diffractometer
    """
    def __init__(self, u = 0, v = 0, w = 0.01, x = 0, y = 0):
        super(ResolutionPowder2D, self).__init__()
        self._p_u = None
        self._p_v = None
        self._p_w = None
        self._p_x = None
        self._p_y = None

        self._p_tan_th = None
        self._p_tan_th_sq = None
        self._p_cos_th = None
        self._p_hg = None
        self._p_hl = None
        self._p_hpv = None
        self._p_eta = None
        self._p_ag = None
        self._p_bg = None
        self._p_al = None
        self._p_bl = None
        
        self._refresh(u, v, w, x, y)
        
    def __repr__(self):
        lsout = """ResolutionPowder2D: 
 u: {:}\n v: {:}\n w: {:}\n x: {:}\n y: {:}""".format(self.get_val("u"),  
 self.get_val("v"), self.get_val("w"), self.get_val("x"), self.get_val("y"))
        return lsout

    def _refresh(self, u, v, w, x, y):
        if not(isinstance(u, type(None))):
            self._p_u = u
        if not(isinstance(v, type(None))):
            self._p_v = v
        if not(isinstance(w, type(None))):
            self._p_w = w
        if not(isinstance(x, type(None))):
            self._p_x = x
        if not(isinstance(y, type(None))):
            self._p_y = y
            
    def set_val(self, u=None, v=None, w=None, x=None, y=None):
        self._refresh(u, v, w, x, y)
        
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
u, v, w are coefficients to describe the Gauss part of peak shape

x, y are coefficients to describe the Lorentz part of peak shape
        """
        print(lsout)

    
    def _calc_tancos(self, th_hkl):
        """
        tth_hkl in radianas
        calculate tangenth (theta)
        """
        self._p_t_th = numpy.tan(th_hkl)
        self._p_t_th_sq = self._p_t_th**2
        res = numpy.cos(th_hkl)

        self._p_c_th = res
        self._p_ic_th = 1./res
        
    def _calc_hg(self, i_g = 0.):
        """
        ttheta in radians, could be array
        gauss size
        """
        u, v, w = 1.*self._p_u, 1.*self._p_v, 1.*self._p_w
        res_sq = (u*self._p_t_th_sq + v*self._p_t_th + w + i_g*self._p_ic_th**2)
        self._p_hg = numpy.sqrt(res_sq)
        
    def _calc_hl(self):
        """
        ttheta in radians, could be array
        lorentz site
        """
        x, y = 1.*self._p_x, 1.*self._p_y
        self._p_hl = x*self._p_t_th + y*self._p_ic_th


    def _calc_hpveta(self):
        """
        ttheta in radians, could be array
        pseudo-Voight function
        """
        hg = self._p_hg
        hl = self._p_hl

        hg_2, hl_2 = hg**2, hl**2
        hg_3, hl_3 = hg_2*hg, hl_2*hl
        hg_4, hl_4 = hg_3*hg, hl_3*hl
        hg_5, hl_5 = hg_4*hg, hl_4*hl
        c_2, c_3, c_4, c_5 = 2.69269, 2.42843, 4.47163, 0.07842
        hpv = (hg_5 + c_2*hg_4*hl + c_3*hg_3*hl_2 + 
                       c_4*hg_2*hl_3 + c_5*hg*hl_4 + hl_5)**0.2
        hh = hl*1./hpv
        self._p_hpv = hpv 
        self._p_eta = 1.36603*hh - 0.47719*hh**2 + 0.11116*hh**3

    def _calc_agbg(self):
        hpv = self._p_hpv

        self._p_ag = (2./hpv)*(numpy.log(2.)/numpy.pi)**0.5
        self._p_bg = 4*numpy.log(2)/(hpv**2)
        
    def _calc_albl(self):
        hpv = self._p_hpv
        self._p_al = 2./(numpy.pi*hpv )
        self._p_bl = 4./(hpv**2)
    
    def calc_resolution(self, tth_hkl, i_g = 0.):
        """
        Calculate parameters for tth
        tth_hkl in degrees
        """
        self._calc_tancos(0.5*tth_hkl*numpy.pi/180.)
        self._calc_hg(i_g = i_g)
        self._calc_hl()
        self._calc_hpveta()
        self._calc_agbg()
        self._calc_albl()

        a_g = self._p_ag  
        b_g = self._p_bg  
        a_l = self._p_al  
        b_l = self._p_bl 
        h_g = self._p_hg  
        h_l = self._p_hl 
        h_pv = self._p_hpv 
        eta = self._p_eta
        
        return h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l


    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_u, Variable), 
                   isinstance(self._p_v, Variable), 
                   isinstance(self._p_w, Variable), 
                   isinstance(self._p_x, Variable), 
                   isinstance(self._p_y, Variable)])
        return res        

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_u, Variable):
            l_variable.append(self._p_u)
        if isinstance(self._p_v, Variable):
            l_variable.append(self._p_v)
        if isinstance(self._p_w, Variable):
            l_variable.append(self._p_w)
        if isinstance(self._p_x, Variable):
            l_variable.append(self._p_x)
        if isinstance(self._p_y, Variable):
            l_variable.append(self._p_y)
        return l_variable

class FactorLorentzPowder2D(dict):
    """
    Lorentz Factor for one dimensional powder diffraction
    """
    def __init__(self):
        super(FactorLorentzPowder2D, self).__init__()
        dd= {}
        self.update(dd)
        
    def __repr__(self):
        lsout = """FactorLorentzPowder2D:"""
        return lsout

    def _refresh(self):
        print("'_refresh' is not introduced for FactorLorentzPD")
            
    def set_val(self):
        print("'set_val' is not introduced for FactorLorentzPD")
        
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
no parameters
        """
        print(lsout)
        
    
    def calc_f_lorentz(self, tth):
        """
        Lorentz factor
        tth should be in degrees
        """
        tth_rad = tth*numpy.pi/180.
        factor_lorentz = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
        return factor_lorentz 



class AsymmetryPowder2D(dict):
    """
    Asymmetry of the diffractometer
    """
    def __init__(self, p1 = 0., p2 = 0., p3 = 0., p4 = 0.):
        super(AsymmetryPowder2D, self).__init__()
        self._p_p1 = None
        self._p_p2 = None
        self._p_p3 = None
        self._p_p4 = None
        self._refresh(p1, p2, p3, p4)
        
    def __repr__(self):
        lsout = """AsymmetryPowder2D:\n p1: {:}\n p2: {:}\n p3: {:}
 p4: {:}""".format(self.get_val("p1"),  self.get_val("p2"),  
 self.get_val("p3"),  self.get_val("p4"))
        return lsout

    def _refresh(self, p1, p2, p3, p4):
        if not(isinstance(p1, type(None))):
            self._p_p1 = p1
        if not(isinstance(p2, type(None))):
            self._p_p2 = p2
        if not(isinstance(p3, type(None))):
            self._p_p3 = p3
        if not(isinstance(p4, type(None))):
            self._p_p4 = p4
            
    def set_val(self, p1=None, p2=None, p3=None, p4=None):
        self._refresh(p1, p2, p3, p4)
        
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
p1, p2, p3, p4 are coefficients to describe the assymetry shape for all 
               reflections like in FullProf
        """
        print(lsout)

        
    def _func_fa(self, tth):
        """
        for assymmetry correction
        """ 
        return 2*tth*numpy.exp(-tth**2)
        
    def _func_fb(self, tth):
        """
        for assymmetry correction
        """ 
        return 2.*(2.*tth**2-3.)* self._func_fa(tth)
        
    def calc_asymmetry(self, tth, tth_hkl):
        """
        Calculate asymmetry coefficients for  on the given list ttheta for 
        bragg reflections flaced on the position ttheta_hkl
        tth and tth_hkl in degrees
        
        IMPORTANT: THERE IS MISTAKE (look page 54 in FullProf Manual)
        """
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl, indexing="ij")
        np_zero = numpy.zeros(tth_2d.shape, dtype = float)
        np_one = numpy.ones(tth_2d.shape, dtype = float)
        val_1, val_2 = np_zero, np_zero
        
        
        p1, p2 = self.get_val("p1"), self.get_val("p2")
        p3, p4 = self.get_val("p3"), self.get_val("p4")
        flag_1, flag_2 = False, False
        if ((p1!= 0.)|(p3!= 0.)):
            flag_1 = True
            fa = self._func_fa(tth)
        if ((p2!= 0.)|(p4!= 0.)):
            flag_2 = True
            fb = self._func_fb(tth)
            
        flag_3, flag_4 = False, False
        if ((p1!= 0.)|(p2!= 0.)):
            if flag_1:
                val_1 += p1*fa
                flag_3 = True
            if flag_2:
                val_1 += p2*fb
                flag_3 = True
            if flag_3:
                c1 = 1./numpy.tanh(0.5*tth_hkl)
                c1_2d = numpy.meshgrid(tth, c1, indexing="ij")[1]
                val_1 *= c1_2d

        if ((p3!= 0.)|(p4!= 0.)):
            if flag_1:
                val_2 += p3*fa
                flag_4 = True
            if flag_2:
                val_2 += p4*fb
                flag_4 = True
            if flag_4:
                c2 = 1./numpy.tanh(tth_hkl)
                c2_2d = numpy.meshgrid(tth, c2, indexing="ij")[1]
                val_2 *= c2_2d

        asymmetry_2d = np_one+val_1+val_2
        return asymmetry_2d
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_p1, Variable), 
                   isinstance(self._p_p2, Variable), 
                   isinstance(self._p_p3, Variable), 
                   isinstance(self._p_p4, Variable)])
        return res        
    
    def get_variables(self):
        l_variable = []
        if isinstance(self._p_p1, Variable):
            l_variable.append(self._p_p1)
        if isinstance(self._p_p2, Variable):
            l_variable.append(self._p_p2)
        if isinstance(self._p_p3, Variable):
            l_variable.append(self._p_p3)
        if isinstance(self._p_p4, Variable):
            l_variable.append(self._p_p4)
        return l_variable
    
        
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
        
class SetupPowder2D(dict):
    """
    Class to describe characteristics of powder diffractometer
    """
    def __init__(self, label="exp", wave_length=1.4, zero_shift=0., 
                 resolution=ResolutionPowder2D(), factor_lorentz=FactorLorentzPowder2D(), 
                 asymmetry=AsymmetryPowder2D(), beam_polarization=BeamPolarization(),
                 background=BackgroundPowder2D()):
        super(SetupPowder2D, self).__init__()
        self._p_label = None
        self._p_wave_length = None
        self._p_zero_shift = None
        self._p_resolution = None
        self._p_factor_lorentz = None
        self._p_asymmetry = None
        self._p_beam_polarization = None
        self._p_background = None
        
        self._refresh(label, wave_length, zero_shift, resolution, 
                      factor_lorentz, asymmetry, beam_polarization, background)

    def __repr__(self):
        lsout = """SetupPowder2D:\n label: {:}\n wave_length: {:}
 zero_shift: {:}\n{:}\n{:}\n{:}\n{:}\n{:}""".format(self._p_label, 
 self._p_wave_length, self._p_zero_shift, self._p_resolution, 
 self._p_factor_lorentz, self._p_asymmetry, self._p_beam_polarization, 
 self._p_background)
        return lsout

    def _refresh(self, label, wave_length, zero_shift, resolution, 
                 factor_lorentz, asymmetry, beam_polarization, background):
        if label is not None:
            self._p_label = label
        if wave_length is not None:
            self._p_wave_length = wave_length
        if zero_shift is not None:
            self._p_zero_shift = zero_shift
        if resolution is not None:
            self._p_resolution = resolution
        if factor_lorentz is not None:
            self._p_factor_lorentz = factor_lorentz
        if asymmetry is not None:
            self._p_asymmetry = asymmetry
        if beam_polarization is not None:
            self._p_beam_polarization = beam_polarization
        if background is not None:
            self._p_background = background
            
    def set_val(self, label=None, wave_length=None, zero_shift=None, 
                resolution=None, factor_lorentz=None, asymmetry=None, 
                beam_polarization=None, background=None):
        
        self._refresh(label, wave_length, zero_shift, resolution, 
                      factor_lorentz, asymmetry, beam_polarization, background)
        
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
label is just to make labe
wave_length is to describe wave_length in angstrems
zero_shift is to describe zeroshift in degrees

resolution is a class to describe resolution of powder diffractometer
factor_lorentz is a class to describe factor Lorentz
asymmetry is a class to descibe the asymmetry
beam_polarization is a class to describe beam polarization
background  is Background class
        """
        print(lsout)

    def _gauss_pd(self, tth_3d):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self._p_ag, self._p_bg
        val_1 = bg*tth_3d**2
        val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
        self._p_gauss_pd = ag*val_2
        
    def _lor_pd(self, tth_3d):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self._p_al, self._p_bl
        self._p_lor_pd = al*1./(1.+bl*tth_3d**2)
    
    def calc_shape_profile(self, tth, phi, tth_hkl, i_g=0.):
        """
        calculate profile in the range ttheta for reflections placed on 
        ttheta_hkl with i_g parameter by defoult equal to zero
        
        tth, phi, tth_hkl in degrees

        """
        
        resolution = self._p_resolution
        
        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
                                                  tth_hkl, i_g=i_g)

        
        tth_3d, phi_3d, tth_hkl_3d = numpy.meshgrid(tth, phi, tth_hkl, indexing="ij")

        self._p_ag = numpy.meshgrid(tth, phi, a_g, indexing="ij")[2]
        self._p_bg = numpy.meshgrid(tth, phi, b_g, indexing="ij")[2]
        self._p_al = numpy.meshgrid(tth, phi, a_l, indexing="ij")[2]
        self._p_bl = numpy.meshgrid(tth, phi, b_l, indexing="ij")[2]
        eta_3d = numpy.meshgrid(tth, phi, eta, indexing="ij")[2]
        self._p_eta = eta_3d 

        self._gauss_pd(tth_3d-tth_hkl_3d)
        self._lor_pd(tth_3d-tth_hkl_3d)
        g_pd_3d = self._p_gauss_pd 
        l_pd_3d = self._p_lor_pd
        
        profile_3d = eta_3d * l_pd_3d + (1.-eta_3d) * g_pd_3d
        
        return profile_3d
    
    def calc_profile(self, tth, phi, tth_hkl, i_g):
        """
        tth and tth_hkl in degrees
        """
        zero_shift = 1*self._p_zero_shift
        tth_zs = tth-zero_shift
        np_shape_3d = self.calc_shape_profile(tth_zs, phi, tth_hkl, i_g=i_g)
        asymmetry = self.get_val("asymmetry")
        #dimension (tth, hkl)
        np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl)
        
        #dimension (tth, phi, hkl)
        np_ass_3d = np_ass_2d[:, numpy.newaxis,:]*numpy.ones(phi.size, dtype=float)[numpy.newaxis, :, numpy.newaxis]
        
        factor_lorentz = self.get_val("factor_lorentz")
        np_lor_1d = factor_lorentz.calc_f_lorentz(tth_zs)
        
        
        np_lor_3d = numpy.meshgrid(np_lor_1d, phi, tth_hkl, indexing="ij")[0]
        
        
        profile_3d = np_shape_3d*np_ass_3d*np_lor_3d
        return profile_3d 
    

    def calc_background(self, tth, phi):
        """
        estimates background points on the given ttheta positions
        """                   
        background = self._p_background
        int_bkgd = background.interpolate_by_points(tth, phi)
        return int_bkgd 
    
    def is_variable(self):
        """
        without extinction
        """
        beam_polarization = self.get_val("beam_polarization")
        resolution = self.get_val("resolution")
        asymmetry = self.get_val("asymmetry")
        res = any([isinstance(self._p_zero_shift, Variable), 
                   beam_polarization.is_variable(), 
                   resolution.is_variable(), 
                   asymmetry.is_variable()])
        return res   
    
    def get_variables(self):
        l_variable = []
        if isinstance(self._p_zero_shift, Variable):
            l_variable.append(self._p_zero_shift)

        background = self.get_val("background")
        l_var = background.get_variables()
        l_variable.extend(l_var)

        beam_polarization = self.get_val("beam_polarization")
        l_var = beam_polarization.get_variables()
        l_variable.extend(l_var)

        resolution = self.get_val("resolution")
        l_var = resolution.get_variables()
        l_variable.extend(l_var)

        asymmetry = self.get_val("asymmetry")
        l_var = asymmetry.get_variables()
        l_variable.extend(l_var)

        return l_variable
