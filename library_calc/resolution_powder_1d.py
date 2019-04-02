"""
define resolution for the powder diffractometer along ttheta
"""
__author__ = 'ikibalin'
__version__ = "2019_03_29"
import numpy


class ResolutionPD(dict):
    """
    Resoulution of the diffractometer
    """
    def __init__(self, u = 0, v = 0, w = 0.01, i_g = 0, x = 0, y = 0):
        super(ResolutionPD, self).__init__()
        self._p_tan_tth = None
        self._p_tan_tth_sq = None
        self._p_cos_tth = None
        self._p_hg = None
        self._p_hl = None
        self._p_hpv = None
        self._p_eta = None
        self._p_ag = None
        self._p_bg = None
        self._p_al = None
        self._p_bl = None

        dd= {"u":u, "v":v, "w":w, "i_g":i_g, "x":x, "y":y}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Resolution: 
 U {:}\n V {:}\n W {:}\n Ig {:}\n X {:}\n Y {:}""".format(self["u"],  
 self["v"], self["w"], self["i_g"], self["x"], self["y"])
        return lsout
    
    def _calc_tancos(self, tth_hkl):
        self._p_t_tth = numpy.tan(tth_hkl)
        self._p_t_tth_sq = self._p_t_tth**2
        res = numpy.cos(tth_hkl)
        self._p_c_tth = res
        self._p_ic_tth = 1./res
        
    def _calc_hg(self):
        """
        ttheta in radians, could be array
        gauss size
        """
        res_sq = (self["u"]*self._p_t_tth_sq + self["v"]*self._p_t_tth + 
                  self["w"]*self._p_c_tth + self["i_g"]*self._p_ic_tth**2)
        self._p_hg = numpy.sqrt(res_sq)
        
    def _calc_hl(self):
        """
        ttheta in radians, could be array
        lorentz site
        """
        self._p_hl = self["x"]*self._p_t_tth + self["y"]*self._p_ic_tth


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
    
    def calc_model(self, tth_hkl, u = None, v = None, w = None, i_g = None, 
                    x = None, y = None):
        """
        Calculate parameters for tth
        """
        if u != None:
            self["u"] = u
        if v != None:
            self["v"] = v
        if w != None:
            self["w"] = w
        if i_g != None:
            self["i_g"] = i_g
        if x != None:
            self["x"] = x
        if y != None:
            self["y"] = y
            
        self._calc_tancos(tth_hkl)
        self._calc_hg()
        self._calc_hl()
        self._calc_hpveta()
        self._calc_agbg()
        self._calc_albl()
        
        d_out = dict(a_g = self._p_ag, b_g = self._p_bg, 
                     a_l = self._p_al, b_l = self._p_bl,
                     h_g = self._p_hg, h_l = self._p_hl,
                     h_pv = self._p_hpv, eta = self._p_eta,
                     tth_hkl = tth_hkl)
        self.update(d_out)
        return 


class AsymmetryPD(dict):
    """
    Asymmetry of the diffractometer
    """
    def __init__(self, p1 = 0, p2 = 0, p3 = 0, p4 = 0):
        super(AsymmetryPD, self).__init__()
        dd= {"p1":p1, "p2":p2, "p3":p3, "p4":p4}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Asymmetry: \n p1: {:}\n p2: {:}\n p3: {:}
 p4: {:}""".format(self["p1"],  self["p2"],  self["p3"],  self["p4"])
        return lsout
        
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
        
    def calc_model(self, tth, tth_hkl, p1 = None, p2 = None, p3 = None, p4 = None):
        if p1 != None:
            self["p1"] = p1
        if p2 != None:
            self["p2"] = p2
        if p3 != None:
            self["p3"] = p3
        if p4 != None:
            self["p4"] = p4

        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl)
        np_zero = numpy.zeros(tth_2d.shape, dtype = float)
        np_one = numpy.ones(tth_2d.shape, dtype = float)
        val_1, val_2 = np_zero, np_zero

        p1, p2, p3, p4 = self["p1"], self["p2"], self["p3"], self["p4"]
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
                c1_2d = numpy.meshgrid(tth, c1)[1]
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
                c2_2d = numpy.meshgrid(tth, c2)[1]
                val_2 *= c2_2d

        res_2d = np_one+val_1+val_2
        d_out = dict(asymmetry = res_2d, tth = tth, tth_hkl = tth_hkl)
        self.update(d_out)
        return




class PeakShapePD(dict):
    """
    Shape of the peaks for powder diffraction measurements
    """
    def __init__(self, resolution = ResolutionPD()):
        super(PeakShapePD, self).__init__()
        dd= {"resolution":resolution}
        self._p_ag = None
        self._p_bg = None
        self._p_al = None
        self._p_bl = None
        self._eta = None
        self._p_gauss_pd = None
        self._p_lor_pd = None
        self.update(dd)
        
    def __repr__(self):
        lsout = """Shape of the profile: \n {:}""".format(self["resolution"])
        return lsout
        
    def _gauss_pd(self, tth):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self._p_ag, self._p_bg
        self._p_gauss_pd = ag*numpy.exp(-bg*tth**2)
        
    def _lor_pd(self, tth):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self._p_al, self._p_bl
        self._p_lor_pd = al*1./(1.+bl*tth**2)
    
    def calc_model(self, tth, tth_hkl, resolution = None):
        """
        pseudo voight function
        """
        if resolution != None :
            self["resolution"] = resolution
            
        resolution = self["resolution"]            
        resolution.calc_model(tth_hkl)
        
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl)
        
        self._p_ag = resolution["a_g"]
        self._p_bg = resolution["b_g"]
        self._p_al = resolution["a_l"]
        self._p_bl = resolution["b_l"]
        eta_2d = resolution["eta"]
        self._p_eta = eta_2d 

        g_pd_2d = self._gauss_pd(tth_2d-tth_hkl_2d)
        l_pd_2d = self._lor_pd(tth_2d-tth_hkl_2d)
        
        res_2d = eta_2d * l_pd_2d + (1.-eta_2d) * g_pd_2d
        
        d_out = dict(profile = res_2d, 
                     tth = tth, tth_hkl = tth_hkl)
        self.update(d_out)


class FactorLorentzPD(dict):
    """
    Lorentz Factor for one dimensional powder diffraction
    """
    def __init__(self):
        super(FactorLorentzPD, self).__init__()
        dd= {}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Lorentz factor for one dimensional powder diffraction. """
        return lsout
        
    
    def calc_model(self, tth):
        """
        Lorentz factor
        tth should be in radians
        """
        factor_lorentz = 1./(numpy.sin(tth)*numpy.sin(0.5*tth))
        d_out = dict(tth = tth, factor_lorentz = factor_lorentz)
        self.update(d_out)





class PeakProfilePD(dict):
    """
    Peak Profile
    """
    def __init__(self, peak_shape = PeakShapePD(), asymmetry = AsymmetryPD(),
                 factor_lorentz = FactorLorentzPD(), i_g = 0.):
        super(PeakProfilePD, self).__init__()
        dd= {"peak_shape": peak_shape, "asymmetry": asymmetry, 
             "factor_lorentz": factor_lorentz, "i_g": i_g}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Peak profile: \n {:}\n {:}""".format(self["peak_shape"],  
                                                   self["asymmetry"])
        return lsout
    
    def calc_model(self, tth, tth_hkl, peak_shape = None, asymmetry = None, 
                   factor_lorentz = None):
        if peak_shape != None:
            self["peak_shape"] = peak_shape 
        if asymmetry != None:
            self["asymmetry"] = asymmetry
        if factor_lorentz != None:
            self["factor_lorentz"] = factor_lorentz
            
        peak_shape, asymmetry = self["peak_shape"], self["asymmetry"]
        factor_lorentz = self["factor_lorentz"]
        peak_shape["resolution"].calc_model(i_g = )
        factor_lorentz.calc_model()

        
    def pvoight_pd(self, tth, tth_hkl, resolution = None, assymetry = None):
        """
        pseudo voight function
        """
        if resolution != None :
            self["resolution"] = resolution
        if assymetry != None :
            self["assymetry"] = assymetry
            
        d_param = self["resolution"].calc_param(tth_hkl)
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl)
        self._p_ag = d_param["ag"]
        self._p_bg = d_param["bg"]
        self._p_al = d_param["al"]
        self._p_bl = d_param["bl"]
        eta_2d = d_param["eta"]
        self._p_eta = eta_2d 


        assym_2d = self["assymetry"].calc_assym(tth, tth_hkl)
        
        g_pd_2d = self.gauss_pd(tth_2d-tth_hkl_2d)
        l_pd_2d = self.lor_pd(tth_2d-tth_hkl_2d)
        res_2d = eta_2d * l_pd_2d + (1.-eta_2d) * g_pd_2d
        return res_2d*assym_2d
            

    

cc = Variable(2,True,"sqdf sfd")
bb = Resolution(W=2,U=cc,V=3)    
cc[0]=74
bb["U"]
bb["U"]
print (bb.calc_hg(41))
    
if (__name__ == "__main__"):
  pass
