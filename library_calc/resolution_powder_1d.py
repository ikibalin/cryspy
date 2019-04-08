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
    def __init__(self, u = 0, v = 0, w = 0.01, x = 0, y = 0):
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

        dd= {"u":u, "v":v, "w":w, "x":x, "y":y}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Resolution: 
 U {:}\n V {:}\n W {:}\n X {:}\n Y {:}""".format(self["u"],  
 self["v"], self["w"], self["x"], self["y"])
        return lsout
    
    def _calc_tancos(self, tth_hkl):
        self._p_t_tth = numpy.tan(tth_hkl)
        self._p_t_tth_sq = self._p_t_tth**2
        res = numpy.cos(tth_hkl)
        self._p_c_tth = res
        self._p_ic_tth = 1./res
        
    def _calc_hg(self, i_g = 0.):
        """
        ttheta in radians, could be array
        gauss size
        """
        res_sq = (self["u"]*self._p_t_tth_sq + self["v"]*self._p_t_tth + 
                  self["w"]*self._p_c_tth + i_g*self._p_ic_tth**2)
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
    
    def calc_model(self, tth_hkl, i_g = 0., u = None, v = None, w = None, 
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
        if x != None:
            self["x"] = x
        if y != None:
            self["y"] = y
            
        self._calc_tancos(tth_hkl)
        self._calc_hg(i_g = i_g)
        self._calc_hl()
        self._calc_hpveta()
        self._calc_agbg()
        self._calc_albl()
        
        d_out = dict(a_g = self._p_ag, b_g = self._p_bg, 
                     a_l = self._p_al, b_l = self._p_bl,
                     h_g = self._p_hg, h_l = self._p_hl,
                     h_pv = self._p_hpv, eta = self._p_eta)
        self.update(d_out)
        return 

    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["u", "v", "w", "x", "y"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]

        if refresh:
            tth = self["tth"]
            tth_hkl = self["tth_hkl"]
            i_g = self["i_g"]
            self.calc_model(tth, tth_hkl, i_g=i_g)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = ResolutionPD(u = self["u"], v = self["v"], w = self["w"], 
                               x = self["x"], y = self["y"])
        llab = ["tth_hkl", "i_g"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new
        



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

    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        #keys = d_vals.keys()
        #llab = []
        #llab_in = [(hh in keys) for hh in llab]
        #for lab, cond_h in zip(llab, llab_in):
        #    if cond_h:
        #        self[lab] = d_vals[lab]

        if refresh:
            tth = self["tth"]
            self.calc_model(tth)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = FactorLorentzPD()
        llab = ["tth"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new


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

        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl, indexing="ij")
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

        res_2d = np_one+val_1+val_2
        d_out = dict(asymmetry = res_2d, tth = tth, tth_hkl = tth_hkl)
        self.update(d_out)
        return
    
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]


        if refresh:
            tth = self["tth"]
            tth_hkl = self["tth_hkl"]
            self.calc_model(tth, tth_hkl)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = AsymmetryPD(p1 = self["p1"], p2 = self["p2"], 
                              p3 = self["p3"], p4 = self["p4"])
        llab = ["tth", "tth_hkl"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new



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
        self._p_eta = None
        self._p_gauss_pd = None
        self._p_lor_pd = None
        self.update(dd)
        
    def __repr__(self):
        lsout = """Shape of the profile: \n {:}""".format(self["resolution"])
        return lsout
        
    def _gauss_pd(self, tth_2d):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self._p_ag, self._p_bg
        self._p_gauss_pd = ag*numpy.exp(-bg*tth_2d**2)
        
    def _lor_pd(self, tth_2d):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self._p_al, self._p_bl
        self._p_lor_pd = al*1./(1.+bl*tth_2d**2)
    
    def calc_model(self, tth, tth_hkl, i_g = 0., resolution = None):
        """
        pseudo voight function
        """
        if resolution != None :
            self["resolution"] = resolution
            
        resolution = self["resolution"]            
        resolution.calc_model(tth_hkl, i_g = i_g)
        
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl, indexing="ij")
        self._p_ag = numpy.meshgrid(tth, resolution["a_g"], indexing="ij")[1]
        self._p_bg = numpy.meshgrid(tth, resolution["b_g"], indexing="ij")[1]
        self._p_al = numpy.meshgrid(tth, resolution["a_l"], indexing="ij")[1]
        self._p_bl = numpy.meshgrid(tth, resolution["b_l"], indexing="ij")[1]
        eta_2d = numpy.meshgrid(tth, resolution["eta"], indexing="ij")[1]
        self._p_eta = eta_2d 

        self._gauss_pd(tth_2d-tth_hkl_2d)
        self._lor_pd(tth_2d-tth_hkl_2d)
        g_pd_2d = self._p_gauss_pd 
        l_pd_2d = self._p_lor_pd
        res_2d = eta_2d * l_pd_2d + (1.-eta_2d) * g_pd_2d
        
        
        d_out = dict(profile = res_2d)
        self.update(d_out)
        

    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["resolution", "asymmetry", "factor_lorentz"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]


        llab = ["u", "v", "w", "x", "y"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["resolution"].set_vals(d_val_as, refresh)

        if refresh:
            tth = self["tth"]
            tth_hkl = self["tth_hkl"]
            i_g = self["i_g"]
            self.calc_model(tth, tth_hkl, i_g=i_g)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = PeakShapePD(resolution = self["resolution"])
        llab = ["tth", "tth_hkl", "i_g", "profile"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new
        

class PeakProfilePD(dict):
    """
    Peak Profile
    """
    def __init__(self, tth_hkl = 0.0, i_g = 0, peak_shape = PeakShapePD(), 
                asymmetry = AsymmetryPD(), factor_lorentz = FactorLorentzPD()):
        super(PeakProfilePD, self).__init__()
        dd= {"tth_hkl": tth_hkl, "i_g": i_g, "peak_shape": peak_shape, 
             "asymmetry": asymmetry, "factor_lorentz": factor_lorentz}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Peak profile: \n tth_hkl: {:}\n i_g {:}\n {:}\n {:}""".format(
                self["tth_hkl"],self["i_g"],self["peak_shape"], self["asymmetry"])
        return lsout
    
    def calc_model(self, tth, tth_hkl = None, i_g = None, peak_shape = None, 
                   asymmetry = None, factor_lorentz = None):
        if tth_hkl != None:
            self["tth_hkl"] = tth_hkl
        if i_g != None:
            self["i_g"] = i_g
        if peak_shape != None:
            self["peak_shape"] = peak_shape 
        if peak_shape != None:
            self["peak_shape"] = peak_shape 
        if asymmetry != None:
            self["asymmetry"] = asymmetry
        if factor_lorentz != None:
            self["factor_lorentz"] = factor_lorentz
            
        i_g, tth_hkl = self["i_g"], self["tth_hkl"]
        peak_shape, asymmetry = self["peak_shape"], self["asymmetry"]
        factor_lorentz = self["factor_lorentz"]
        peak_shape.calc_model(tth, tth_hkl, i_g = i_g)
        asymmetry.calc_model(tth, tth_hkl)
        factor_lorentz.calc_model(tth)
        
        np_shape_2d = peak_shape["profile"]
        np_ass_2d = asymmetry["asymmetry"]
        np_lor_1d = factor_lorentz["factor_lorentz"]
        np_lor_2d = numpy.meshgrid(np_lor_1d, tth_hkl, indexing="ij")[0]
        
        np_res_2d = np_shape_2d*np_ass_2d*np_lor_2d
        
        d_out = dict(profile = np_res_2d, tth = tth)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["tth_hkl", "i_g", "peak_shape", "asymmetry", "factor_lorentz"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]

        llab = ["resolution", "u", "v", "w", "x", "y"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["peak_shape"].set_vals(d_val_as, refresh)

        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["asymmetry"].set_vals(d_val_as, refresh)
            
        if refresh:
            tth = self["tth"]
            self.calc_model(tth)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = PeakProfilePD(tth_hkl = self["tth_hkl"], i_g = self["i_g"], 
                                peak_shape = self["peak_shape"], 
                                asymmetry= self["asymmetry"], 
                                factor_lorentz= self["factor_lorentz"])
        llab = ["profile", "tth"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new
             
#import Variable
#v_u = Variable.Variable(0.2)
#v_i_g_1=Variable.Variable(0.1)
#v_i_g_2=Variable.Variable(0.2)
#v_u.value = 0.777
#v_i_g_1.value = 0.222
#v_i_g_2.value = 0.333

#rr_1 = ResolutionPD(u = v_u, w = v_i_g_1)

#rr_2 = rr_1.soft_copy()

        
if (__name__ == "__main__"):
  pass
