"""
define resolution for the powder diffractometer along ttheta
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import numpy


class UnitCell(dict):
    """
    Unit Cell parameters
    """
    def __init__(self, a = 1.0, b = 1.0, c = 1.0, alpha = 1.0, beta = 90.0, 
                 gamma= 90., singony = "Monoclinic"):
        super(UnitCell, self).__init__()
        self._p_cos_a = None
        self._p_cos_b = None
        self._p_cos_g = None
        self._p_cos_a_sq = None
        self._p_cos_b_sq = None
        self._p_cos_g_sq = None
        self._p_sin_a = None
        self._p_sin_b = None
        self._p_sin_g = None
        self._p_sin_a_sq = None
        self._p_sin_b_sq = None
        self._p_sin_g_sq = None
        
        self._p_ia = None
        self._p_ib = None
        self._p_ic = None
        self._p_ialpha = None
        self._p_ibeta = None
        self._p_igamma = None        

        self._p_cos_ia = None
        self._p_cos_ib = None
        self._p_cos_ig = None
        self._p_cos_ia_sq = None
        self._p_cos_ib_sq = None
        self._p_cos_ig_sq = None
        self._p_sin_ia = None
        self._p_sin_ib = None
        self._p_sin_ig = None
        self._p_sin_ia_sq = None
        self._p_sin_ib_sq = None
        self._p_sin_ig_sq = None
        
        self._p_vol = None
        self._p_ivol = None
        self._p_m_b = None
        self._p_m_ib = None

        dd= {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma,
             "singony": singony}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Unit cell: \n a: {:}\n b: {:}\n c: {:}\n alpha: {:}
 beta: {:}\n gamma: {:}\n singony: {:}""".format(self["a"], self["b"], 
                 self["c"], self["alpha"], self["beta"], self["gamma"], self["singony"])
        return lsout
    
    def _constr_singony(self):
        singony = self["singony"]
        if singony == "Cubic":
            self["b"] = self["a"]
            self["c"] = self["a"]
            self["alpha"] = 90.
            self["beta"] = 90.
            self["gamma"] = 90.
        elif singony == "Hexagonal":
            self["b"] = self["a"]
            self["alpha"] = 90.
            self["beta"] = 90.
            self["gamma"] = 120.        
        elif singony == "Trigonal":
            self["b"] = self["a"]
            self["c"] = self["a"]
        elif singony == "Tetragonal":
            self["b"] = self["a"]
            self["alpha"] = 90.
            self["beta"] = 90.
            self["gamma"] = 90.
        elif singony == "Orthorhombic":
            self["alpha"] = 90.
            self["beta"] = 90.
            self["gamma"] = 90.
        elif singony == "Monoclinic":
            self["alpha"] = 90.
            self["gamma"] = 90.
            
    def _calc_cos_abc(self):
        rad=numpy.pi/180.
        self._p_cos_a = numpy.cos(self["alpha"]*rad)
        self._p_cos_b = numpy.cos(self["beta"]*rad)
        self._p_cos_g = numpy.cos(self["gamma"]*rad)
        
        self._p_sin_a = numpy.sin(self["alpha"]*rad)
        self._p_sin_b = numpy.sin(self["beta"]*rad)
        self._p_sin_g = numpy.sin(self["gamma"]*rad)
        
        self._p_cos_a_sq = self._p_cos_a**2
        self._p_cos_b_sq = self._p_cos_b**2
        self._p_cos_g_sq = self._p_cos_g**2

        self._p_sin_a_sq = 1.-self._p_cos_a_sq
        self._p_sin_b_sq = 1.-self._p_cos_b_sq
        self._p_sin_g_sq = 1.-self._p_cos_g_sq
        
    def _calc_cos_iabc(self):
        rad=numpy.pi/180.
        self._p_cos_ia = numpy.cos(self._p_ialpha*rad)
        self._p_cos_ib = numpy.cos(self._p_ibeta*rad)
        self._p_cos_ig = numpy.cos(self._p_igamma*rad)
        
        self._p_sin_ia = numpy.sin(self._p_ialpha*rad)
        self._p_sin_ib = numpy.sin(self._p_ibeta*rad)
        self._p_sin_ig = numpy.sin(self._p_igamma*rad)
        
        self._p_cos_ia_sq = self._p_cos_ia**2
        self._p_cos_ib_sq = self._p_cos_ib**2
        self._p_cos_ig_sq = self._p_cos_ig**2

        self._p_sin_a_sq = 1.-self._p_cos_a_sq
        self._p_sin_b_sq = 1.-self._p_cos_b_sq
        self._p_sin_g_sq = 1.-self._p_cos_g_sq

    def _calc_volume(self):
        a = self["a"]
        b = self["b"]
        c = self["c"]
        c_a = self._p_cos_a
        c_b = self._p_cos_b
        c_g = self._p_cos_g
        c_a_sq = self._p_cos_a_sq
        c_b_sq = self._p_cos_b_sq
        c_g_sq = self._p_cos_g_sq
        vol = a*b*c*(1.-c_a_sq-c_b_sq-c_g_sq+2.*c_a*c_b*c_g)**0.5
        self._p_vol = vol
        
    
    def _calc_iucp(self):
        """
        calculate inverse unit cell
        """
        irad = 180./numpy.pi

        a = self["a"]
        b = self["b"]
        c = self["c"]
        c_a = self._p_cos_a
        c_b = self._p_cos_b
        c_g = self._p_cos_g
        s_a = self._p_sin_a
        s_b = self._p_sin_b
        s_g = self._p_sin_g
        vol = self._p_vol
        
        self._p_ialpha = numpy.arccos((c_b*c_g-c_a)/(s_b*s_g))*irad
        self._p_ibeta = numpy.arccos((c_g*c_a-c_b)/(s_g*s_a))*irad
        self._p_igamma = numpy.arccos((c_a*c_b-c_g)/(s_a*s_b))*irad

        self._p_ia = b*c*s_a/vol
        self._p_ib = c*a*s_b/vol
        self._p_ic = a*b*s_g/vol


    def _calc_m_b(self):
        """
        calculate matrix B 
        """
        c = self["c"] 

        ia = self._p_ia 
        ib = self._p_ib 
        ic = self._p_ic 
        
        c_a = self._p_cos_a
        
        #ic_a = self._p_cos_ia 
        ic_b = self._p_cos_ib 
        ic_g = self._p_cos_ig 
        #is_a = self._p_sin_ia 
        is_b = self._p_sin_ib 
        is_g = self._p_sin_ig 
        
        self._p_m_b = numpy.array([[ia,  ib*ic_g,  ic*ic_b],
            [0.,  ib*is_g, -ic*is_b*c_a],
            [0.,       0.,  1./c]], dtype = float)

    def _calc_m_ib(self):
        """
        calculate inverse B matrix 
        """
        x1 = self._p_m_b[0,0]
        x2 = self._p_m_b[1,1]
        x3 = self._p_m_b[2,2]
        x4 = self._p_m_b[0,1]
        x5 = self._p_m_b[0,2]
        x6 = self._p_m_b[1,2]
        #B=[[x1,x4,x5],
        #   [0.,x2,x6],
        #   [0.,0.,x3]]
        #it shuld be checked
        #iB=numpy.linalg.inv(B)
        y1 = 1./x1
        y2 = 1./x2
        y3 = 1./x3
        y4 = -1*x4*1./(x1*x2)
        y6 = -1*x6*1./(x2*x3)
        y5 = (x4*x6-x2*x5)*1./(x1*x2*x3)
        
        self._p_m_ib = numpy.array([[y1,y4,y5],[0.,y2,y6],[0.,0.,y3]], 
                                   dtype = float)
            
                
        
    def _refresh(self, a, b, c, alpha, beta, gamma, singony):
        """
        refresh variables
        """
        if a != None:
            self["a"] = a
        if b != None:
            self["b"] = b
        if c != None:
            self["c"] = c
        if alpha != None:
            self["alpha"] = alpha
        if beta != None:
            self["beta"] = beta
        if gamma != None:
            self["gamma"] = gamma
        if singony != None:
            self["singony"] = singony
        self._constr_singony()
    
    def calc_sthovl(self, h, k, l, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None):
        """
        calculate sin(theta)/lambda for list of hkl reflections
        """
        cond = any([hh != None for hh in [a, b, c, alpha, beta, gamma, singony]])
        if cond:
            self.calc_model(a, b, c, alpha, beta, gamma, singony)
        a = self["a"]
        b = self["b"]
        c = self["c"]
        c_a = self._p_cos_a
        c_b = self._p_cos_b
        c_g = self._p_cos_g
        c_a_sq = self._p_cos_a_sq
        c_b_sq = self._p_cos_b_sq
        c_g_sq = self._p_cos_g_sq
        s_a_sq = self._p_sin_a_sq
        s_b_sq = self._p_sin_b_sq
        s_g_sq = self._p_sin_g_sq

        A=( 1. - c_a_sq - c_b_sq - c_g_sq + 2.*c_a*c_b*c_g)
        B1 = (s_a_sq*(h*1./a)**2+s_b_sq*(k*1./b)**2+s_g_sq*(l*1./c)**2)
        B2 = 2.*(k*l*c_a)/(b*c)+2.*(h*l*c_b)/(a*c)+2.*(h*k*c_g)/(a*b)
        #it should be checked, I am not sure
        B = B1-B2
        inv_d = (B*1./A)**0.5
        return 0.5*inv_d
        
    def calc_model(self, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None):
        self._refresh(a, b, c, alpha, beta, gamma, singony)
        self._calc_cos_abc()
        self._calc_volume()
        self._calc_iucp()
        self._calc_cos_iabc()
        self._calc_m_b()
        self._calc_m_ib()
        volume = self._p_vol
        ia = self._p_ia 
        ib = self._p_ib 
        ic = self._p_ic 
        ialpha = self._p_ialpha 
        ibeta = self._p_ibeta
        igamma = self._p_igamma 
        matrix_B = self._p_m_b
        matrix_iB = self._p_m_ib
        d_out = dict(volume=volume, matrix_B=matrix_B, matrix_iB=matrix_iB,
                     ia=ia, ib=ib, ic=ic,
                     ialpha=ialpha, ibeta=ibeta, igamma=igamma)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["a", "b", "c", "alpha", "beta", "gamma", "singony"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
        self._constr_singony()
                
            
        if refresh:
            #tth = self["tth"]
            self.calc_model()

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = UnitCell(a=self["a"], b= self["b"], c = self["c"], 
                           alpha = self["alpha"], beta=self["beta"], 
                           gamma=self["gamma"], singony=self["singony"])
        
        llab = ["volume", "matrix_B", "ia", "ib", "ic", "ialpha", "ibeta", 
                "igamma"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new


class CrystSymmetry(dict):
    """
    Crystall symmetry
    """
    def __init__(self, spgr_name = "P1", spgr_number = 1, spgr_choise = "1"):
        super(CrystSymmetry, self).__init__()
        dd= {"spgr_name": spgr_name, "spgr_number": spgr_number, 
             "spgr_choise": spgr_choise}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Space group: \n name: {:}\n number: {:}\n choise: {:}
""".format(self["spgr_name"], self["spgr_number"],self["spgr_choise"])
        return lsout
    
    def calc_model(self, spgr_name = None, spgr_number = None, 
                   spgr_choise = None):
        if spgr_name != None:
            self["spgr_name"] = spgr_name
        if spgr_number != None:
            self["spgr_number"] = spgr_number
        if spgr_choise != None:
            self["spgr_choise"] = spgr_choise

        d_out = dict()
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["spgr_name", "spgr_number", "spgr_choise"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
        self._constr_singony()
                
            
        if refresh:
            #tth = self["tth"]
            self.calc_model()

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = UnitCell(spgr_name=self["spgr_name"], 
                           spgr_number= self["spgr_number"], 
                           spgr_choise = self["spgr_choise"])
        #llab = ["fn", "sft"]
        #keys = self.keys()
        #llab_in = [(hh in keys) for hh in llab]
        #for lab, cond_h in zip(llab, llab_in):
        #    if cond_h:
        #        obj_new[lab] = self[lab]
        return obj_new

class DebyeWaller(dict):
    """
    DebyeWaller
    """
    def __init__(self, beta_11 = 0., beta_22 = 0., beta_33 = 0., 
                 beta_12 = 0., beta_13 = 0., beta_23 = 0., b_iso = 0.):
        super(DebyeWaller, self).__init__()
        dd= {"beta_11": beta_11, "beta_22": beta_22, "beta_33": beta_33,
             "beta_12": beta_12, "beta_13": beta_13, "beta_23": beta_23,
             "b_iso": b_iso}
        self.update(dd)
        

class Magnetism(dict):
    """
    Magnetism
    """
    def __init__(self, kappa = 1.0, factor_lande = 2.0, type = "",
                 chi_11 = 0., chi_22 = 0., chi_33 = 0., 
                 chi_12 = 0., chi_13 = 0., chi_23 = 0.):
        super(Magnetism, self).__init__()
        dd= {"kappa": kappa, "factor_lande ": factor_lande , 
             "chi_11": chi_11, "chi_22": chi_22, "chi_33": chi_33,
             "chi_12": chi_12, "chi_13": chi_13, "chi_23": chi_23}
        self.update(dd)

    
class Atom(dict):
    """
    Atom
    """
    def __init__(self, magnetism = Magnetism(), 
                 debye_waller = DebyeWaller(), position = Position(),
                 type = ""):
        super(Atom, self).__init__()
        dd= {"crystal_symmetry": crystal_symmetry, "unit_cell": unit_cell, 
             "atom": atom}
        self.update(dd)

    
class Crystal(dict):
    """
    Crystal
    """
    def __init__(self, crystal_symmetry = CrystSymmetry(), 
                 unit_cell = UnitCell(), atom = Atom()):
        super(Crystal, self).__init__()
        dd= {"crystal_symmetry": crystal_symmetry, "unit_cell": unit_cell, 
             "atom": atom}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Phase: \n crystal symmetry: {:}\n unit cell {:}
 atom {:}""".format(self["crystal_symmetry"], self["unit_cell"], self["atom"])
        return lsout
    
    def calc_model(self, reflection, crystal_symmetry = None, unit_cell = None,
                   atom = None):
        if crystal_symmetry != None:
            self["crystal_symmetry"] = crystal_symmetry 
        if unit_cell != None:
            self["unit_cell"] = unit_cell 
        if atom != None:
            self["atom"] = atom 
            
        np_fn_1d = None #complex 1D numpy array of nuclear structure factors over list of hkl
        np_sft_2d = None #complex 2D numpy array of tensor structure factors /11, 22, 33, 12, 13, 23/... over list of hkl 
        d_out = dict(fn = np_fn_1d, sft = np_sft_2d)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["crystal_symmetry", "unit_cell", "atom"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
                
        #redo it
        llab = ["resolution", "u", "v", "w", "x", "y"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["crystal_symmetry"].set_vals(d_val_as, refresh)

        #redo it
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["unit_cell"].set_vals(d_val_as, refresh)

        #redo it
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["atom"].set_vals(d_val_as, refresh)
            
        if refresh:
            tth = self["tth"]
            self.calc_model(tth)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Crystal(crystal_symmetry = self["crystal_symmetry"], 
                          unit_cell = self["unit_cell"], atom = self["atom"])
        llab = ["fn", "sft"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new


             
#import Variable
#v_u = Variable.Variable(0.2)

        
if (__name__ == "__main__"):
  pass
