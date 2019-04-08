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
        d_out = dict(sthovl = 0.5*inv_d)
        self.update(d_out)
        
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
    def __init__(self, spgr_given_name = "P1", spgr_choiсe = "1",
                 f_dir_prog = os.getcwd()):
        super(CrystSymmetry, self).__init__()
        
        if isinstance(spgr_choiсe, float):
            spgr_choise = "{:}".format(int(spgr_choiсe))



        f_itables=os.path.join(f_dir_prog,"itables.txt")
        self._read_el_cards(f_itables)
        
        dd= {"spgr_given_name": spgr_given_name, "spgr_choiсe": spgr_choiсe, 
             "f_dir_prog": f_dir_prog}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Space group: \n name: {:}\n choiсe: {:}
 directory: '{:}'""".format(self["spgr_name"], self["spgr_choiсe"], 
                            self["f_dir_prog"])
        return lsout


    def _read_el_cards(self, f_itables):
        """
        reading information about space grooupe from file fitables to list of cards ldcard
        Info in file fitables:
        
        1 P1               Triclinic
        choice: 1
        centr: false
        pcentr: 0, 0, 0
        symmetry: X,Y,Z
        
        2 P-1              Triclinic
        ...
        """
        fid = open(f_itables, "r")
        lcontent = fid.readlines()
        fid.close()
    
        lcontent = [hh.strip() for hh in lcontent if hh.strip() != ""]
        ldcard = []
        dcard = None
        for hh in lcontent:
            lhelp = hh.split()
            if lhelp[0].isdigit():
                if dcard != None:
                    ldcard.append(dcard)
                dcard = {"number":lhelp[0], "name": lhelp[1], "syngony": lhelp[2]}
            else:
                lhelp = hh.split(":")
                if (lhelp[0].strip() in dcard.keys()):
                    dcard[lhelp[0].strip()].append(lhelp[1].strip())
                else:
                    dcard[lhelp[0].strip()] = [lhelp[1].strip()]
        d_out = dict(spgr_table = ldcard)
        self.update(d_out)

    def _get_symm(self):
        """
        get symmetry from space group
        """
        
        spgr_choiсe = self["spgr_choiсe"]
        
        spgr_given_name = self["spgr_given_name"]
        

        if spgr_given_name.isdigit():
            spgr_n = spgr_given_name
            spgr_name = ""
        else:
            spgr_n = ""
            spgr_name = spgr_given_name
        
        spgr_table = self["spgr_table"]
        
        for dcard in spgr_table:
            if (((dcard["number"] == spgr_n)|(dcard["name"] == spgr_name))&(dcard["choiсe"][0] == spgr_choiсe)):
                flag = True
                break
        if (not flag):
            print("Space groupe is not found")
            return
        
        flag = False
            
        lelsymm = []
        for ssymm in dcard["symmetry"]:
            lelsymm.append(cfunc.tr2elsymm(ssymm))
        centr = dcard["centr"][0]=="true"
        pcentr = [float(hh) for hh in dcard["pcentr"][0].split(",")]
        fletter = dcard["name"][0]
        spgr = dcard["name"]
        number = int(dcard["number"])
        if (fletter == "P"):
            lorig = [(0, 0, 0)]
        elif fletter == "C":
            lorig = [(0, 0, 0), (0.5, 0.5, 0)]
        elif fletter == "I":
            lorig = [(0, 0, 0), (0.5, 0.5, 0.5)]
        elif fletter == "F":
            lorig = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        elif (fletter == "R"):
            if n_choise == "1":
                lorig = [(0, 0, 0), (0.66667, 0.33333, 0.33333), (0.33334, 0.66666, 0.66666)]
            else:
                lorig = [(0, 0, 0)]
        else:
            print("Undefined syngony")

            
        d_symm = {"centr":centr, "el_symm":lelsymm, "orig":lorig, 
                  "p_centr":pcentr, "spgr_name":spgr, "spgr_number":number}

        self.update(d_symm)


    def calc_model(self, spgr_given_name = None, spgr_choiсe = None,
                   f_dir_prog = None):
        
        if spgr_given_name != None:
            self["spgr_given_name"] = spgr_given_name 
        if spgr_choiсe != None:
            self["spgr_choiсe"] = spgr_choiсe
        if f_dir_prog != None:
            self["f_dir_prog"] = f_dir_prog
            f_itables=os.path.join(f_dir_prog,"itables.txt")
            self._read_el_cards(f_itables)            

        self._get_symm()
        d_out = dict()
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["spgr_given_name", "spgr_choiсe", "f_dir_prog"]
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
        obj_new = CrystSymmetry(spgr_given_name=self["spgr_given_name"], 
                                spgr_choiсe= self["spgr_choiсe"], 
                                f_dir_prog = self["f_dir_prog"])
        llab = ["fn", "sft"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new

    
    def _trans_str_to_el_symm(str1):
        """
        transform string to element of symmetry: (x,y,-z) -> 0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
        """
        str2="".join(str1.split(" "))
        lhelp1,lhelp2,lhelp3=[],[],[]
        lhelp1=[hh for hh in str2.split('(') if hh!=""]
        [lhelp2.extend(hh.split(')')) for hh in lhelp1 if hh!=""]
        [lhelp3.extend(hh.split(',')) for hh in lhelp2 if hh!=""]
        lAx=['x','y','z']
        lelsymm=[]
        for hh in lhelp3:
            elsymmh=[0.0,0,0,0]
            strh=hh
            for inum,Ax in enumerate(lAx):
                if (strh.find(Ax)!=-1):
                    if (strh.find("+"+Ax)!=-1):
                        elsymmh[inum+1]=1
                        strh="".join(strh.split("+"+Ax))
                    elif (strh.find("-"+Ax)!=-1):
                        elsymmh[inum+1]=-1
                        strh="".join(strh.split("-"+Ax))
                    else:
                        elsymmh[inum+1]=1
                        strh="".join(strh.split(Ax))
            if (strh==""):
                pass
            elif (strh.find("/")!=-1):
                lhelp1=strh.split("/")
                elsymmh[0]=float(lhelp1[0])/float(lhelp1[1])
            else:
                elsymmh[0]=float(strh)
            lelsymm.append(elsymmh)
        elsymm=[]
        [elsymm.extend(hh) for hh in lelsymm]
        return elsymm
    
  
class Position(dict):
    """
    Position of atom(s) in unit cell.
    """
    def __init__(self, x = 0., y = 0., z = 0.):
        super(Position, self).__init__()
        dd= {"x": x, "y": y, "z": z}
        self.update(dd)


    def __repr__(self):
        lsout = """Position: \n xyz: {:} {:} {:}""".format(
                self["x"], self["y"], self["z"])
        return lsout


    def _calc_phase(self, h, k, l, r_11, r_22, r_33, r_12, r_13, r_23):
        """
        calculate phase: exp(-2 pi i * (h*x+k*y+l*z))
        """
        x, y, z = self["x"], self["y"], self["z"]
        
        #redo: there should be also element of symmetry, and output is 3D instead of 2D
        np_x, np_h = numpy.meshgrid(x, h, indexing="ij")
        np_y, np_k = numpy.meshgrid(y, k, indexing="ij")
        np_z, np_l = numpy.meshgrid(z, l, indexing="ij")
        
        phase = numpy.exp(2*numpy.pi*1j*(np_x*np_h + np_y*np_k + np_z*np_l))
        
        d_out = dict(phase = phase)
        self.update(d_out)
        

    def calc_model(self, h, k, l, x = None, y = None, z = None):
        
        if x != None:
            self["x"] = x
        if y != None:
            self["y"] = y
        if z != None:
            self["z"] = z
        
        self._calc_phase(self, h, k, l)
        
        d_out = dict(h=h, k=k, l=l)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["x", "y", "z"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
            
        if refresh:
            h = self["h"] 
            k = self["k"] 
            l = self["l"] 
            
            self.calc_model(h, k, l)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Position(x=self["x"], y=self["y"], z=self["z"])
        llab = ["h", "k", "l", "phase"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
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

    def __repr__(self):
        lsout = """Debye Waller: \n beta: {:}\n choiсe: {:}
 directory: '{:}'""".format(self["spgr_name"], self["spgr_choiсe"], 
                            self["f_dir_prog"])
        return lsout


    def _calc_dwf_iso(self, sthovl):
        """
        isotropic harmonic Debye-Waller factor
        """
        b_iso = self["b_iso"]
        sthovl_sq = sthovl**2
        np_biso, np_sthovl_sq = numpy.meshgrid(b_iso, sthovl_sq, indexing="ij")
        dwf_iso = numpy.exp(-np_biso*np_sthovl_sq)
        d_out = dict(dwf_iso = dwf_iso)
        self.update(d_out)


    def _calc_dwf_aniso(self, h, k, l):
        """
        anisotropic harmonic Debye-Waller factor
        """
        beta_11, beta_22 = self["beta_11"], self["beta_22"] 
        beta_33, beta_12 = self["beta_33"], self["beta_12"]
        beta_13, beta_23 = self["beta_13"], self["beta_23"]
        
        h_sq, k_sq, l_sq = h**2, k**2, l**2
        hk, hl, kl = h*k, h*l, k*l
        
        dwf_aniso = numpy.exp(-1*(np_beta_11*np_h_sq + np_beta_22*np_k_sq + 
                                  np_beta_33*np_l_sq + 2.*np_beta_12*np_hk + 
                                  2.*np_beta_13*np_hl + 2.*np_beta_23*np_kl))
        
        d_out = dict(dwf_aniso = dwf_aniso)
        self.update(d_out)

        
    def calc_model(self, h, k, l, sthovl, beta_11 = None, beta_22 = None, 
                   beta_33 = None, beta_12 = None, beta_13 = None, beta_23 = None, b_iso = None):
        
        if beta_11 != None:
            self["beta_11"] = beta_11 
        if beta_22 != None:
            self["beta_22"] = beta_22 
        if beta_33 != None:
            self["beta_33"] = beta_33
        if beta_12 != None:
            self["beta_12"] = beta_12 
        if beta_13 != None:
            self["beta_13"] = beta_13 
        if beta_23 != None:
            self["beta_23"] = beta_23 
        
        self._calc_dwf_iso(self, sthovl)
        self._calc_dwf_aniso(self, h, k, l)
        
        d_out = dict(h=h, k=k, l=l, sthovl=sthovl)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["beta_11", "beta_22", "beta_33", "beta_12", "beta_13", 
                "beta_23"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
            
        if refresh:
            h = self["h"] 
            k = self["k"] 
            l = self["l"] 
            sthovl = self["sthovl"]
            
            self.calc_model(h, k, l, sthovl)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = DebyeWaller(beta_11 = self["beta_11"], 
                              beta_22 = self["beta_22"], 
                   beta_33 = self["beta_33"], beta_12 = self["beta_12"], 
                   beta_13 = self["beta_13"], beta_23 = self["beta_23"], 
                   b_iso = self["b_iso"])
        llab = ["dwf_iso", "dwf_aniso", "h", "k", "l", "sthovl"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new



class Magnetism(dict):
    """
    Magnetism
    """
    def __init__(self, kappa = 1.0, factor_lande = 2.0, type = "",
                 chi_11 = 0., chi_22 = 0., chi_33 = 0., 
                 chi_12 = 0., chi_13 = 0., chi_23 = 0.):
        super(Magnetism, self).__init__()
        dd= {"kappa": kappa, "factor_lande": factor_lande , 
             "chi_11": chi_11, "chi_22": chi_22, "chi_33": chi_33,
             "chi_12": chi_12, "chi_13": chi_13, "chi_23": chi_23}
        self.update(dd)
        
    
    def _calc_chi_loc(ia, ib, ic, matrix_ib):
        """
        representation of chi in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        expressions are taken from international tables
        matrix_ib is inversed matrix B
        ia, ib, ic is inversed unit cell parameters (it can be estimated from matrix matrix_ib)

        X = B x, x = iB X
        xT*CHI*x = XT iBT CHI iB X
    
        output chiLOC = iBT CHI iB
        """
        matrix_chi = numpy.array(
                [[self["chi_11"], self["chi_12"], self["chi_13"]],
                 [self["chi_12"], self["chi_22"], self["chi_23"]],
                 [self["chi_13"], self["chi_23"], self["chi_33"]]], 
                 dtype = float)
        #mchi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
        #[a,b,c,alpha,beta,gamma]=ucp
        y1 = matrix_ib[0,0]
        y2 = matrix_ib[1,1]
        y3 = matrix_ib[2,2]
        y4 = matrix_ib[0,1]
        y5 = matrix_ib[0,2]
        y6 = matrix_ib[1,2]
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
        matrix_ib_norm = matrix_ib
        matrix_ib_norm[:,0] *= ia
        matrix_ib_norm[:,1] *= ib
        matrix_ib_norm[:,2] *= ic
        
        matrix_ibt_norm = matrix_ib_norm.transpose()
        #it is not compatible with case, vhen chi_ij is 1D array 
        ibt_chi = numpy.matmul(matrix_ibt_norm, matrix_chi)
        matrix_chi_loc = numpy.matmul(ibt_chi, matrix_ib_norm)
        d_out = dict(matrix_chi_loc = matrix_chi_loc)
        self.update(d_out)
    
    def calc_model(self, ia, ib, ic, matrix_ib, chi_11=None, chi_22=None, 
                   chi_33=None, chi_12=None, chi_13=None, chi_23=None):
        
        if chi_11 != None:
            self["chi_11"] = chi_11 
        if chi_22 != None:
            self["chi_22"] = chi_22 
        if chi_33 != None:
            self["chi_33"] = chi_33
        if chi_12 != None:
            self["chi_12"] = chi_12 
        if chi_13 != None:
            self["chi_13"] = chi_13 
        if chi_23 != None:
            self["chi_23"] = chi_23 
        
        self._calc_chi_loc(ia, ib, ic, matrix_ib)
        d_out = dict(ia=ia, ib=ib, ic=ic, matrix_ib=matrix_ib)
        self.update(d_out)

    
    def calc_chi_rot(matrix_chi, elsymm):
        """
        calculate R*chi*RT
        rotation of chi by element of symmetry
        """
        [b1,r11,r12,r13,b2,r21,r22,r23,b3,r31,r32,r33]=elsymm
        matrix_r = numpy.array([[r11, r12, r13], [r21, r22, r23], 
                               [r31, r32, r33]], dtype=float)
        matrix_rt = matrix_r.transpose()
        r_chi = numpy.matmul(matrix_r, matrix_chi)
        
        matrix_chi_rot = numpy.matmul(r_chi, matrix_rt)
        return matrix_chi_rot 

    
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["chi_11", "chi_22", "chi_33", "chi_12", "chi_13", "chi_23"]
                
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
        self._constr_singony()
            
        if refresh:
            ia = self["ia"] 
            ib = self["ib"] 
            ic = self["ic"] 
            matrix_ib = self["matrix_ib"]
            self.calc_model(ia, ib, ic, matrix_ib)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Magnetism(chi_11=self["chi_11"], chi_22=self["chi_22"], 
                            chi_33=self["chi_33"], chi_12=self["chi_12"], 
                            chi_13=self["chi_13"], chi_23=self["chi_23"])
        llab = ["matrix_chi_loc", "ia", "ib", "ic", "matrix_ib"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new



    
    
class Atom(dict):
    """
    Atom
    """
    def __init__(self, atom_name="", atom_nucl_type="", b_scat=0.,
                 position=Position(), debye_waller=DebyeWaller(),
                 magnetism=Magnetism()):
        super(Atom, self).__init__()
        dd= {"atom_name": atom_name, "atom_nucl_type": atom_nucl_type, 
             "b_scat": b_scat, "position": position, 
             "debye_waller": debye_waller, "magnetism": magnetism,  }
        self.update(dd)
        
    def __repr__(self):
        lsout = """Atom: \n name: {:}\n nucl_type: {:}, b_scat: {}\n position {:}
 debye_waller: {}\n magnetism: {}""".format(self["atom_name"], 
 self["atom_nucl_type"],  self["b_scat"], self["position"], 
 self["debye_waller"], self["magnetism"])
        return lsout
    
    
    
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

    
    def calc_fn(self, h, k, l, crystal_symmetry = None, unit_cell = None, 
                atom = None):
        """
        calculate nuclear structure factor
        """
        if crystal_symmetry != None:
            self["crystal_symmetry"] = crystal_symmetry
        if unit_cell != None:
            self["unit_cell"] = unit_cell 
        if atom != None:
            self["atom"] = atom 

        crystal_symmetry = self["crystal_symmetry"]
        unit_cell = self["unit_cell"]
        atom = self["atom"] 
        
        unit_cell.calc_sthovl(h, k, l)
        sthovl = unit_cell["sthovl"]
        
        crystal_symmetry.calc_model()

        lelsymm = crystal_symmetry["elsymm"]
        lorig = crystal_symmetry["orig"]
        centr = crystal_symmetry["centr"]
        
        
        
        
        lFhelp=[0]
        lFhkl = []
        for hkl,sthovl in zip(lhkl,lsthovl):
            for cat in self.atom:
                occ,mult=cat.occ[0],cat.mult
                coord=[cat.coordx[0],cat.coordy[0],cat.coordz[0]]
                bscat=cat.bscat[0]
                AA=mult*occ*bscat
                if (AA!=0.):
                    for isymm,symm in enumerate(lelsymm):
                        hkls=[hkl[0]*symm[1]+hkl[1]*symm[5]+hkl[2]*symm[9], hkl[0]*symm[2]+hkl[1]*symm[6]+hkl[2]*symm[10], hkl[0]*symm[3]+hkl[1]*symm[7]+hkl[2]*symm[11]]
                        #if (cat.modedw=='iso'):
                        biso = cat.biso[0]
                        DWF_i = math.exp(-biso*sthovl**2)
                        #elif (cat.modedw=='aniso'):
                        #MISTAKE: introduce constraints on the betas
                        beta11,beta22=cat.beta11[0],cat.beta22[0]
                        beta33,beta12=cat.beta33[0],cat.beta12[0]
                        beta13,beta23=cat.beta13[0],cat.beta23[0]
                        #check it, not sure about 2 in front of beta12, beta13 and beta 23
                        #it is not correct, introduce equivalent reflections
                        DWF_a = math.exp(-1*(beta11*hkls[0]**2+beta22*hkls[1]**2+beta33*hkls[2]**2+
                                             2.*beta12*hkls[0]*hkls[1]+2.*beta13*hkls[0]*hkls[2]+2.*beta23*hkls[1]*hkls[2]))
                        #else:
                            #print 'Unknown thermal parameters for atom {0[name]:}'.format(datom)
                        #    print 'Program is stopped'
                        #    quit()
                        chelp = hkls[0]*coord[0]+hkl[0]*symm[0]+hkls[1]*coord[1]+hkl[1]*symm[4]+hkls[2]*coord[2]+hkl[2]*symm[8]
                        lFhelp.append(AA * DWF_i * DWF_a * cmath.exp(2*math.pi*1j*chelp))
            Fhklas=sum(lFhelp)*1./len(lelsymm)
            lFhelp=[]
            for orig in lorig:
                lFhelp.append(Fhklas*cmath.exp( 2*math.pi*1j* (hkl[0]*orig[0]+hkl[1]*orig[1]+hkl[2]*orig[2])))
            Fhklas=sum(lFhelp)*1./len(lFhelp)
            lFhelp=[]
            if (centr):
                orig = crystal_symmetry["p_centr"]
                Fhkl = (Fhklas+Fhklas.conjugate()*cmath.exp(2*2*math.pi*1j* (hkl[0]*orig[0]+hkl[1]*orig[1]+hkl[2]*orig[2])))*0.5
            else:
                Fhkl=Fhklas
            lFhkl.append(Fhkl)
            
            
            
            
        return lFhkl


    def calc_sft(self):
        """
        calculate structure factor tensor
        """
        pass
    
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
