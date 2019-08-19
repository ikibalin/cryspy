"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_atom_site import AbstractAtomSite


from neupy.f_crystal.cl_fract import Fract
from neupy.f_crystal.cl_magnetism import Magnetism
from neupy.f_crystal.cl_adp import ADP


class AtomSite(dict):
    """
    AtomSite
    """
    def __init__(self):
        super(AtomSite, self).__init__()
        self._p_fract = None
        self._p_adp = None
        self._p_magnetism = None
        self._p_b_scat = None
        self._p_occupation = None
        self._flag_refresh = False
        self._list_atom_type = []
    
    def __repr__(self):
        lsout = """AtomSite:\n number of atoms: {:}""".format(len(self._list_atom_type))
        for iatom, atom_type in enumerate(self._list_atom_type):
            lsout += "\n"+70*"*"+"\n {:}. ".format(iatom+1)+atom_type.__repr__()
        return lsout
    
    
    def _refresh(self):
        print("The option '_refresh' is not valiable")
        pass
    
    
    def set_val(self):
        print("The option 'set_val' is not valiable")
        pass
    
    
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

b_scat is amplitude scattering  
occupation is occupation factor  
fract  is fraction of atoms

        """
        print(lsout)
        
    def add_atom(self, atom):
        self._list_atom_type.append(atom)
        self._flag_refresh = True
        #self._form_arrays(cell)
    
    def del_atom(self, ind):
        self._list_atom_type.pop(ind)
        self._flag_refresh = True        
        #self._form_arrays(cell)

    def replace_atom(self, ind, atom):
        self._list_atom_type.pop(ind)
        self._list_atom_type.insert(ind, atom)
        self._flag_refresh = True
        #self._form_arrays(cell_)

    def _form_arrays(self, cell):
        lb_scat, locc = [], []
        lx, ly, lz = [], [], []
        lu_11, lu_22, lu_33 = [], [], []
        lu_12, lu_13, lu_23 = [], [], []
        lb_iso = []
        lchi_11, lchi_22, lchi_33 = [], [], []
        lchi_12, lchi_13, lchi_23 = [], [], []
        lkappa, lfactor_lande = [], []
        lj0_A, lj0_a, lj0_B, lj0_b, lj0_C, lj0_c = [], [], [], [], [], []
        lj0_D = []
        lj2_A, lj2_a, lj2_B, lj2_b, lj2_C, lj2_c = [], [], [], [], [], []
        lj2_D = []
        for atom_type in self._list_atom_type:
            lb_scat.append(1.*atom_type.get_val("b_scat"))
            locc.append(1.*atom_type.get_val("occupation"))
            lx.append(1.*atom_type.get_val("x"))
            ly.append(1.*atom_type.get_val("y"))
            lz.append(1.*atom_type.get_val("z"))
            adp_type = atom_type.get_val("adp_type")
            if adp_type == "uiso":
                lb_iso.append(1.*atom_type.get_val("b_iso"))
                lu_11.append(0.)
                lu_22.append(0.)
                lu_33.append(0.)
                lu_12.append(0.)
                lu_13.append(0.)
                lu_23.append(0.)
            else:
                lu_11.append(1.*atom_type.get_val("u_11"))
                lu_22.append(1.*atom_type.get_val("u_22"))
                lu_33.append(1.*atom_type.get_val("u_33"))
                lu_12.append(1.*atom_type.get_val("u_12"))
                lu_13.append(1.*atom_type.get_val("u_13"))
                lu_23.append(1.*atom_type.get_val("u_23"))
                lb_iso.append(0.)
            chi_type = atom_type.get_val("chi_type")
            if chi_type == "ciso":
                atom_type.apply_chi_iso(cell)
            lchi_11.append(1.*atom_type.get_val("chi_11"))
            lchi_22.append(1.*atom_type.get_val("chi_22"))
            lchi_33.append(1.*atom_type.get_val("chi_33"))
            lchi_12.append(1.*atom_type.get_val("chi_12"))
            lchi_13.append(1.*atom_type.get_val("chi_13"))
            lchi_23.append(1.*atom_type.get_val("chi_23"))
            lkappa.append(1.*atom_type.get_val("kappa"))
            lfactor_lande.append(1.*atom_type.get_val("factor_lande"))
            lj0_A.append(atom_type.get_val("j0_A"))
            lj0_a.append(atom_type.get_val("j0_a")) 
            lj0_B.append(atom_type.get_val("j0_B")) 
            lj0_b.append(atom_type.get_val("j0_b")) 
            lj0_C.append(atom_type.get_val("j0_C")) 
            lj0_c.append(atom_type.get_val("j0_c"))
            lj0_D.append(atom_type.get_val("j0_D"))
            lj2_A.append(atom_type.get_val("j2_A"))
            lj2_a.append(atom_type.get_val("j2_a")) 
            lj2_B.append(atom_type.get_val("j2_B")) 
            lj2_b.append(atom_type.get_val("j2_b")) 
            lj2_C.append(atom_type.get_val("j2_C")) 
            lj2_c.append(atom_type.get_val("j2_c"))
            lj2_D.append(atom_type.get_val("j2_D"))

        np_b_scat = numpy.array(lb_scat, dtype=complex)
        np_occ = numpy.array(locc, dtype=float)
        np_x = numpy.array(lx, dtype=float)
        np_y = numpy.array(ly, dtype=float)
        np_z = numpy.array(lz, dtype=float)
        fract = Fract(x=np_x, y=np_y, z=np_z)

        np_u_11 = numpy.array(lu_11, dtype=float)
        np_u_22 = numpy.array(lu_22, dtype=float)
        np_u_33 = numpy.array(lu_33, dtype=float)
        np_u_12 = numpy.array(lu_12, dtype=float)
        np_u_13 = numpy.array(lu_13, dtype=float)
        np_u_23 = numpy.array(lu_23, dtype=float)
        np_b_iso = numpy.array(lb_iso, dtype=float)
        adp = ADP(u_11=np_u_11, u_22=np_u_22, u_33=np_u_33,
                  u_12=np_u_12, u_13=np_u_13, u_23=np_u_23,
                  b_iso = np_b_iso)

        np_chi_11 = numpy.array(lchi_11, dtype=float)
        np_chi_22 = numpy.array(lchi_22, dtype=float)
        np_chi_33 = numpy.array(lchi_33, dtype=float)
        np_chi_12 = numpy.array(lchi_12, dtype=float)
        np_chi_13 = numpy.array(lchi_13, dtype=float)
        np_chi_23 = numpy.array(lchi_23, dtype=float)
        np_kappa = numpy.array(lkappa, dtype=float)
        np_factor_lande = numpy.array(lfactor_lande, dtype=float)
        np_j0_A = numpy.array(lj0_A, dtype=float)
        np_j0_a = numpy.array(lj0_a, dtype=float)
        np_j0_B = numpy.array(lj0_B, dtype=float)
        np_j0_b = numpy.array(lj0_b, dtype=float)
        np_j0_C = numpy.array(lj0_C, dtype=float)
        np_j0_c = numpy.array(lj0_c, dtype=float)
        np_j0_D = numpy.array(lj0_D, dtype=float)
        np_j2_A = numpy.array(lj2_A, dtype=float)
        np_j2_a = numpy.array(lj2_a, dtype=float)
        np_j2_B = numpy.array(lj2_B, dtype=float)
        np_j2_b = numpy.array(lj2_b, dtype=float)
        np_j2_C = numpy.array(lj2_C, dtype=float)
        np_j2_c = numpy.array(lj2_c, dtype=float)
        np_j2_D = numpy.array(lj2_D, dtype=float)
        

        magnetism = Magnetism(chi_11=np_chi_11, chi_22=np_chi_22, 
            chi_33=np_chi_33, chi_12=np_chi_12, chi_13=np_chi_13, 
            chi_23=np_chi_23, kappa = np_kappa, factor_lande = np_factor_lande,
            j0_A=np_j0_A, j0_a=np_j0_a, j0_B=np_j0_B, j0_b=np_j0_b, 
            j0_C=np_j0_C, j0_c=np_j0_c, j0_D=np_j0_D, j2_A=np_j2_A, 
            j2_a=np_j2_a, j2_B=np_j2_B, j2_b=np_j2_b, j2_C=np_j2_C, 
            j2_c=np_j2_c, j2_D=np_j2_D)
        
        self._p_b_scat = np_b_scat
        self._p_occupation = np_occ
        self._p_fract = fract
        self._p_adp = adp 
        self._p_magnetism = magnetism
    
    def calc_sf(self, space_group, cell, h, k, l, d_map={}):
        """
        calculate nuclear structure factor
        """
        #if d_map == {}:
        #    d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = d_map["out"]
        #    return f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33
        if (self.is_variable()|self._flag_refresh):
            self._form_arrays(cell)
        #sthovl = cell.calc_sthovl(h, k, l)

        fract = self._p_fract
        adp = self._p_adp
        magnetism = self._p_magnetism
        b_scat = self._p_b_scat
        occupation = self._p_occupation
        x, y, z = fract.get_val("x"), fract.get_val("y"), fract.get_val("z")
        atom_multiplicity = space_group.calc_atom_mult(x, y, z)
        occ_mult = occupation*atom_multiplicity 
        
        #d_phase = d_map["phase"]
        #if not(d_phase["flag"]|(d_phase["out"] is None)):
        #    phase_3d = d_phase["out"] 
        #else:
        phase_3d = fract.calc_phase(space_group, h, k, l)#3d object
        #d_phase["out"] = phase_3d 
        
        #d_adp = d_map["adp"]
        #if not(d_adp["flag"]|(d_adp["out"] is None)):
        #    dwf_3d = d_adp["out"] 
        #else:
        dwf_3d = adp.calc_dwf(space_group, cell, h, k, l)
        #    d_adp["out"] = dwf_3d 
        
        #d_magnetism = d_map["magnetism"]
        #if not(d_magnetism["flag"]|(d_magnetism["out"] is None)):
        #    ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33 = d_magnetism["out"] 
        #else:
        ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33 =  magnetism.calc_form_factor_tensor(space_group, cell, h, k, l)
        #    d_magnetism["out"] = (ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33)
        
        hh = phase_3d*dwf_3d
        
        phase_2d = hh.sum(axis=2)

        ft_11 = (ff_11*hh).sum(axis=2)
        ft_12 = (ff_12*hh).sum(axis=2)
        ft_13 = (ff_13*hh).sum(axis=2)
        ft_21 = (ff_21*hh).sum(axis=2)
        ft_22 = (ff_22*hh).sum(axis=2)
        ft_23 = (ff_23*hh).sum(axis=2)
        ft_31 = (ff_31*hh).sum(axis=2)
        ft_32 = (ff_32*hh).sum(axis=2)
        ft_33 = (ff_33*hh).sum(axis=2)
        
        b_scat_2d = numpy.meshgrid(h, b_scat, indexing="ij")[1]
        occ_mult_2d = numpy.meshgrid(h, occ_mult, indexing="ij")[1]
        
        lel_symm = space_group.get_val("el_symm")
        lorig = space_group.get_val("orig")
        centr = space_group.get_val("centr")

        #calculation of nuclear structure factor        
        hh = phase_2d * b_scat_2d * occ_mult_2d
        f_hkl_as = hh.sum(axis=1)*1./len(lel_symm)
        
        
        orig_x = [hh[0] for hh in lorig]
        orig_y = [hh[1] for hh in lorig]
        orig_z = [hh[2] for hh in lorig]
        
        np_h, np_orig_x = numpy.meshgrid(h, orig_x, indexing = "ij")
        np_k, np_orig_y = numpy.meshgrid(k, orig_y, indexing = "ij")
        np_l, np_orig_z = numpy.meshgrid(l, orig_z, indexing = "ij")
        
        np_orig_as = numpy.exp(2*numpy.pi*1j*(np_h*np_orig_x+np_k*np_orig_y+np_l*np_orig_z))
        f_hkl_as = f_hkl_as*np_orig_as.sum(axis=1)*1./len(lorig)

        if (centr):
            orig = space_group.get_val("p_centr")
            f_nucl = 0.5*(f_hkl_as+f_hkl_as.conjugate()*numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2])))
        else:
            f_nucl = f_hkl_as

        #calculation of structure factor tensor
        sft_as_11 = (ft_11 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_12 = (ft_12 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_13 = (ft_13 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_21 = (ft_21 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_22 = (ft_22 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_23 = (ft_23 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_31 = (ft_31 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_32 = (ft_32 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)
        sft_as_33 = (ft_33 * occ_mult_2d).sum(axis=1)*1./len(lel_symm)

        sft_as_11 = sft_as_11 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_12 = sft_as_12 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_13 = sft_as_13 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_21 = sft_as_21 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_22 = sft_as_22 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_23 = sft_as_23 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_31 = sft_as_31 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_32 = sft_as_32 * np_orig_as.sum(axis=1)*1./len(lorig)
        sft_as_33 = sft_as_33 * np_orig_as.sum(axis=1)*1./len(lorig)
    
        if (centr):
            orig = space_group.get_val("p_centr")
            hh = numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2]))
            sft_11 = 0.5*(sft_as_11+sft_as_11.conjugate()*hh)
            sft_12 = 0.5*(sft_as_12+sft_as_12.conjugate()*hh)
            sft_13 = 0.5*(sft_as_13+sft_as_13.conjugate()*hh)
            sft_21 = 0.5*(sft_as_21+sft_as_21.conjugate()*hh)
            sft_22 = 0.5*(sft_as_22+sft_as_22.conjugate()*hh)
            sft_23 = 0.5*(sft_as_23+sft_as_23.conjugate()*hh)
            sft_31 = 0.5*(sft_as_31+sft_as_31.conjugate()*hh)
            sft_32 = 0.5*(sft_as_32+sft_as_32.conjugate()*hh)
            sft_33 = 0.5*(sft_as_33+sft_as_33.conjugate()*hh)
        else:
            sft_11, sft_12, sft_13 = sft_as_11, sft_as_12, sft_as_13
            sft_21, sft_22, sft_23 = sft_as_21, sft_as_22, sft_as_23
            sft_31, sft_32, sft_33 = sft_as_31, sft_as_32, sft_as_33            
            

        d_map["out"] = f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33
        return f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33


    def plot_map(self):
        b_variable = self.is_variable()   
        d_map = {"flag": b_variable, "out":None}
        
        b_variable_phase = self.is_variable_phase()   
        d_phase = {"flag": b_variable_phase, "out":None}
        b_variable_adp = self.is_variable_adp()   
        d_adp = {"flag": b_variable_adp, "out":None}
        b_variable_magnetism = self.is_variable_magnetism()   
        d_magnetism = {"flag": b_variable_magnetism, "out":None}
        d_map.update({"phase": d_phase, "adp": d_adp, 
                      "magnetism": d_magnetism})
        
        return d_map
    
    def is_variable_phase(self):
        res = any([atom_type.is_variable_phase() for atom_type in 
                   self._list_atom_type])
        return res

    def is_variable_adp(self):
        res = any([atom_type.is_variable_adp() for atom_type in 
                   self._list_atom_type])
        return res

    def is_variable_magnetism(self):
        res = any([atom_type.is_variable_magnetism() for atom_type in 
                   self._list_atom_type])
        return res
    
    def is_variable(self):
        res = any([atom_type.is_variable() for atom_type in 
                   self._list_atom_type])
        return res
    
    def get_variables(self):
        l_variable = []
        for atom_type in self._list_atom_type:
            l_var = atom_type.get_variables()    
            l_variable.extend(l_var)
        return l_variable
    
    def apply_constraint(self, space_group, cell):
        for atom_type in self._list_atom_type:
            atom_type.apply_constraint(space_group, cell)
        self._form_arrays(cell)
