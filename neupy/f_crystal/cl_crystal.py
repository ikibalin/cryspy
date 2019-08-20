"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_crystal import AbstractCrystal

from neupy.f_crystal.cl_space_group import SpaceGroup
from neupy.f_crystal.cl_cell import Cell
from neupy.f_crystal.cl_atom_site import AtomSite
from neupy.f_crystal.cl_extinction import Extinction

from neupy.f_crystal.cl_magnetism import calc_mRmCmRT
    
class Crystal(AbstractCrystal):
    """
    Crystal
    """
    def __init__(self, name=None, space_group=SpaceGroup(), cell=Cell(), 
                 atom_site=AtomSite(), extinction=Extinction(), i_g=0.):
        super(Crystal, self).__init__()
        
        self._p_name = None
        self._p_space_group = None
        self._p_cell = None
        self._p_atom_site = None
        self._p_extinction = None
        self._p_i_g = None
        self._refresh(name, space_group, cell, atom_site, extinction, i_g)
        
    def __repr__(self):
        lsout = """Crystal: \n name: {:}\n i_g: {:}\n{:}\n{:}\n{:}
{:}""".format(self._p_name, self._p_i_g, self._p_space_group, self._p_cell, 
                         self._p_atom_site, self._p_extinction)
        return lsout

    def _refresh(self, name, space_group, cell, atom_site, extinction, i_g):
        if name is not None:
            self._p_name = name
        if cell is not None:
            self._p_cell = cell
        if space_group is not None:
            self._p_space_group = space_group
            if self._p_cell is not None:
                self._p_cell.set_val(singony=space_group.get_val("singony"))
        if atom_site is not None:
            self._p_atom_site = atom_site
        if extinction is not None:
            self._p_extinction = extinction
        if i_g is not None:
            self._p_i_g = i_g

    def set_val(self, name=None, space_group=None, cell=None, atom_site=None, 
                extinction=None, i_g=None):
        self._refresh(name, space_group, cell, atom_site, extinction, i_g)
        
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
    
name is the name of the Crystal
space_group is the space group
cell is the unit cell parameters
atome_sie is the atom site
extinction is the extinction
i_g is the parameter to described broadening of Bragg reflection due to the 
       particle size
        """
        print(lsout)

    def load_from_rcif(self, f_name):
        from neupy import (RCif, conv_data_to_crystal)
        rcif = RCif(f_name)
        p_glob = rcif.glob
        l_data = p_glob["data"]
        l_crystal = conv_data_to_crystal(l_data)
        if len(l_crystal) != 0:
            cryst = l_crystal[0]

            self.set_val(name=cryst.get_val("name"), 
                         space_group=cryst.get_val("space_group"), 
                         cell=cryst.get_val("cell"), 
                         atom_site=cryst.get_val("atom_site"), 
                         extinction=cryst.get_val("extinction"), 
                         i_g=cryst.get_val("i_g"))
        else:
            print("In file '{:}' phase is not found".format(f_name))
        return

    def calc_fn(self, l_hkl, f_print=False):
        np_h = numpy.array([hh[0] for hh in l_hkl], dtype=int)
        np_k = numpy.array([hh[1] for hh in l_hkl], dtype=int)
        np_l = numpy.array([hh[2] for hh in l_hkl], dtype=int)
        f_nucl = self.calc_sf(np_h, np_k, np_l, f_print=f_print)[0]
        return f_nucl

    def calc_sf(self, h, k, l, d_map={}, f_print=False):
        """
        calculate structure factors (nuclear and components of susceptibility)

        output:

        f_nucl, s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33, d_info_out

        f_nucl is a nuclear structure factor
        s_ij is a component of structure factor tensor in Cartezian crystallographic system
        d_info_out is a dictionary with additional information
        """
        #if d_map == {}:
        #    d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    f_nucl, s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = d_map["out"]
        #    return f_nucl, s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33
        
        cell = self._p_cell
        space_group = self._p_space_group
        atom_site = self._p_atom_site
        #d_sf = d_map["sf"]
        #if not(d_sf["flag"]|(d_sf["out"] is None)):
        #    f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = d_sf["out"]
        #else:
        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = atom_site.calc_sf(space_group, cell, h, k, l)#, d_sf
        #sft_ij form the structure factor tensor in local coordinate system (ia, ib, ic)
        #chi in 10-12 cm; chim in muB (it is why here 0.2695)
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = self._orto_matrix(
                cell,
                sft_11*0.2695, sft_12*0.2695, sft_13*0.2695, 
                sft_21*0.2695, sft_22*0.2695, sft_23*0.2695, 
                sft_31*0.2695, sft_32*0.2695, sft_33*0.2695)
        d_info_out = {"h": h, "k": k, "l": l,
                      "f_nucl": f_nucl, "sft_11": sft_11, "sft_12": sft_12, 
                      "sft_13": sft_13, "sft_21": sft_21, "sft_22": sft_22, 
                      "sft_23": sft_23, "sft_31": sft_31, "sft_32": sft_32, 
                      "sft_33": sft_33}
        if f_print:
            ls_out = ["   h   k   l     Re(f_nucl)    Im(f_nucl)"]
            try:
                ls_out.extend(["{:4}{:4}{:4} {:14.3f}{:14.3f}".format(h_1, h_2, h_3, h_4, h_5) 
                    for h_1, h_2, h_3, h_4, h_5 in zip(h, k, l, f_nucl.real, f_nucl.imag)])
            except TypeError:
                ls_out.append("{:4}{:4}{:4} {:14.3f}{:14.3f}".format(h, k, l, float(f_nucl.real), float(f_nucl.imag)))
            print("\n".join(ls_out))
        return f_nucl, s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33, d_info_out


    def _orto_matrix(self, cell, l_11, l_12, l_13, l_21, l_22, l_23, l_31, 
                     l_32, l_33):
        """
        rewrite matrix l_ij defined in coordinate (ia, ib, ic) to matrix s_ij, 
        which is denined in Chartesian coordinate system, such as:
        x||ia, y in blane (ia, ib), z perpendicular to that plane.
        ...
        
        ...
        representation of chi in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        expressions are taken from international tables
        matrix_ib is inversed matrix B
        ia, ib, ic is inversed unit cell parameters (it can be estimated from matrix matrix_ib)

        X = B x, x = iB X
        xT*CHI*x = XT iBT CHI iB X
    
        output chiLOC = iBT CHI iB
        """
        m_ib = cell.get_val("m_ib")
        ia, ib, ic = cell.get_val("ia"), cell.get_val("ib"), cell.get_val("ic")
        """
        matrix_chi = numpy.array(
                [[self["chi_11"], self["chi_12"], self["chi_13"]],
                 [self["chi_12"], self["chi_22"], self["chi_23"]],
                 [self["chi_13"], self["chi_23"], self["chi_33"]]], 
                 dtype = float)
        #mchi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
        #[a,b,c,alpha,beta,gamma]=ucp
        y1 = m_ib[0,0]
        y2 = m_ib[1,1]
        y3 = m_ib[2,2]
        y4 = m_ib[0,1]
        y5 = m_ib[0,2]
        y6 = m_ib[1,2]
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
        """
        m_ib_norm = numpy.copy(m_ib)
        m_ib_norm[:,0] *= ia
        m_ib_norm[:,1] *= ib
        m_ib_norm[:,2] *= ic
        
        m_ibt_norm = m_ib_norm.transpose()
        
        r11, r12, r13 = m_ibt_norm[0, 0], m_ibt_norm[0, 1], m_ibt_norm[0, 2]
        r21, r22, r23 = m_ibt_norm[1, 0], m_ibt_norm[1, 1], m_ibt_norm[1, 2]
        r31, r32, r33 = m_ibt_norm[2, 0], m_ibt_norm[2, 1], m_ibt_norm[2, 2]
        
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = calc_mRmCmRT(
                r11, r12, r13, r21, r22, r23, r31, r32, r33,
                l_11, l_12, l_13, l_21, l_22, l_23, l_31, l_32, l_33)        
        """
        ibt_chi = numpy.matmul(m_ibt_norm, matrix_chi)
        matrix_chi_loc = numpy.matmul(ibt_chi, m_ib_norm)
        d_out = dict(matrix_chi_loc = matrix_chi_loc)
        self.update(d_out)
        """
        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33
    
    def calc_extinction(self, h, k, l, fp_sq, fm_sq, fpm_sq, wave_length):
        """
        extinction correction coefficients for polarized neutrons
        """
        cell = self._p_cell
        extinction = self._p_extinction
        yp = extinction.calc_extinction(cell, h, k, l, fp_sq, wave_length)
        ym = extinction.calc_extinction(cell, h, k, l, fm_sq, wave_length)
        ypm = extinction.calc_extinction(cell, h, k, l, fpm_sq, wave_length)
        return yp, ym, ypm
        
    def plot_map(self):
        b_variable = self.is_variable()         
        d_map = {"flag": b_variable, "out":None}
        atom_site = self.get_val("atom_site")
        d_sf = atom_site.plot_map()
        d_map.update({"sf":d_sf})
        return d_map
        
    def is_variable_sf(self):
        cell = self._p_cell 
        atom_site = self._p_atom_site 
        res = any([cell.is_variable(), atom_site.is_variable()])
        return res
    
    def is_variable_extinction(self):
        extinction = self._p_extinction 
        res = extinction.is_variable()
        return res
    
    def is_variable(self):
        """
        without extinction
        """
        i_g = self._p_i_g 
        res = any([self.is_variable_sf(), 
                   isinstance(i_g, Variable)])
        return res

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_i_g, Variable):
            l_variable.append(self._p_i_g)
            
        cell = self._p_cell 
        l_var = cell.get_variables()
        l_variable.extend(l_var)
        
        extinction = self._p_extinction 
        l_var = extinction.get_variables()
        l_variable.extend(l_var)
        
        atom_site = self.get_val("atom_site")
        l_var = atom_site.get_variables()
        l_variable.extend(l_var)

        return l_variable
    
    def apply_constraint(self):
        space_group = self.get_val("space_group")
        cell = self.get_val("cell")
        cell.apply_constraint()
        atom_site = self.get_val("atom_site")
        atom_site.apply_constraint(space_group, cell)
        
    def print_report(self):
        s_out = "{:}".format(self)
        return s_out
    
    def add_atom(self, atom):
        atom_site = self._p_atom_site
        atom_site._list_atom_type.append(atom)
        atom_site._flag_refresh = True
        #self._form_arrays(cell)
    
    def del_atom(self, ind):
        atom_site = self._p_atom_site
        atom_site._list_atom_type.pop(ind)
        atom_site._flag_refresh = True        

    def calc_atoms_in_cell(self):
        atom_site = self.get_val("atom_site")
        space_group = self.get_val("space_group")
        l_atom_type = atom_site._list_atom_type
        l_name = []
        l_frac = []
        for atom_type in l_atom_type:
            name = atom_type.get_val("name")
            x, y, z = 1.*atom_type.get_val("x"), 1.*atom_type.get_val("y"), 1.*atom_type.get_val("z")
            x_s, y_s, z_s, mult_a = space_group.calc_xyz_mult(x, y, z)
            l_name.append(name)
            l_frac.append([(hh_1, hh_2, hh_3) for hh_1, hh_2, hh_3 in zip(x_s, y_s, z_s)])
        return l_name, l_frac
        
if (__name__ == "__main__"):
  pass

