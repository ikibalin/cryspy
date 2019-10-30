
"""
define class Diffrn which describes the single diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_09"
import os
import numpy
from pycifstar import Data, Loop

from cryspy.f_common.cl_fitable import Fitable
from cryspy.f_experiment.cl_beam_polarization import BeamPolarization
from cryspy.f_experiment.cl_orient_matrix import OrientMatrix
from cryspy.f_experiment.f_single.cl_diffrn_refln import DiffrnRefln
from cryspy.f_experiment.cl_extinction import Extinction



class Diffrn(object):
    """
    Class to describe information about single diffraction measurements
    """
    def __init__(self, label="mono", 
                 extinction = Extinction(), diffrn_refln=DiffrnRefln(), orient_matrix=OrientMatrix(),
                 beam_polarization=BeamPolarization(), 
                 wavelength=1.4, field=1.0, phase_label="phase1",
                 ):
        super(Diffrn, self).__init__()
        self.__label = None
        self.__extinction = None
        self.__diffrn_refln = None
        self.__orient_matrix = None
        self.__beam_polarization = None

        self.__diffrn_radiation_wavelength = None
        self.__diffrn_ambient_field = None
        self.__refln = None
        self.__phase_label = None

        self.label = label
        self.extinction = extinction
        self.diffrn_refln = diffrn_refln
        self.beam_polarization = beam_polarization
        self.wavelength = wavelength
        self.phase_label = phase_label
        
        orient_matrix.wavelength = wavelength
        self.orient_matrix = orient_matrix
        
        self.field = field

    @property
    def label(self):
        return self.__label
    @label.setter
    def label(self, x: str):
        self.__label = str(x)

    @property
    def phase_label(self):
        return self.__phase_label
    @phase_label.setter
    def phase_label(self, x: str):
        self.__phase_label = str(x)

    @property
    def extinction(self):
        return self.__extinction
    @extinction.setter
    def extinction(self, x):
        if isinstance(x, Extinction):
            x_in = x
        elif isinstance(x, str):
            x_in = Extinction()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Extinction")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for extinction is not recognized")
        self.__extinction = x_in

    @property
    def diffrn_refln(self):
        return self.__diffrn_refln
    @diffrn_refln.setter
    def diffrn_refln(self, x):
        if isinstance(x, DiffrnRefln):
            x_in = x
        elif isinstance(x, str):
            x_in = DiffrnRefln()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to DiffrnRefln")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for diffrn_refln is not recognized")
        self.__diffrn_refln = x_in

    @property
    def orient_matrix(self):
        return self.__orient_matrix
    @orient_matrix.setter
    def orient_matrix(self, x):
        if isinstance(x, OrientMatrix):
            x_in = x
        elif isinstance(x, str):
            x_in = OrientMatrix()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to OrientMatrix")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for orient_matrix is not recognized")
        self.__orient_matrix = x_in

    @property
    def beam_polarization(self):
        return self.__beam_polarization
    @beam_polarization.setter
    def beam_polarization(self, x):
        if isinstance(x, BeamPolarization):
            x_in = x
        elif isinstance(x, str):
            x_in = BeamPolarization()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to BeamPolarization")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for beam_polarization is not recognized")
        self.__beam_polarization = x_in

    @property
    def wavelength(self):
        return self.__diffrn_radiation_wavelength
    @wavelength.setter
    def wavelength(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__diffrn_radiation_wavelength = x_in

    @property
    def field(self):
        return self.__diffrn_ambient_field
    @field.setter
    def field(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__diffrn_ambient_field = x_in

    @property
    def refln(self):
        return self.__refln
    @refln.setter
    def refln(self, x):
        self.__refln = x


    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    def __repr__(self):
        ls_out = ["Diffrn"]
        if self.label is not None:
            ls_out.append("label: {:}".format(self.label))
        if self.wavelength is not None:
            ls_out.append("\nwavelength: {:.3f}".format(float(self.wavelength)))
        if self.field is not None:
            ls_out.append("field: {:.3f}".format(float(self.field)))
        if self.extinction is not None:
            ls_out.append("\n"+str(self.extinction))
        if self.beam_polarization is not None:
            ls_out.append("\n"+str(self.beam_polarization))
        if self.orient_matrix is not None:
            ls_out.append("\n"+str(self.orient_matrix))
        if self.diffrn_refln is not None:
            ls_out.append("\n"+str(self.diffrn_refln))
        if self.refln is not None:
            ls_out.append("\n"+str(self.refln))
        return "\n".join(ls_out)

    
    def calc_iint_u_d_flip_ratio(self, h, k, l, l_crystal):
        """
        calculate intensity for the given diffraction angle
        """
        l_label = [_.label for _ in l_crystal]
        if self.phase_label in l_label:
            ind = l_label.index(self.phase_label)
            crystal = l_crystal[ind]
        elif self.label in l_label:
            self.phase_label = self.label
            ind = l_label.index(self.label)
            crystal = l_crystal[ind]
        else:
            self._show_message("Crystal not found.\nThe first one is taken.")
            crystal = l_crystal[0]
            self.phase_label = crystal.label
        wavelength = float(self.wavelength)
        field_z = self.field
        beam_polarization = self.beam_polarization
        extinction = self.extinction
        orient_matrix = self.orient_matrix

        orientation = orient_matrix.u 
        
        field_vec = numpy.array([0., 0., field_z], dtype=float)

        phi_d, chi_d, omega_d = 0., 0., 0.
        m_phi_d = numpy.array([[ numpy.cos(phi_d), numpy.sin(phi_d),0.],
                               [-numpy.sin(phi_d), numpy.cos(phi_d),0.],
                               [         0.,       0.,    1.]], dtype=float)

        m_omega_d = numpy.array([[ numpy.cos(omega_d), numpy.sin(omega_d),0.],
                                 [-numpy.sin(omega_d), numpy.cos(omega_d),0.],
                                 [       0.,       0.,    1.]], dtype=float)

        m_chi_d = numpy.array([[ numpy.cos(chi_d), 0., numpy.sin(chi_d)],
                               [               0., 1.,               0.],
                               [-numpy.sin(chi_d), 0., numpy.cos(chi_d)]], dtype=float)
        
        m_u_d = numpy.matmul(m_omega_d, numpy.matmul(m_chi_d, numpy.matmul(m_phi_d, 
                             orientation)))

        field_loc = numpy.matmul(m_u_d.transpose(), field_vec)
        field_norm = ((field_loc**2).sum())**0.5
        
        p_u = float(beam_polarization.polarization)
        p_d = (2.*float(beam_polarization.efficiency)-1.)*p_u
        
        e_u_loc = field_loc/field_norm

        refln = crystal.calc_sf(h, k, l)
        f_nucl, sft_11, sft_12, sft_13 = refln.f_nucl, refln.sft_11, refln.sft_12, refln.sft_13
        sft_21, sft_22, sft_23 = refln.sft_21, refln.sft_22, refln.sft_23
        sft_31, sft_32, sft_33 = refln.sft_31, refln.sft_32, refln.sft_33
        
        cell = crystal.cell
        k_1, k_2, k_3 = cell.calc_k_loc(h, k, l)
        
        mag_1 = sft_11*field_loc[0] + sft_12*field_loc[1] + sft_13*field_loc[2]
        mag_2 = sft_21*field_loc[0] + sft_22*field_loc[1] + sft_23*field_loc[2]
        mag_3 = sft_31*field_loc[0] + sft_32*field_loc[1] + sft_33*field_loc[2]
        
        #vector product k x mag x k        
        mag_p_1 = (k_3*mag_1 - k_1*mag_3)*k_3 - (k_1*mag_2 - k_2*mag_1)*k_2
        mag_p_2 = (k_1*mag_2 - k_2*mag_1)*k_1 - (k_2*mag_3 - k_3*mag_2)*k_3
        mag_p_3 = (k_2*mag_3 - k_3*mag_2)*k_2 - (k_3*mag_1 - k_1*mag_3)*k_1
        
        mag_p_sq = abs(mag_p_1*mag_p_1.conjugate() + mag_p_2*mag_p_2.conjugate() + 
                       mag_p_3*mag_p_3.conjugate())
        
        mag_p_e_u = mag_p_1*e_u_loc[0]+mag_p_2*e_u_loc[1]+mag_p_3*e_u_loc[2]
        
        f_nucl_sq = abs(f_nucl)**2
        mag_p_e_u_sq = abs(mag_p_e_u*mag_p_e_u.conjugate())
        fnp = (mag_p_e_u*f_nucl.conjugate()+mag_p_e_u.conjugate()*f_nucl).real
        fp_sq = f_nucl_sq + mag_p_sq + fnp
        fm_sq = f_nucl_sq + mag_p_sq - fnp
        fpm_sq = mag_p_sq - mag_p_e_u_sq

        #extinction correction     
        yp = extinction.calc_extinction(cell, h, k, l, fp_sq, wavelength)
        ym = extinction.calc_extinction(cell, h, k, l, fm_sq, wavelength)
        ypm = extinction.calc_extinction(cell, h, k, l, fpm_sq, wavelength)

        pppl = 0.5*((1+p_u)*yp+(1-p_u)*ym)
        ppmin= 0.5*((1-p_d)*yp+(1+p_d)*ym)
        pmpl = 0.5*((1+p_u)*yp-(1-p_u)*ym)
        pmmin= 0.5*((1-p_d)*yp-(1+p_d)*ym)
        """
        print("   h   k   l  f_nucl f_m_p_sq   f_np  fpm_sq")
        for h1, k1, l1, f_n_sq, f_m_sq, f_np, f_pm_sq in zip(h, k, l, f_nucl_sq, mag_p_e_u_sq, fnp, fpm_sq):
            print(" {:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, f_n_sq, f_m_sq, f_np, f_pm_sq))

        print("   h   k   l      yp      ym     ypm      pppl   ppmin    pmpl   pmmin")
        for h1, k1, l1, y_p, y_m, y_pm, p_ppl, p_pmin, p_mpl, p_mmin in zip(
                h, k, l, yp, ym, ypm, pppl, ppmin, pmpl, pmmin):
            print(" {:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f}   {:7.3f} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, y_p, y_m, y_pm, p_ppl, p_pmin, p_mpl, p_mmin))
        """
        #integral intensities and flipping ratios
        iint_u = (f_nucl_sq+mag_p_e_u_sq)*pppl + pmpl*fnp + ypm*fpm_sq
        iint_d = (f_nucl_sq+mag_p_e_u_sq)*ppmin + pmmin*fnp + ypm*fpm_sq

        flip_ratio = iint_u/iint_d
        
        """
        d_info_out = {"iint_u": iint_u, "iint_d": iint_d, 
                      "flip_ratio": flip_ratio}
        d_info_out.update(d_info_sf)      
        print("   h   k   l  iint_u  iint_d flip_ratio")
        for h1, k1, l1, i_u, i_d, f_r in zip(h, k, l, iint_u, iint_d, flip_ratio):
            print("{:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, i_u, i_d, f_r))
        """
            
        return iint_u, iint_d, flip_ratio, refln
    
    def calc_chi_sq(self, l_crystal):
        """
        calculate chi square
        """
        diffrn_refln = self.diffrn_refln
        h, k, l = diffrn_refln.h, diffrn_refln.k, diffrn_refln.l
        fr_exp = diffrn_refln.fr
        fr_sigma = diffrn_refln.fr_sigma
        int_u_mod, int_d_mod, fr_mod, refln = self.calc_iint_u_d_flip_ratio(
                                               h, k, l, l_crystal)
        self.diffrn_refln.fr_calc = fr_mod
        self.diffrn_refln.intensity_up_calc = int_u_mod
        self.diffrn_refln.intensity_down_calc = int_d_mod

        self.refln = refln

        chi_sq = ((fr_mod-fr_exp)/fr_sigma)**2
        chi_sq_val = (chi_sq[numpy.logical_not(numpy.isnan(chi_sq))]).sum()
        n = numpy.logical_not(numpy.isnan(chi_sq)).sum()
        return chi_sq_val, n
    
    @property
    def is_variable(self):
        """
        without extinction
        """
        l_bool = []
        if self.beam_polarization is not None:
            l_bool.append(self.beam_polarization.is_variable)
        if self.extinction is not None:
            l_bool.append(self.extinction.is_variable)
        if self.wavelength.refinement:
            l_bool.append(self.wavelength)
        res = any(l_bool)
        return res

    def get_variables(self):
        l_variable = []
        if self.beam_polarization is not None:
            l_variable.extend(self.beam_polarization.get_variables())
        if self.extinction is not None:
            l_variable.extend(self.extinction.get_variables())
        if self.wavelength.refinement: l_variable.append(self.wavelength)
        return l_variable


    @property
    def to_cif(self):
        ls_out = []
        ls_out.append("data_{:}".format(self.label))
        str_1 = self.params_to_cif
        if str_1 != "":
            ls_out.append(str_1)
        str_2 = self.data_to_cif
        if str_2 != "":
            ls_out.append(str_2)
        str_3 = self.calc_to_cif
        if str_3 != "":
            ls_out.append(str_3)
        return "\n".join(ls_out)

    @property
    def params_to_cif(self):
        ls_out = []
        if self.extinction is not None:
            ls_out.append("\n"+self.extinction.to_cif)
        if self.beam_polarization is not None:
            ls_out.append("\n"+self.beam_polarization.to_cif)
        if self.phase_label is not None:
            ls_out.append("\n"+"_diffrn_phase_label {:}".format(self.phase_label))
        if self.wavelength is not None:
            ls_out.append("\n"+"_diffrn_radiation_wavelength {:}".format(self.wavelength.print_with_sigma))
        if self.field is not None:
            ls_out.append("_diffrn_ambient_field {:}".format(self.field))
        if self.orient_matrix is not None:
            ls_out.append("\n"+self.orient_matrix.to_cif)
        return "\n".join(ls_out)

    @property
    def data_to_cif(self):
        ls_out = []
        if self.diffrn_refln is not None:
            ls_out.append("\n"+self.diffrn_refln.to_cif)
        return "\n".join(ls_out)

    @property
    def calc_to_cif(self):
        ls_out = []
        #it is doubled with data_to_cif as there is calculated values.  Should it be separated???
        if self.diffrn_refln is not None:
            ls_out.append("\n"+self.diffrn_refln.to_cif)
        if self.refln is not None:
            ls_out.append("\n"+self.refln.to_cif)
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        self.label = cif_data.name
        cif_values = cif_data.items
        if cif_values is not None:
            if cif_values.is_prefix("_diffrn_phase_label"):
                self.phase_label = cif_values["_diffrn_phase_label"]
            if cif_values.is_prefix("_diffrn_radiation_wavelength"):
                self.wavelength = cif_values["_diffrn_radiation_wavelength"]
            if cif_values.is_prefix("_diffrn_ambient_field"):
                self.field = float(cif_values["_diffrn_ambient_field"])
            if cif_values.is_prefix("_diffrn_radiation"):
                self.beam_polarization = str(cif_values)
            if cif_values.is_prefix("_diffrn_orient_matrix_UB"):
                self.orient_matrix = str(cif_values)
            if cif_values.is_prefix("_refine_ls_extinction"):
                self.extinction = str(cif_values)
        if cif_data.is_prefix("_diffrn_refln"): self.diffrn_refln = str(cif_data["_diffrn_refln"])
        return True

