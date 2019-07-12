"""
define classes to describe calculated data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy

from f_crystal.cl_crystal import *
from f_common.cl_variable import *


class CalculatedDataSingle(dict):
    """
    Calculate the model data for single crystal in polarized neutron diffraction experiment
    """
    def __init__(self, name=None, field=1, orientation=numpy.array([[1., 0., 0.], 
                 [0., 1., 0.], [0., 0., 1.]], dtype=float)):
        """
        field is magnetic field in global coordinate system
        orientation is transfer matrix from local coordinate system to global one
        """
        super(CalculatedDataSingle, self).__init__()
        self._p_name = None
        self._p_field = None
        self._p_orientation = None
        orientation
        self._refresh(name, field, orientation)

    def __repr__(self):
        lsout = """CalculatedDataSingle:\n name: {}\n field: {:}
 orientation: {:}""".format(self._p_name,
                self._p_field, self._p_orientation)
        return lsout

    def _refresh(self, name, field, orientation):
        if name is not None:
            self._p_name = name
        if field is not None:
            self._p_field = field
        if orientation is not None:
            self._p_orientation = orientation
            
    def set_val(self, name=None, field=None, orientation=None):
        self._refresh(name, field, orientation)
        
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
name is the name of Calculated data
scale is the scale factor for crystal
field is magnetic field in global coordinate system
orientation is transfer matrix from local coordinate system to global one
        """
        print(lsout)
    
    def calc_iint_u_d_flip_ratio(self, h, k, l, beam_polarization, wave_length, crystal):
        """
        calculate the integral intensity for h, k, l reflections
        wave_length is needed only for extinction correction
        """
        field_z = 1.*self._p_field
        orientation = self._p_orientation
        
        field = numpy.array([0., 0., field_z], dtype=float)

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
        field_loc = numpy.matmul(m_u_d.transpose(), field)
        field_norm = ((field_loc**2).sum())**0.5
        
        p_u = 1.*beam_polarization.get_val("p_u")
        p_d = (2.*beam_polarization.get_val("flipper_efficiency")-1)*p_u
        
        e_u_loc = field_loc/field_norm


        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33, d_info_sf = crystal.calc_sf(h, k, l)
        
        
        cell = crystal.get_val("cell")
        
        k_1, k_2, k_3 = cell.calc_k_loc(h, k, l)
        
        
        mag_1 = sft_11*field_loc[0] + sft_12*field_loc[1] + sft_13*field_loc[2]
        mag_2 = sft_21*field_loc[0] + sft_22*field_loc[1] + sft_23*field_loc[2]
        mag_3 = sft_31*field_loc[0] + sft_32*field_loc[1] + sft_33*field_loc[2]
        
        #vector product k x mag x k        
        mag_p_1 = (k_3*mag_1 - k_1*mag_3)*k_3 - (k_1*mag_2 - k_2*mag_1)*k_2
        mag_p_2 = (k_1*mag_2 - k_2*mag_1)*k_1 - (k_2*mag_3 - k_3*mag_2)*k_3
        mag_p_3 = (k_2*mag_3 - k_3*mag_2)*k_2 - (k_3*mag_1 - k_1*mag_3)*k_1
        
        mag_p_sq = abs(mag_p_1*mag_p_1.conjugate() + 
                       mag_p_2*mag_p_2.conjugate() + 
                       mag_p_3*mag_p_3.conjugate())
        
        mag_p_e_u = mag_p_1*e_u_loc[0]+mag_p_2*e_u_loc[1]+mag_p_3*e_u_loc[2]

        
        f_nucl_sq = abs(f_nucl)**2
        mag_p_e_u_sq = abs(mag_p_e_u*mag_p_e_u.conjugate())
        fnp = (mag_p_e_u*f_nucl.conjugate()+mag_p_e_u.conjugate()*f_nucl).real
        fp_sq = f_nucl_sq + mag_p_sq + fnp
        fm_sq = f_nucl_sq + mag_p_sq - fnp
        fpm_sq = mag_p_sq - mag_p_e_u_sq


        #extinction correction        
        yp, ym, ypm =crystal.calc_extinction(h, k, l, fp_sq, fm_sq, fpm_sq, 
                                             wave_length)

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
        
        d_info_out = {"iint_u": iint_u, "iint_d": iint_d, 
                      "flip_ratio": flip_ratio}
        d_info_out.update(d_info_sf)      
        """
        print("   h   k   l  iint_u  iint_d flip_ratio")
        for h1, k1, l1, i_u, i_d, f_r in zip(h, k, l, iint_u, iint_d, flip_ratio):
            print("{:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, i_u, i_d, f_r))
        """
        return iint_u, iint_d, flip_ratio, d_info_out
    
    
    def is_variable(self):
        """
        without extinction
        """
        res = False
        return res        
    
    def get_variables(self):
        l_variable = []
        return l_variable


if (__name__ == "__main__"):
  pass

