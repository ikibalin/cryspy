"""
define classes to describe calculated data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy

from crystal import *
from variable import *


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
        
        p_u = beam_polarization.get_val("p_u")
        p_d = beam_polarization.get_val("p_d")
        
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



class CalculatedDataPowder1D(dict):
    """
    Calculate the model data for 1D powder diffraction experiment
    """
    def __init__(self, name=None, scale=1., field=1.):
        super(CalculatedDataPowder1D, self).__init__()
        self._p_name = None
        self._p_scale = None
        self._p_field = None
        self._refresh(name, scale, field)

    def __repr__(self):
        lsout = """CalculatedDataPowder1D: \n name: {:}\n scale: {:}
 field: {:}""".format(self._p_name, self._p_scale, self._p_field)
        return lsout

    def _refresh(self, name, scale, field):
        if name is not None:
            self._p_name = name
        if scale is not None:
            self._p_scale = scale
        if field is not None:
            self._p_field = field
            
    def set_val(self, name=None, scale=None, field=None):
        self._refresh(name, scale, field)
        
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
name is the name of CalculatedData1D
scale is the scale factor for crystal
field is the value of magnetic field applied along vertical direction in Tesla
crystal is the definition of crystal 
        """
        print(lsout)
    
    def calc_iint(self, h, k, l, beam_polarization, crystal):
        """
        calculate the integral intensity for h, k, l reflections
        """
        field = 1.*self._p_field
        p_u = beam_polarization.get_val("p_u")
        p_d = beam_polarization.get_val("p_d")


        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33, d_info_cry = crystal.calc_sf(h, k, l)
        
        cell = crystal.get_val("cell")
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, 
                sft_33)
        fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        fm_p_field = field*0.5*(th_11+th_22) 
        cross = 2.*(f_nucl.real*fm_p_field.real+f_nucl.imag*fm_p_field.imag)
        #lkloc=[cfunc.calck(hkl,mB) for hkl in lhkl]

        iint_u = abs(f_nucl*f_nucl.conjugate()) + fm_p_sq + p_u*cross
                 

        iint_d = abs(f_nucl*f_nucl.conjugate()) + fm_p_sq - p_d*cross
        d_info_out = {"iint_u": iint_u, "iint_d": iint_d}   
        d_info_out.update(d_info_cry)
        
        #I_p, I_m = self.calc_extinc_powder(h, k, l, fn, fm_perp_eup, fm_p_sq,ext, p_up, p_down, ucp, wave_length)
        #
        #print("   h   k   l fn_real fn_imag")
        #for h1, k1, l1, f, c11, c12, c13, c21, c22, c23, c31, c32, c33 in zip(h, k, l, f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33):
        #    print(""" {:3} {:3} {:3} {:7.3f} {:7.3f}
        #    {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i
        #    {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i
        #    {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i""".format(
        #            h1, k1, l1, f.real, f.imag, c11.real, c11.imag, c12.real, 
        #            c12.imag, c13.real, c13.imag, c21.real, c21.imag, c22.real, 
        #            c22.imag, c23.real, c23.imag, c31.real, c31.imag, c32.real, 
        #            c32.imag, c33.real, c33.imag))
        return iint_u, iint_d, d_info_out
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_scale, Variable)])
        return res        

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_scale, Variable):
            l_variable.append(self._p_scale)
        return l_variable

    
class CalculatedDataPowder2D(dict):
    """
    Calculate the model data for 2D powder diffraction experiment
    """
    def __init__(self, name=None, scale=1., field=1.):
        super(CalculatedDataPowder2D, self).__init__()
        self._p_name = None
        self._p_scale = None
        self._p_field = None
        self._refresh(name, scale, field)

    def __repr__(self):
        lsout = """CalculatedDataPowder2D:\n name: {:}\n scale: {:}
 field: {:}""".format(self._p_name, self._p_scale, self._p_field)
        return lsout

    def _refresh(self, name, scale, field):
        if name is not None:
            self._p_name = name
        if scale is not None:
            self._p_scale = scale
        if field is not None:
            self._p_field = field
            
    def set_val(self, name=None, scale=None, field=None):
        self._refresh(name, scale, field)
        
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
name is the name of CalculatedDataPowder2D
scale is the scale factor for crystal
field is the value of magnetic field applied along vertical direction in Tesla
        """
        print(lsout)
    
    def calc_for_iint(self, h, k, l, crystal, d_map={}):
        """
        calculate the integral intensity for h, k, l reflections
        Output: f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = d_map["out"]
        #    return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin
        
        
        field = 1.*self._p_field

        #d_sf = d_map["sf"]
        #if not(d_sf["flag"]|(d_sf["out"] is None)):
        #    f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = d_sf["out"]
        #else:
        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33, d_info_cry = crystal.calc_sf(h, k, l)
        
        cell = crystal.get_val("cell")
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, 
                sft_33)
        f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
        f_m_p_sin_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        f_m_p_cos_sq = (field**2)*abs(th_13*th_13.conjugate()+th_23*th_23.conjugate())
        f_m_p_field = 0.5*field*(th_11+th_22) 
        cross_sin = 2.*(f_nucl.real*f_m_p_field.real+f_nucl.imag*f_m_p_field.imag)
        #print("   h   k   l f_nucl_sq f_m_p_sin_sq f_m_p_cos_sq cross_sin")
        #for h_1, k_1, l_1, hh_1, hh_2, hh_3, hh_4  in zip(h, k, l, f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin):
        #    print(""" {:3} {:3} {:3} {:9.3f} {:9.3f} {:9.3f} {:9.3f}""".format(
        #            h_1, k_1, l_1, hh_1, hh_2, hh_3, hh_4))        
        d_info_out = {"f_nucl_sq": f_nucl_sq, "f_m_p_sin_sq": f_m_p_sin_sq, 
                      "f_m_p_cos_sq": f_m_p_cos_sq, "cross_sin": cross_sin}   
        d_info_out.update(d_info_cry)        
        return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, d_info_out 

    def plot_map(self):
        b_variable = self.is_variable()       
        d_map = {"flag": b_variable, "out":None}
        return d_map

    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_scale, Variable)])
        return res        
    
    def get_variables(self):
        l_variable = []
        if isinstance(self._p_scale, Variable):
            l_variable.append(self._p_scale)
        return l_variable
    
if (__name__ == "__main__"):
  pass

