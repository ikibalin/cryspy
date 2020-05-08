import numpy 

def calc_mRmCmRT(r11, r12, r13, r21, r22, r23, r31, r32, r33,
                 c11, c12, c13, c21, c22, c23, c31, c32, c33):
    """
    calculate matrix multiplication R*C*RT, when matrices are expressed through 
    its component and can be expressed as nD-array
    """
    rc_11, rc_12 = r11*c11+r12*c21+r13*c31, r11*c12+r12*c22+r13*c32
    rc_13 = r11*c13+r12*c23+r13*c33
    rc_21, rc_22 = r21*c11+r22*c21+r23*c31, r21*c12+r22*c22+r23*c32
    rc_23 = r21*c13+r22*c23+r23*c33
    rc_31, rc_32 = r31*c11+r32*c21+r33*c31, r31*c12+r32*c22+r33*c32
    rc_33 = r31*c13+r32*c23+r33*c33

    #dimension (atoms, symmetry)
    rcrt_11 = (rc_11*r11+rc_12*r12+rc_13*r13)
    rcrt_12 = (rc_11*r21+rc_12*r22+rc_13*r23)
    rcrt_13 = (rc_11*r31+rc_12*r32+rc_13*r33)

    rcrt_21 = (rc_21*r11+rc_22*r12+rc_23*r13)
    rcrt_22 = (rc_21*r21+rc_22*r22+rc_23*r23)
    rcrt_23 = (rc_21*r31+rc_22*r32+rc_23*r33)

    rcrt_31 = (rc_31*r11+rc_32*r12+rc_33*r13)
    rcrt_32 = (rc_31*r21+rc_32*r22+rc_33*r23)
    rcrt_33 = (rc_31*r31+rc_32*r32+rc_33*r33)
    return rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, rcrt_32, rcrt_33



def calc_phase_by_hkl_xyz_rb(h, k, l, x, y, z, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2, b_3):
    np_h, np_x, np_r_11 = numpy.meshgrid(h, x, r_11, indexing="ij")
    np_k, np_y, np_r_22 = numpy.meshgrid(k, y, r_22, indexing="ij")
    np_l, np_z, np_r_33 = numpy.meshgrid(l, z, r_33, indexing="ij")

    np_r_12 = numpy.meshgrid(h, x, r_12, indexing="ij")[2]
    np_r_13 = numpy.meshgrid(k, y, r_13, indexing="ij")[2]
    np_r_23 = numpy.meshgrid(l, z, r_23, indexing="ij")[2]
    np_r_21 = numpy.meshgrid(h, x, r_21, indexing="ij")[2]
    np_r_31 = numpy.meshgrid(k, y, r_31, indexing="ij")[2]
    np_r_32 = numpy.meshgrid(l, z, r_32, indexing="ij")[2]
    np_b_1 = numpy.meshgrid(l, z, b_1, indexing="ij")[2]
    np_b_2 = numpy.meshgrid(l, z, b_2, indexing="ij")[2]
    np_b_3 = numpy.meshgrid(l, z, b_3, indexing="ij")[2]

    np_x_s = np_x*np_r_11 + np_y*np_r_12 + np_z*np_r_13 + np_b_1
    np_y_s = np_x*np_r_21 + np_y*np_r_22 + np_z*np_r_23 + np_b_2
    np_z_s = np_x*np_r_31 + np_y*np_r_32 + np_z*np_r_33 + np_b_3
    hh = (2*numpy.pi*1j*(np_h*np_x_s + np_k*np_y_s+ np_l*np_z_s)).astype(complex)
    phase_3d = numpy.exp(hh)
    return phase_3d


def calc_power_dwf_iso(b_iso, sthovl):
    """
isotropic harmonic Debye-Waller factor
    """
    sthovl_sq = sthovl**2
    b_iso_2d, sthovl_sq_2d = numpy.meshgrid(sthovl_sq, b_iso, indexing="ij")
    
    power_dwf_iso_2d = b_iso_2d*sthovl_sq_2d
    return power_dwf_iso_2d

def calc_power_dwf_aniso(h, k, l, beta, 
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33):
    """
anisotropic harmonic Debye-Waller factor

h,k,l is 1D (temporary solution)
    """
    b_11, b_22, b_33 = beta[:, 0], beta[:, 1], beta[:, 2]
    b_12, b_13, b_23 = beta[:, 3], beta[:, 4], beta[:, 5]

    np_h, np_b_11, np_r_11 = numpy.meshgrid(h, b_11, r_11, indexing="ij")
    np_k, np_b_22, np_r_22 = numpy.meshgrid(k, b_22, r_22, indexing="ij")
    np_l, np_b_33, np_r_33 = numpy.meshgrid(l, b_33, r_33, indexing="ij")
    np_h, np_b_12, np_r_12 = numpy.meshgrid(h, b_12, r_12, indexing="ij")
    np_h, np_b_13, np_r_13 = numpy.meshgrid(h, b_13, r_13, indexing="ij")
    np_h, np_b_23, np_r_23 = numpy.meshgrid(h, b_23, r_23, indexing="ij")
    np_r_21 = numpy.meshgrid(h, b_23, r_21, indexing="ij")[2]
    np_r_31 = numpy.meshgrid(h, b_23, r_31, indexing="ij")[2]
    np_r_32 = numpy.meshgrid(h, b_23, r_32, indexing="ij")[2]
    
    np_h_s = np_h*np_r_11 + np_k*np_r_21 + np_l*np_r_31
    np_k_s = np_h*np_r_12 + np_k*np_r_22 + np_l*np_r_32
    np_l_s = np_h*np_r_13 + np_k*np_r_23 + np_l*np_r_33
    
    power_dwf_aniso = (np_b_11*np_h_s**2 + np_b_22*np_k_s**2 + 
                       np_b_33*np_l_s**2 + 2.*np_b_12*np_h_s*np_k_s + 
                      2.*np_b_13*np_h_s*np_l_s + 2.*np_b_23*np_k_s*np_l_s)
    
    return power_dwf_aniso 
    
def calc_dwf(cell, h, k, l, b_iso, beta, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33):
    """
calculate Debye-Waller factor
    """
    sthovl = cell.calc_sthovl(h, k, l)
    #dimensions (hkl, atoms in assymmetric unit cell)
    power_iso_2d = calc_power_dwf_iso(b_iso, sthovl)

    #dimensions (hkl, atoms in assymmetric unit cell, el.symmetry)
    power_aniso_3d = calc_power_dwf_aniso(h, k, l, beta, 
                      r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
    power_3d = power_iso_2d[:, :, numpy.newaxis] + power_aniso_3d
    dwf_3d = numpy.exp(-power_3d)
    return dwf_3d



def calc_form_factor_tensor_susceptibility(chi_11, chi_22, chi_33, chi_12, chi_13, chi_23, space_group_symop, form_factor, cell, h, k, l):
    """
give components of form factor tensor for susceptibility:
    fft_11, fft_12, fft_13
    fft_21, fft_22, fft_23
    fft_31, fft_32, fft_33
    
in 3 dimension (hkl, atoms, symmetry elements)            
    """
    ff = numpy.array(form_factor, dtype=float)
    sthovl = cell.calc_sthovl(h, k, l)
    #dimension (hkl, atoms)

    r_11 = numpy.array(space_group_symop.r_11, dtype=float)
    r_12 = numpy.array(space_group_symop.r_12, dtype=float)
    r_13 = numpy.array(space_group_symop.r_13, dtype=float)
    r_21 = numpy.array(space_group_symop.r_21, dtype=float)
    r_22 = numpy.array(space_group_symop.r_22, dtype=float)
    r_23 = numpy.array(space_group_symop.r_23, dtype=float)
    r_31 = numpy.array(space_group_symop.r_31, dtype=float)
    r_32 = numpy.array(space_group_symop.r_32, dtype=float)
    r_33 = numpy.array(space_group_symop.r_33, dtype=float)
    
    chi_21, chi_31, chi_32 = chi_12, chi_13, chi_23 
    
    c11, r11 = numpy.meshgrid(chi_11, r_11, indexing="ij")
    c22, r22 = numpy.meshgrid(chi_22, r_22, indexing="ij")
    c33, r33 = numpy.meshgrid(chi_33, r_33, indexing="ij")
    c12, r12 = numpy.meshgrid(chi_12, r_12, indexing="ij")
    c13, r13 = numpy.meshgrid(chi_13, r_13, indexing="ij")
    c23, r23 = numpy.meshgrid(chi_23, r_23, indexing="ij")
    c21, r21 = numpy.meshgrid(chi_21, r_21, indexing="ij")
    c31, r31 = numpy.meshgrid(chi_31, r_31, indexing="ij")
    c32, r32 = numpy.meshgrid(chi_32, r_32, indexing="ij")
    
    rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, rcrt_32, rcrt_33 = calc_mRmCmRT(
            r11, r12, r13, r21, r22, r23, r31, r32, r33,
            c11, c12, c13, c21, c22, c23, c31, c32, c33)        

    #dimension (hkl, atoms, symmetry)
    fft_11 = ff[:,:,numpy.newaxis]*rcrt_11[numpy.newaxis, :,:]
    fft_12 = ff[:,:,numpy.newaxis]*rcrt_12[numpy.newaxis, :,:]
    fft_13 = ff[:,:,numpy.newaxis]*rcrt_13[numpy.newaxis, :,:]
    fft_21 = ff[:,:,numpy.newaxis]*rcrt_21[numpy.newaxis, :,:]
    fft_22 = ff[:,:,numpy.newaxis]*rcrt_22[numpy.newaxis, :,:]
    fft_23 = ff[:,:,numpy.newaxis]*rcrt_23[numpy.newaxis, :,:]
    fft_31 = ff[:,:,numpy.newaxis]*rcrt_31[numpy.newaxis, :,:]
    fft_32 = ff[:,:,numpy.newaxis]*rcrt_32[numpy.newaxis, :,:]
    fft_33 = ff[:,:,numpy.newaxis]*rcrt_33[numpy.newaxis, :,:]

    
    
    #ortogonalization should be done
    
    return fft_11, fft_12, fft_13, fft_21, fft_22, fft_23, fft_31, fft_32, fft_33


