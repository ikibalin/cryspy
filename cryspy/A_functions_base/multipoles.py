import numpy
from cryspy.A_functions_base.matrix_operations import \
    calc_m_v, calc_vector_product_v1_v2, calc_m1_m2_inv_m1, calc_inv_m
from cryspy.A_functions_base.orbital_functions import calc_rho_normalized
from cryspy.A_functions_base.symmetry_elements import apply_symm_elem_to_xyz
from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters

INV_PI = 1./numpy.pi
INV_PI_SQRT = numpy.sqrt(INV_PI)

A_PLUS = numpy.sqrt((30+numpy.sqrt(480))/70.)
A_MINUS = numpy.sqrt((30-numpy.sqrt(480))/70.)
N_ANG = 0.5*INV_PI/(
    14*(numpy.power(A_MINUS, 5)-numpy.power(A_PLUS, 5))-
    20*(numpy.power(A_MINUS, 3)-numpy.power(A_PLUS, 3))-
    6*(A_MINUS-A_PLUS))

# Order lm:
# 00, 11+, 11-, 10, 20 ,21+, 21-, 22+, 22-, 30, 31+, 31-, 32+, 32-, 33+, 33-, 40, 41+, 41-, 42+, 42-, 43+, 43-, 44+, 44-
NORMALIZATION_DENSITY_FUNCTION = numpy.array([
    0.25*INV_PI, # l = 0
    INV_PI, INV_PI, INV_PI, # l = 1
    0.125*3*numpy.sqrt(3)*INV_PI, 0.75, 0.75, 0.75, 0.75, # l = 2
    10.*INV_PI/13., 1./(numpy.arctan(2.)+2.8-0.25*numpy.pi), 1./(numpy.arctan(2.)+2.8-0.25*numpy.pi), 1., 1., 4.*INV_PI/3., 4.*INV_PI/3., # l = 3
    N_ANG, 735./(512*numpy.sqrt(7)-196), 735./(512*numpy.sqrt(7)-196), 0.25*105.*numpy.sqrt(7.)/(136.+28.*numpy.sqrt(7.)), 0.25*105.*numpy.sqrt(7.)/(136.+28.*numpy.sqrt(7.)), 1.25, 1.25, 15./32., 15./32. # l = 4
], dtype=float)

NORMALIZATION_WAFE_FUNCTION = numpy.sqrt(INV_PI*numpy.array([
    0.25, # l = 0
    0.75, 0.75, 0.75, # l = 1
    5./16., 15./4., 15./4., 15./4., 15./4., # l = 2
    7./16., 21./32., 21./32., 105./16., 105./16., 35./32., 35./32.,# l = 3
    9./256., 45./32., 45./32., 45./64., 45./64., 315./32., 315./32., 315./256., 315./256. # l = 4
], dtype=float))


# def calc_real_spherical_harmonic_functions(Plm, xyz, l_max:int = 4, flag_density_function: bool=True):
#     """
#     Order lm in Plm parameters:
#     00, 11+, 11-, 10, 20 ,21+, 21-, 22+, 22-, 30, 31+, 31-, 32+, 32-, 33+, 33-, 40, 41+, 41-, 42+, 42-, 43+, 43-, 44+, 44-
    
#     Calculate 
#     sum_{m=-l}^{l} P_lm y_lm for l = 0, 1, 2, 3, 4
    
#     Output size [l_max, n_xyz]
#     Order of Plm parameters 
#     """
#     n_xyz = tuple([hh for hh in range(1, len(xyz.shape), 1)])

#     n_lm = Plm.shape[0]
#     if flag_density_function:
#         coeff = NORMALIZATION_DENSITY_FUNCTION
#     else:
#         coeff = NORMALIZATION_WAFE_FUNCTION
#     y_l_0 = numpy.stack([numpy.ones(shape=xyz.shape[1:], dtype=float)], axis=0)
#     res_l_0 = numpy.sum(numpy.expand_dims(coeff[:1]*Plm[:1], axis=n_xyz) * y_l_0, axis=0)
#     l_res = [res_l_0]
#     if l_max >= 1:
#         x, y, z = xyz[0], xyz[1], xyz[2]
#         y_l_1 = numpy.stack([x, y, z], axis=0)
#         res_l_1 = numpy.sum(numpy.expand_dims(coeff[1:4]*Plm[1:4], axis=n_xyz) * y_l_1, axis=0)
#         l_res.append(res_l_1)
#     if l_max >= 2:
#         x_sq, y_sq, z_sq = numpy.square(x), numpy.square(y), numpy.square(z)
#         x_sq_minus_y_sq = x_sq - y_sq
#         y_l_2 = numpy.stack([2.*z_sq - x_sq - y_sq, z*x, z*y, x*y, 0.5*x_sq_minus_y_sq], axis=0)
#         res_l_2 = numpy.sum(numpy.expand_dims(coeff[4:9]*Plm[4:9], axis=n_xyz) * y_l_2, axis=0)
#         l_res.append(res_l_2)
#     if l_max >= 3:
#         x_cube, y_cube, z_cube = x_sq*x, y_sq*y, z_sq*z
#         x_sq_plus_y_sq = x_sq + y_sq
#         y_l_3 = numpy.stack([
#             2.*z_cube - 3.*z*x_sq_plus_y_sq, 
#             x*(4.*z_sq-x_sq_plus_y_sq),
#             y*(4.*z_sq-x_sq_plus_y_sq),
#             z*x_sq_minus_y_sq,
#             2.*x*y*z,
#             x_cube-3.*x*y_sq,
#             y_cube-3.*y*x_sq,
#         ], axis=0)
#         res_l_3 = numpy.sum(numpy.expand_dims(coeff[9:16]*Plm[9:16], axis=n_xyz) * y_l_3, axis=0)
#         l_res.append(res_l_3)
#     if l_max >= 4:
#         y_l_4 = numpy.stack([
#             8.*z_sq*z_sq-24.*z_sq*x_sq_plus_y_sq + 3.*x_sq_plus_y_sq*x_sq_plus_y_sq,
#             x*(4.*z_cube-3.*z*x_sq_plus_y_sq),
#             y*(4.*z_cube-3.*z*x_sq_plus_y_sq),
#             x_sq_minus_y_sq*(6.*z_sq-x_sq_plus_y_sq),
#             2*x*y*(6.*z_sq-x_sq_plus_y_sq),
#             z*y_l_3[5], 
#             z*y_l_3[6],
#             x_sq*x_sq-6.*x_sq*y_sq+y_sq*y_sq,
#             4.*(x_cube*y-x*y_cube)
#         ], axis=0)
#         res_l_4 = numpy.sum(numpy.expand_dims(coeff[16:25]*Plm[16:25], axis=n_xyz) * y_l_4, axis=0)
#         l_res.append(res_l_4)
#     res = numpy.stack(l_res, axis=0)
#     return res




def realsphharm2df(n,m,theta,phi):
    """Real spherical harmonics taken from [International tables] normalized for density functions\n
    n is order or real spherical harmonics up to 4\n
    m is in range [-n:n]\n
    theta  -- polar angle in range [0:180]\n
    phi -- azimuthal angle in range [0:360]"""
    one = numpy.ones_like(theta)
    x=numpy.sin(theta)*numpy.cos(phi)
    y=numpy.sin(theta)*numpy.sin(phi)
    z=numpy.cos(theta)
    if (n,m)==(0,0):                                                          
        res=one                       *(1./(4.*numpy.pi))
    elif (n,m)==(1,1):
        res=x                         *(1./(numpy.pi))
    elif (n,m)==(1,-1):
        res=y                         *(1./(numpy.pi))
    elif (n,m)==(1,0):
        res=z                         *(1./(numpy.pi))
    elif (n,m)==(2,0):
        res=(3.*z**2-1.)              *((3.*3.**0.5)/(8.*numpy.pi))
    elif (n,m)==(2,1):
        res=(z*x)                     *(3./4.)
    elif (n,m)==(2,-1):
        res=(z*y)                     *(3./4.)
    elif (n,m)==(2,2):
        res=0.5*(x**2-y**2)           *(3./4.)
    elif (n,m)==(2,-2):
        res=(x*y)                     *(3./4.)
    elif (n,m)==(3,0):
        res=(5.*z**3-3.*z)            *(10./(13.*numpy.pi))
    elif (n,m)==(3,1):
        res=x*(5.*z**2-1.)            *(0.32033)
    elif (n,m)==(3,-1):
        res=y*(5.*z**2-1.)            *(0.32033)
    elif (n,m)==(3,2):
        res=z*(x**2-y**2)             *(1.)
    elif (n,m)==(3,-2):
        res=2*x*y*z                   *(1.)
    elif (n,m)==(3,3):
        res=(x**3-3*x*y**2)           *(4./(3.*numpy.pi))
    elif (n,m)==(3,-3):
        res=(-y**3+3*y*x**2)          *(4./(3.*numpy.pi))
    elif (n,m)==(4,0):
        res=(35.*z**4-30.*z**2+3.)    *(0.06942)
    elif (n,m)==(4,1):
        res=x*(7.*z**3-3.*z)          *(735./(512.*7.**0.5+196.))
    elif (n,m)==(4,-1):
        res=y*(7.*z**3-3.*z)          *(735./(512.*7.**0.5+196.))
    elif (n,m)==(4,2):
        res=(x**2-y**2)*(7.*z**2-1.)  *((105.*7.**0.5)/(4.*(136.+28.*7.**0.5)))
    elif (n,m)==(4,-2):
        res=2*x*y*(7.*z**2-1.)        *((105.*7.**0.5)/(4.*(136.+28.*7.**0.5)))
    elif (n,m)==(4,3):
        res=z*(x**3-3*x*y**2)         *(5./4.)
    elif (n,m)==(4,-3):
        res=z*(-y**3+3*y*x**2)        *(5./4.)
    elif (n,m)==(4,4):
        res=(x**4-6.*x**2*y**2+y**4)  *(15./32.)
    elif (n,m)==(4,-4):
        res=(4.*x**3*y-4.*x*y**3)     *(15./32.)
    else:
        res=[]
    return res

def realsphharm2wfXYZ(n,m,dcosxyz):
    """Real spherical harmonics taken from [International tables] normalized for wave functions\n
    n is order or real spherical harmonics up to 4
    m is in range [-n:n]
    theta  -- polar angle in range [0:180]\n
    phi -- azimuthal angle in range [0:360]"""
    x,y,z=dcosxyz[0],dcosxyz[1],dcosxyz[2]
    one = numpy.ones_like(x)
    if ((n,m)==(0,0)):
        res=one                        /((4.*numpy.pi)**0.5)
    elif (n==1):
        if (m==1):
            res=x                         *(3./(4*numpy.pi))**0.5
        elif (m==-1):
            res=y                         *(3./(4*numpy.pi))**0.5
        elif (m==0):
            res=z                         *(3./(4*numpy.pi))**0.5
        else:
            res=None
    elif (n==2):
        if (m==0):
            res=(3.*z**2-1.)              *(5./(16*numpy.pi))**0.5
        elif (m==1):
            res=(z*x)                     *(15./(4.*numpy.pi))**0.5
        elif (m==-1):
            res=(z*y)                     *(15./(4.*numpy.pi))**0.5
        elif (m==2):
            res=0.5*(x**2-y**2)           *(15./(4.*numpy.pi))**0.5
        elif (m==-2):
            res=(x*y)                     *(15./(4.*numpy.pi))**0.5
        else:
            res=None
    elif (n==3):
        if (m==0):
            res=(5.*z**3-3.*z)            *(7./(16.*numpy.pi))**0.5
        elif (m==1):
            res=x*(5.*z*z-1.)            *(21./(32.*numpy.pi))**0.5
        elif (m==-1):
            res=y*(5.*z*z-1.)            *(21./(32.*numpy.pi))**0.5
        elif (m==2):
            res=z*(x*x-y*y)             *(105./(16.*numpy.pi))**0.5
        elif (m==-2):
            res=2.*x*y*z                   *(105./(16.*numpy.pi))**0.5
        elif (m==3):
            res=(x*x*x-3.*x*y*y)           *(35./(32.*numpy.pi))**0.5
        elif (m==-3):
            res=(-y*y*y+3.*y*x*x)          *(35./(32.*numpy.pi))**0.5
        else:
            res=None
    elif (n==4):
        if (m==0):
            res=(35.*z*z*z*z-30.*z*z+3.)    *(9./(256.*numpy.pi))**0.5
        elif (m==1):
            res=x*(7.*z*z*z-3.*z)          *(45./(32.*numpy.pi))**0.5
        elif (m==-1):
            res=y*(7.*z*z*z-3.*z)          *(45./(32.*numpy.pi))**0.5
        elif (m==2):
            res=(x*x-y*y)*(7.*z*z-1.)  *(45./(64.*numpy.pi))**0.5
        elif (m==-2):
            res=2*x*y*(7.*z*z-1.)        *(45./(64.*numpy.pi))**0.5
        elif (m==3):
            res=z*(x*x*x-3.*x*y*y)         *(315./(32.*numpy.pi))**0.5
        elif (m==-3):
            res=z*(-y*y*y+3.*y*x*x)        *(315./(32.*numpy.pi))**0.5
        elif (m==4):
            res=(x*x*x*x-6.*x*x*y*y+y*y*y*y)  *(315./(256.*numpy.pi))**0.5
        elif (m==-4):
            res=(4.*x*x*x*y-4.*x*y*y*y)     *(315./(256.*numpy.pi))**0.5
        else:
            res=None
    else:
        res=None
    return res

def realsphharm2wf(n,m,theta,phi):
    """Real spherical harmonics taken from [International tables] normalized for wave functions\n
    n is order or real spherical harmonics up to 4
    m is in range [-n:n]
    theta  -- polar angle in range [0:180]\n
    phi -- azimuthal angle in range [0:360]"""
    x=numpy.sin(theta)*numpy.cos(phi)
    y=numpy.sin(theta)*numpy.sin(phi)
    z=numpy.cos(theta)
    one = numpy.ones_like(x)
    if ((n,m)==(0,0)):
        res=one                        /((4.*numpy.pi)**0.5)
    elif (n==1):
        if (m==1):
            res=x                         *(3./(4*numpy.pi))**0.5
        elif (m==-1):
            res=y                         *(3./(4*numpy.pi))**0.5
        elif (m==0):
            res=z                         *(3./(4*numpy.pi))**0.5
        else:
            res=None
    elif (n==2):
        if (m==0):
            res=(3.*z**2-1.)              *(5./(16*numpy.pi))**0.5
        elif (m==1):
            res=(z*x)                     *(15./(4.*numpy.pi))**0.5
        elif (m==-1):
            res=(z*y)                     *(15./(4.*numpy.pi))**0.5
        elif (m==2):
            res=0.5*(x**2-y**2)           *(15./(4.*numpy.pi))**0.5
        elif (m==-2):
            res=(x*y)                     *(15./(4.*numpy.pi))**0.5
        else:
            res=None
    elif (n==3):
        if (m==0):
            res=(5.*z**3-3.*z)            *(7./(16.*numpy.pi))**0.5
        elif (m==1):
            res=x*(5.*z*z-1.)            *(21./(32.*numpy.pi))**0.5
        elif (m==-1):
            res=y*(5.*z*z-1.)            *(21./(32.*numpy.pi))**0.5
        elif (m==2):
            res=z*(x*x-y*y)             *(105./(16.*numpy.pi))**0.5
        elif (m==-2):
            res=2.*x*y*z                   *(105./(16.*numpy.pi))**0.5
        elif (m==3):
            res=(x*x*x-3.*x*y*y)           *(35./(32.*numpy.pi))**0.5
        elif (m==-3):
            res=(-y*y*y+3.*y*x*x)          *(35./(32.*numpy.pi))**0.5
        else:
            res=None
    elif (n==4):
        if (m==0):
            res=(35.*z*z*z*z-30.*z*z+3.)    *(9./(256.*numpy.pi))**0.5
        elif (m==1):
            res=x*(7.*z*z*z-3.*z)          *(45./(32.*numpy.pi))**0.5
        elif (m==-1):
            res=y*(7.*z*z*z-3.*z)          *(45./(32.*numpy.pi))**0.5
        elif (m==2):
            res=(x*x-y*y)*(7.*z*z-1.)  *(45./(64.*numpy.pi))**0.5
        elif (m==-2):
            res=2*x*y*(7.*z*z-1.)        *(45./(64.*numpy.pi))**0.5
        elif (m==3):
            res=z*(x*x*x-3.*x*y*y)         *(315./(32.*numpy.pi))**0.5
        elif (m==-3):
            res=z*(-y*y*y+3.*y*x*x)        *(315./(32.*numpy.pi))**0.5
        elif (m==4):
            res=(x*x*x*x-6.*x*x*y*y+y*y*y*y)  *(315./(256.*numpy.pi))**0.5
        elif (m==-4):
            res=(4.*x*x*x*y-4.*x*y*y*y)     *(315./(256.*numpy.pi))**0.5
        else:
            res=None
    else:
        res=None
    return res


def cartesian_to_spherical(x, y, z):
    x = numpy.asarray(x)
    y = numpy.asarray(y)
    z = numpy.asarray(z)

    r = numpy.sqrt(x**2 + y**2 + z**2)

    # Initialize theta and phi with zeros
    theta = numpy.zeros_like(r)
    phi = numpy.zeros_like(r)

    # Mask for non-zero radius
    mask = r > 0

    # Compute only where r > 0
    theta[mask] = numpy.arccos(numpy.clip(z[mask] / r[mask], -1.0, 1.0))
    phi[mask] = numpy.arctan2(y[mask], x[mask])

    return r, theta, phi



def calc_transformation_matrix_from_global_to_local(
        index_atom_origin_ax1:int, index_atom_end_ax1:int, label_ax1:str, 
        index_atom_origin_ax2:int, index_atom_end_ax2:int, label_ax2:str, 
        atom_fract_xyz:numpy.ndarray, matrix_m:numpy.ndarray):
    """To transform coordinates by this matrix the matrix should be inversed.
    matrix_m is matrix M. Dimension is [9,]. The order is [11,12,13,21,22,23,31,32,33]
    atom_fract_xyz is the position of atoms. Dimension is [3, N_atom] 
    index_atom_origin_ax1
    index_atom_end_ax1
    label_ax1
    index_atom_origin_ax2
    index_atom_end_ax2
    label_ax2
    Transfomation matrix gives coordinates of new vectors X',Y',Z' in notation of global coordinate sytem XYZ
    """    
    l_flag_ax1 = [hh in label_ax1.lower() for hh in ['x', 'y', 'z']]
    l_flag_ax2 = [hh in label_ax2.lower() for hh in ['x', 'y', 'z']]
    if not any(l_flag_ax1) or not any(l_flag_ax2):
        raise AttributeError("label_ax1 and label_ax2 should be 'X', 'Y' or 'Z'")

    axis_1_coord_o = atom_fract_xyz[:,index_atom_origin_ax1]
    axis_1_coord_f = atom_fract_xyz[:,index_atom_end_ax1]

    axis_2_coord_o = atom_fract_xyz[:,index_atom_origin_ax2]
    axis_2_coord_f = atom_fract_xyz[:,index_atom_end_ax2]

    axis_1_vector = calc_m_v(matrix_m, axis_1_coord_f-axis_1_coord_o, flag_m=False, flag_v=False)[0]
    axis_2_vector = calc_m_v(matrix_m, axis_2_coord_f-axis_2_coord_o, flag_m=False, flag_v=False)[0]
    
    axis_1_vector = axis_1_vector/numpy.sqrt(numpy.sum(numpy.square(axis_1_vector)))
    axis_2_vector = axis_2_vector - numpy.sum(axis_2_vector*axis_1_vector)*axis_1_vector
    axis_2_vector = axis_2_vector/numpy.sqrt(numpy.sum(numpy.square(axis_2_vector)))

    r_g_to_l = numpy.zeros(shape=(9,), dtype=float)
    flag_x, flag_y, flag_z = True, True, True
    if l_flag_ax1[0]:
        r_g_to_l[:3] = axis_1_vector
        flag_x = False
    elif l_flag_ax1[1]:
        r_g_to_l[3:6] = axis_1_vector
        flag_y = False
    else:# "z"
        r_g_to_l[6:9] = axis_1_vector
        flag_z = False

    if (l_flag_ax2[1] and flag_y):
        r_g_to_l[3:6] = axis_2_vector
        flag_y = False
    elif (l_flag_ax2[2] and flag_z):
        r_g_to_l[6:9] = axis_2_vector
        flag_z =False
    else:# "x"
        r_g_to_l[:3] = axis_2_vector
        flag_x =False

    if flag_x:
        r_g_to_l[:3] = calc_vector_product_v1_v2(r_g_to_l[3:6], r_g_to_l[6:9])[0]
    elif flag_y:
        r_g_to_l[3:6] = calc_vector_product_v1_v2(r_g_to_l[6:9], r_g_to_l[:3])[0]
    else:
        r_g_to_l[6:9] = calc_vector_product_v1_v2(r_g_to_l[:3], r_g_to_l[3:6])[0]
    return r_g_to_l



def calc_multipole_density(
        point_atom_symm_elems_auc, symm_elem_auc, unit_cell_parameters, 
        atom_rho_multipole_fract_xyz, atom_rho_multipole_transformation_matrix,
        atom_rho_multipole_lm, atom_rho_multipole_plm, l_n, l_zeta, l_coeff, l_kappa ):
    xyz_atoms_assigned_to_point = apply_symm_elem_to_xyz(point_atom_symm_elems_auc, atom_rho_multipole_fract_xyz) 
    r_d = point_atom_symm_elems_auc[4:13]
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters=unit_cell_parameters, flag_unit_cell_parameters=False)[0]
    r_xyz = calc_m1_m2_inv_m1(m_m, r_d)[0]
    # [9, n_points, n_atoms_multipole]
    np_r_g_to_l_s = calc_m1_m2_inv_m1(r_xyz, atom_rho_multipole_transformation_matrix)[0]

    diff_points_auc_xyz = xyz_atoms_assigned_to_point - numpy.expand_dims(symm_elem_auc[:3, :]/numpy.expand_dims(symm_elem_auc[3, :], axis=0), axis=2)
    flag_g_half = numpy.abs(diff_points_auc_xyz) > 0.5
    flag_sign = numpy.where(diff_points_auc_xyz >= 0.0, -1.0, 1.0)
    translation = numpy.where(flag_g_half, flag_sign, 0.)
    diff_points_auc_xyz =  diff_points_auc_xyz+translation

    diff_points_auc = calc_m_v(m_m,diff_points_auc_xyz, flag_m=False, flag_v=False)[0]
    # !!! probably it should be inversed matrix
    hh = np_r_g_to_l_s
    hh = calc_inv_m(np_r_g_to_l_s)[0]
    diff_points_auc_local_xyz = calc_m_v(hh, diff_points_auc, flag_m=False, flag_v=False)[0]
    r_sph, theta_sph, phi_sph = cartesian_to_spherical(diff_points_auc_local_xyz[0], diff_points_auc_local_xyz[1], diff_points_auc_local_xyz[2])

    i_mult = 0
    
    density_atom = numpy.zeros_like(r_sph)
    for n, zeta, coeff, kappa in zip(l_n, l_zeta, l_coeff, l_kappa):
        population = 1.
        rho_normalized = calc_rho_normalized(r_sph[:, i_mult], coeff, zeta, n, kappa=kappa)
        density_atom[:,i_mult] = 4.*numpy.pi*population*numpy.square(rho_normalized)
        i_mult += 1

    np_orient = numpy.zeros(shape=r_sph.shape, dtype=float)
    for lm, plm in zip(atom_rho_multipole_lm, atom_rho_multipole_plm):
        if numpy.any(numpy.logical_not(numpy.isclose(plm,0)),axis=0):
            np_orient += plm * realsphharm2df(int(lm[0]), int(lm[1]), theta_sph, phi_sph)
    # [N_point, N_multipole]
    multipole_density = density_atom * np_orient
    return multipole_density