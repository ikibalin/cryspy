"""
Functions:
    - gammanu_to_tthphi
    - tthphi_to_gammanu
    - recal_int_to_tthphi_grid
    - recal_int_to_gammanu_grid
    - app_grid
    - app_2d

"""
import math


def gammanu_to_tthphi(gamma, nu):
    """Transfer gamma-nu to ttheta-phi.

    gamma
        is the angle in the equatorial plane

    nu
        is the angle between equatorial plane and scattering neutron
        [-pi/2, pi/2]

    return
        ttheta is diffraction angle
        phi is polar angle
    """
    ttheta = math.acos(math.cos(gamma)*math.cos(nu))
    phi = math.atan2(math.sin(nu), math.cos(nu)*math.sin(gamma))
    return ttheta, phi


def tthphi_to_gammanu(tth, phi):
    """Transfer ttheta-phi to gamma-nu."""
    gamma = math.atan2(math.cos(phi)*math.sin(tth), math.cos(tth))
    nu = math.asin(math.sin(tth)*math.sin(phi))
    return gamma, nu


def recal_int_to_tthphi_grid(l_gamma_grid, l_nu_grid, ll_int_grid,
                             l_ttheta_grid, l_phi_grid):
    l_point = []
    for phi in l_phi_grid:
        for tth in l_ttheta_grid:
            point = tthphi_to_gammanu(tth, phi)
            l_point.append(point)

    l_int_dang_grid = app_grid(ll_int_grid, l_gamma_grid, l_nu_grid, l_point)

    n_tth = len(l_ttheta_grid)
    lint_out = [[l_int_dang_grid[ind_tth+ind_phi*n_tth]
                 for ind_tth, tth in enumerate(l_ttheta_grid)]
                for ind_phi, phi in enumerate(l_phi_grid)]
    return lint_out


def recal_int_to_gammanu_grid(ltth_grid, l_phi_grid, ll_int_grid,
                              l_gamma_grid, l_nu_grid):
    l_point = []
    for nu in l_nu_grid:
        for gamma in l_gamma_grid:
            point = gammanu_to_tthphi(gamma, nu)
            l_point.append(point)

    l_int_dang_grid = app_grid(ll_int_grid, ltth_grid, l_phi_grid, l_point)

    n_gamma = len(l_gamma_grid)
    lint_out = [[l_int_dang_grid[ind_gamma+ind_nu*n_gamma]
                 for ind_gamma, gamma in enumerate(l_gamma_grid)]
                for ind_nu, nu in enumerate(l_nu_grid)]
    return lint_out


def app_grid(mat_xy, x_grid, y_grid, l_point):
    l_res = []
    min_x, max_x = float(x_grid[0]), float(x_grid[-1])
    min_y, max_y = float(y_grid[0]), float(y_grid[-1])
    n_x, n_y = len(x_grid), len(y_grid)
    step_x = (max_x - min_x)/float(n_x-1)
    step_y = (max_y - min_y)/float(n_y-1)
    for ipoint, point in enumerate(l_point):
        try:
            val_x = point[0]
            val_y = point[1]
        except Exception:
            return []

        n_tot_x = (val_x - min_x)/step_x
        n_tot_y = (val_y - min_y)/step_y
        if all([(n_tot_x >= 0.), (n_tot_x < float(n_x-1)),
                (n_tot_y >= 0.), (n_tot_y < float(n_y-1))]):
            nx = n_tot_x % 1.
            ny = n_tot_y % 1.
            ind_x = int(n_tot_x//1.)
            ind_y = int(n_tot_y//1.)
            res = app_2d(mat_xy[ind_y][ind_x], mat_xy[ind_y][ind_x+1],
                         mat_xy[ind_y+1][ind_x], mat_xy[ind_y+1][ind_x+1],
                         nx, ny)
        else:
            res = None
        l_res.append(res)
    return l_res


def app_2d(f11, f12, f21, f22, nx, ny):
    try:
        res2 = f21 - nx*(f21-f22)
        res1 = f11 - nx*(f11-f12)
        res = res1 - ny*(res1-res2)
    except Exception:
        res = None
    return res
