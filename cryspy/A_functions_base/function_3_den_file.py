"""Module define operations with density file: file.den

Functions
---------
    - read_den_file
    - save_to_den_file
    - calc_density_grid_by_den_data
    - read_den_file_grid

"""
import warnings

import numpy

from cryspy.A_functions_base.function_2_sym_elems import \
    form_symm_elems_by_b_i_r_ij, calc_numerators_denominator_for_b_i


def read_den_file(file_name: str):  # noqa: PLR0915
    """Read density file.
    
    Arguments
    ---------
        -   file_name is a file name str

    points_abc: [3] int
    cell_parameters: [6] float
    sym_elems: [13, n_symm_elems]: int
    indexes_xyz: [3, n_points] int
    densities: [3, n_points] float
    centrosymmetry: bool
    np_shift: [4, n_symmelems]  int
    """
    with open(file_name, "r") as fid:
        l_content = fid.readlines()

    number_lines = int(l_content[1])

    hh = l_content[number_lines+2].strip().split()
    rad = float(numpy.pi/180.)
    cell_parameters = numpy.array([float(hh[0]), float(hh[1]), float(hh[2]), 
                       float(hh[3])*rad, float(hh[4])*rad, float(hh[5])*rad],
                                  dtype=float)

    [points_a, points_b, points_c, n_el_symm, centr, n_shift] = [
        int(hh) for hh in l_content[number_lines+3][:-1].split()]
    points_abc = numpy.array((points_a, points_b, points_c), dtype=int)
    centrosymmetry = bool(centr)

    l_ind_xyz = []
    l_dens = []
    for line in l_content[2:number_lines+2]:
        hh = line.strip().split()
        ind_x, ind_y, ind_z = int(hh[0]), int(hh[1]), int(hh[2])
        den_f, den_a = 0., 0.
        den = float(hh[3])
        if len(hh) == 6:
            den_f = float(hh[4])
            den_a = float(hh[5])
        elif den >= 0.:
            den_f = den
        else:  # den < 0 
            den_a = den
        l_ind_xyz.append((ind_x, ind_y, ind_z))
        l_dens.append((den, den_f, den_a))
    indexes_xyz = numpy.array(l_ind_xyz, dtype=int).transpose()
    densities = numpy.array(l_dens, dtype=float).transpose()

    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = [], [], [], [], \
        [], [], [], [], []
    b_1, b_2, b_3 = [], [], []
    for line in l_content[number_lines+4:number_lines+4+n_el_symm]:
        hh = line.replace("-", " -").strip().split()
        r_11.append(int(hh[0]))
        r_12.append(int(hh[3]))
        r_13.append(int(hh[6]))
        
        r_21.append(int(hh[1]))
        r_22.append(int(hh[4]))
        r_23.append(int(hh[7]))

        r_31.append(int(hh[2]))
        r_32.append(int(hh[5]))
        r_33.append(int(hh[8]))

        b_1.append(float(hh[9]))
        b_2.append(float(hh[10]))
        b_3.append(float(hh[11]))

    b_i = (b_1, b_2, b_3)
    r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

    sym_elems = form_symm_elems_by_b_i_r_ij(b_i, r_ij)

    shift_1, shift_2, shift_3 = [], [], []
    for line in l_content[number_lines+4+n_el_symm:
                          number_lines+4+n_el_symm+n_shift]:
        hh = line.strip().split()
        shift_1.append(float(hh[0]))
        shift_2.append(float(hh[1]))
        shift_3.append(float(hh[2]))

    sh_num_x, sh_num_y, sh_num_z, sh_den = calc_numerators_denominator_for_b_i(
        shift_1, shift_2, shift_3)

    np_shift = numpy.stack((sh_num_x, sh_num_y, sh_num_z, sh_den), axis=0)

    return points_abc, cell_parameters, sym_elems, indexes_xyz, densities, \
        centrosymmetry, np_shift


def save_to_den_file(
        file_name: str, points_abc, cell_parameters, sym_elems, indexes_xyz,
        densities, centrosymmetry, shift):
    """Save to file.
    file_name is str
    points_abc: [3] int
    cell_parameters: [6] float
    sym_elems: [13, n_symm_elems]: int
    indexes_xyz: [3, n_points] int
    densities: [3, n_points] float
    centrosymmetry: bool
    np_shift: [4, n_symmelems]  int
    """

    index_x = indexes_xyz[0, :]
    index_y = indexes_xyz[1, :]
    index_z = indexes_xyz[2, :]

    n_symmetry = sym_elems.shape[1]

    n_shape = densities.shape
    if len(n_shape) == 1:
        n_points = n_shape[0]
        n_size = 1 
    else:
        n_size, n_points = n_shape
    if n_size == 3:
        density = densities[0, :]
        density_f = densities[1, :]
        density_a = densities[2, :]
    elif n_size == 2:
        density_f = densities[0, :]
        density_a = densities[1, :]
        density = density_f-density_a
    else:
        density = densities
        density_f = numpy.zeros(density.shape, dtype=float)
        density_a = numpy.zeros(density.shape, dtype=float)
        cond = density >= 0.
        density_f[cond] = density[cond]
        cond_not = numpy.logical_not(cond)
        density_a[cond_not] = density[cond_not]

    ls_out = []
    ls_out.append("Created by CrysPy")
    ls_out.append("{:}".format(n_points))

    for i_x, i_y, i_z, den, den_f, den_a in zip(
            index_x, index_y, index_z, density, density_f, density_a):

        ls_out.append(
            f"{i_x:4} {i_y:4} {i_z:4} {den:15.7f} {den_f:15.7f} {den_a:15.7f}")

    (a, b, c, al, be, ga) = cell_parameters

    ls_out.append(
        f"{a:10.5f}{b:10.5f}{c:10.5f}{al:10.5f}{be:10.5f}{ga:10.5f}")

    n_shift = shift.shape[1]

    ls_out.append(f"{points_abc[0]:5}{points_abc[1]:5}{points_abc[2]:5}\
{n_symmetry:5}{centrosymmetry:5}{n_shift:5}")

    for r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_num_1,  \
            b_num_2, b_num_3, b_den in zip(
            sym_elems[4, :], sym_elems[5, :], sym_elems[6, :], 
            sym_elems[7, :], sym_elems[8, :], sym_elems[9, :], 
            sym_elems[10, :], sym_elems[11, :], sym_elems[12, :], 
            sym_elems[0, :], sym_elems[1, :], sym_elems[2, :],
            sym_elems[3, :]):
        b_1 = float(b_num_1)/float(b_den)
        b_2 = float(b_num_2)/float(b_den)
        b_3 = float(b_num_3)/float(b_den)
        ls_out.append(
            f"{r_11:4}{r_21:4}{r_31:4}  {r_12:4}{r_22:4}{r_32:4}  \
{r_13:4}{r_23:4}{r_33:4}    {b_1:8.5f}{b_2:8.5f}{b_3:8.5f}")
    for orig_1, orig_2, orig_3 in zip(shift[0, :], shift[1, :], shift[2, :]):
        ls_out.append(f"{orig_1:8.4f}{orig_2:8.4f}{orig_3:8.4f}")

    with open(file_name, "w") as fid:
        fid.write("\n".join(ls_out))


def calc_density_grid_by_den_data(
        points_abc, sym_elems, indexes_xyz, density, centrosymmetry: bool,
        shift):
    """Density on the whole grid of the unit cell from the points of a .den file.

    A .den file keeps the density in the points of the asymmetric unit;
    the other points are their images under the symmetry elements, the
    centring translations and, for a centrosymmetric group, the inversion
    (at the origin: the .den format does not keep its position).

    Parameters
    ----------
    points_abc : numpy.ndarray, shape (3,), int
        Number of grid points along a, b, c; point (i, j, k) is at the
        fractional coordinates (i/n_a, j/n_b, k/n_c).
    sym_elems : numpy.ndarray, shape (13, n_symm), int
        Symmetry elements: b_1, b_2, b_3 numerators, their denominator,
        r_11, r_12, ..., r_33.
    indexes_xyz : numpy.ndarray, shape (3, n_point), int
        Grid indexes of the given points.
    density : numpy.ndarray, shape (n_point,)
        Density in the given points (μB/Å³ for MEM magnetization).
    centrosymmetry : bool
        Add the inversion -x, -y, -z to every element.
    shift : numpy.ndarray, shape (4, n_shift), int
        Centring translations: numerators and the denominator.

    Returns
    -------
    numpy.ndarray, shape (n_a, n_b, n_c)
        The density; NaN in points that no element reaches. A warning is
        given when images of the points do not fall on the grid or give
        one point different densities (the file and its symmetry disagree).
    """
    n_abc = numpy.asarray(points_abc, dtype=int)
    sym_elems = numpy.asarray(sym_elems, dtype=int)
    shift = numpy.asarray(shift, dtype=int)
    fract = numpy.asarray(indexes_xyz, dtype=float) / n_abc[:, numpy.newaxis]
    r_ij = sym_elems[4:13].reshape(3, 3, -1).astype(float)
    b_i = sym_elems[0:3].astype(float) / sym_elems[3][numpy.newaxis, :]
    t_i = shift[0:3].astype(float) / shift[3][numpy.newaxis, :]
    signs = (1., -1.) if centrosymmetry else (1.,)

    # [3, n_point, n_symm, n_shift, n_sign]
    rotated = numpy.einsum("ijs,jp->ips", r_ij, fract)
    images = (numpy.array(signs)[numpy.newaxis, numpy.newaxis, numpy.newaxis, numpy.newaxis, :] *
              rotated[:, :, :, numpy.newaxis, numpy.newaxis] +
              b_i[:, numpy.newaxis, :, numpy.newaxis, numpy.newaxis] +
              t_i[:, numpy.newaxis, numpy.newaxis, :, numpy.newaxis])
    index_float = images * n_abc.reshape(3, 1, 1, 1, 1)
    index_int = numpy.rint(index_float)
    if not numpy.allclose(index_float, index_int, atol=1e-4):
        warnings.warn(
            "Symmetry images of the density points do not fall on the grid "
            f"{tuple(n_abc)}: they are taken at the nearest points.",
            stacklevel=2)
    index_int = numpy.mod(index_int.astype(int), n_abc.reshape(3, 1, 1, 1, 1))
    values = numpy.broadcast_to(
        numpy.asarray(density, dtype=float).reshape(-1, 1, 1, 1),
        index_int.shape[1:])

    flat = numpy.ravel_multi_index(
        tuple(index_int.reshape(3, -1)), tuple(n_abc))
    values = values.reshape(-1)
    n_total = int(numpy.prod(n_abc))
    counts = numpy.bincount(flat, minlength=n_total)
    sums = numpy.bincount(flat, weights=values, minlength=n_total)
    with numpy.errstate(invalid="ignore", divide="ignore"):
        grid = sums / counts
    grid[counts == 0] = numpy.nan
    spread = numpy.abs(values - grid[flat])
    scale = numpy.nanmax(numpy.abs(grid)) if numpy.any(counts) else 0.
    if numpy.any(spread > 1e-3 * scale + 1e-9):
        warnings.warn(
            "Symmetry images of the density points give one grid point "
            "different densities: the symmetry of the .den file does not "
            "fit its points (the mean is taken).", stacklevel=2)
    return grid.reshape(tuple(n_abc))


def read_den_file_grid(file_name: str):
    """Read a .den file and give its density on the whole grid of the cell.

    Parameters
    ----------
    file_name : str
        Name of the .den file (e.g. written by MEM in CrysPy).

    Returns
    -------
    cell_parameters : numpy.ndarray, shape (6,)
        a, b, c (Å), alpha, beta, gamma (radians).
    density_grid : numpy.ndarray, shape (n_a, n_b, n_c)
        Total density (the first density column of the file) in the points
        (i/n_a, j/n_b, k/n_c); see :func:`calc_density_grid_by_den_data`.
    """
    points_abc, cell_parameters, sym_elems, indexes_xyz, densities, \
        centrosymmetry, np_shift = read_den_file(file_name)
    density_grid = calc_density_grid_by_den_data(
        points_abc, sym_elems, indexes_xyz, densities[0], centrosymmetry,
        np_shift)
    return cell_parameters, density_grid
