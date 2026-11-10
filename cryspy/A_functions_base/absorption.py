"""Absorption of neutrons in a single crystal.

The crystal is a sphere or an ellipsoid. Its transmission factor for one
reflection,

    A = 1/V int_V exp(-mu (t_in + t_out)) dV,

is the mean over the volume of the attenuation along the path t_in of the
incident beam to a point of the crystal and the path t_out of the diffracted
beam from it (single scattering). It is integrated numerically with a product
Gauss rule in the unit ball, the "gaussian" correction of CIF
(_exptl_absorpt_correction_type; Busing & Levy, Acta Cryst. 10 (1957) 180;
Coppens, in Crystallographic Computing (1970) 255).

Coordinate systems. The directions of the beams are unit vectors in the
Cartesian system of the crystal (x || a*, z || c). The laboratory system of
the orientation matrix U (CCSL notation, matrix_u of the dictionaries) has x
along the incident beam and z vertical: v_lab = U v_ccs at zero angles of the
goniometer. The scattering vector is Q = k_f - k_i, |k| = 1/lambda, in 1/Ang
without 2 pi.

Linear absorption coefficient. mu = n (sigma_abs + sigma_scat) for the
number density n of the atoms; sigma_abs is tabulated for 2200 m/s
(1.798 Ang) and scales as 1/v (proportional to the wavelength), the bound
scattering cross section sigma_scat does not depend on the wavelength
(Sears, Neutron News 3 (1992) 26). mu is calculated from the composition of the
measured crystal; the 1/v law is wrong for the nuclides with resonances at
thermal energies (Cd, Sm, Eu, Gd, Dy, ...): give mu there.

Units: mu in 1/mm, sizes of the crystal in mm, wavelength in Ang,
angles in degrees.
"""
import re
import warnings
import numpy

from cryspy.A_functions_base.database import DATABASE

na = numpy.newaxis

WAVELENGTH_TABULATED = 1.798  # Ang, of the tabulated absorption cross sections

# nuclides whose absorption deviates strongly from the 1/v law for thermal
# and cold neutrons (resonances below ~1 eV)
RESONANT_ELEMENTS = ("Cd", "In", "Sm", "Eu", "Gd", "Dy", "Er", "Hf", "Ir", "Lu")

KEY_ABSORPTION = "absorption cross section for 2200 m/s neutrons (barn)"
KEY_SCATTERING = "total bound scattering cross section (barn)"

RE_NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?|[-+]?\d+\.")


def take_element_symbol(type_symbol: str) -> str:
    """Element or isotope of a type symbol of atom_site: 'Fe3+' -> 'Fe',
    'O2-' -> 'O', '157Gd' -> '157Gd', 'D' -> 'D'."""
    s_symbol = type_symbol.strip()
    match = re.match(r"^(\d*[A-Z][a-z]?)", s_symbol)
    if match is None:
        raise ValueError(f"The type symbol '{type_symbol:}' is not an element.")
    return match.group(1)


def _take_cross_section(key: str, symbol: str) -> float:
    value = DATABASE["Isotopes"].get((key, symbol))
    if value is None:
        raise ValueError(
            f"No neutron cross sections for '{symbol:}' in the database of CrysPy: "
            "give the absorption coefficient _exptl_absorpt_coefficient_mu.")
    if isinstance(value, str):
        # values with the uncertainty, e.g. '49700.(125.)'
        match = RE_NUMBER.search(value.split("(")[0] if not value.startswith("(") else value)
        if match is None:
            raise ValueError(f"Cross section '{value:}' of '{symbol:}' is not a number.")
        value = float(match.group(0))
    value = float(numpy.real(value))
    if numpy.isnan(value):
        raise ValueError(f"No neutron cross section '{key:}' for '{symbol:}'.")
    return value


def take_neutron_cross_sections(type_symbol: str):
    """Absorption (at 1.798 Ang) and total bound scattering cross sections
    (barn per atom) of the element or isotope of the type symbol."""
    symbol = take_element_symbol(type_symbol)
    if symbol == "D":
        symbol = "2H"
    sigma_abs = _take_cross_section(KEY_ABSORPTION, symbol)
    sigma_scat = _take_cross_section(KEY_SCATTERING, symbol)
    return sigma_abs, sigma_scat


def calc_absorption_coefficients_by_composition(
        type_symbols, numbers_in_cell, volume_cell: float, wavelength: float):
    """Linear absorption coefficients of the crystal.

    Parameters
    ----------
    type_symbols : list of str
        Type symbols of the atoms ('Fe3+', 'O').
    numbers_in_cell : array (n_atom,)
        Number of each atom in the unit cell (multiplicity x occupancy).
    volume_cell : float
        Volume of the unit cell, Ang^3.
    wavelength : float
        Wavelength, Ang.

    Returns
    -------
    mu_abs, mu_scat : float
        True absorption and scattering parts of mu, 1/mm; mu = mu_abs + mu_scat.
    """
    numbers_in_cell = numpy.asarray(numbers_in_cell, dtype=float)
    sigma = numpy.array([take_neutron_cross_sections(symbol) for symbol in type_symbols],
                        dtype=float).reshape(-1, 2)
    # barn / Ang^3 = 1e-24 cm^2 / 1e-24 cm^3 = 1/cm; 1/cm = 0.1/mm
    n_sigma_abs = float((numbers_in_cell * sigma[:, 0]).sum()) / volume_cell * 0.1
    n_sigma_scat = float((numbers_in_cell * sigma[:, 1]).sum()) / volume_cell * 0.1
    mu_abs = n_sigma_abs * wavelength / WAVELENGTH_TABULATED
    return mu_abs, n_sigma_scat


def calc_absorption_coefficients_by_dictionary(dict_crystal: dict, wavelength: float):
    """mu of the crystal at the wavelength and at half of it from the
    composition of a crystal dictionary (Crystal.get_dictionary).

    Parameters
    ----------
    dict_crystal : dict
        Keys 'atom_type_symbol', 'atom_multiplicity', 'atom_occupancy' and
        'unit_cell_parameters' (a, b, c in Ang, angles in radians).
    wavelength : float
        Wavelength, Ang.

    Returns
    -------
    None or (mu, mu_half)
        1/mm; None when the dictionary has no atoms. At lambda/2 the true
        absorption halves (1/v law), the scattering part stays.
    """
    keys = ("atom_type_symbol", "atom_multiplicity", "atom_occupancy", "unit_cell_parameters")
    if dict_crystal is None or not all(key in dict_crystal for key in keys):
        return None
    from cryspy.A_functions_base.unit_cell import calc_volume_uc_by_unit_cell_parameters
    volume = float(calc_volume_uc_by_unit_cell_parameters(dict_crystal["unit_cell_parameters"])[0])
    numbers = numpy.asarray(dict_crystal["atom_multiplicity"], dtype=float) * \
        numpy.asarray(dict_crystal["atom_occupancy"], dtype=float)
    mu_abs, mu_scat = calc_absorption_coefficients_by_composition(
        [str(symbol) for symbol in dict_crystal["atom_type_symbol"]], numbers, volume, wavelength)
    return mu_abs + mu_scat, 0.5 * mu_abs + mu_scat


def take_resonant_elements(type_symbols) -> list:
    """Elements of the type symbols that do not follow the 1/v law."""
    l_out = []
    for symbol in type_symbols:
        element = take_element_symbol(symbol).lstrip("0123456789")
        if element in RESONANT_ELEMENTS and element not in l_out:
            l_out.append(element)
    return l_out


def calc_beam_directions_normal_beam(q_ccs, axis_z, wavelength, diffracted_beam: str = "left"):
    """Directions of the incident and diffracted beams of the normal beam
    (lifting counter) geometry.

    The crystal turns about the vertical axis, the incident beam is
    horizontal and the detector can be lifted: this is the geometry of the
    polarized neutron diffractometers (D3, 6T2, 5C1, POLI), where the vertical
    axis is along the magnetic field.

    Parameters
    ----------
    q_ccs : array (3, n)
        Scattering vectors in the Cartesian system of the crystal, 1/Ang
        without 2 pi.
    axis_z : array (3,)
        Vertical axis in the same system (the third row of U, CCSL).
    wavelength : float or array (n,)
        Wavelength, Ang.
    diffracted_beam : str
        'left' or 'right': side of the detector (as calc_angle of
        DiffrnOrientMatrix).

    Returns
    -------
    e_in, e_out : arrays (3, n)
        Unit vectors along the incident and diffracted beams.
    flag_reachable : array (n,) of bool
        False for the reflections out of reach (|Q_horizontal| too small).
    """
    q_ccs = numpy.asarray(q_ccs, dtype=float)
    axis = numpy.asarray(axis_z, dtype=float)
    axis = axis / numpy.linalg.norm(axis)
    wavelength = numpy.asarray(wavelength, dtype=float)
    q_sq = numpy.square(q_ccs).sum(axis=0)
    q_par = (q_ccs * axis[:, na]).sum(axis=0)
    q_perp = q_ccs - q_par[na, :] * axis[:, na]
    q_perp_norm = numpy.sqrt(numpy.square(q_perp).sum(axis=0))
    with numpy.errstate(divide="ignore", invalid="ignore"):
        e_1 = q_perp / q_perp_norm[na, :]
        # e_in . Q = -lambda Q^2 / 2 from |k_i + Q| = |k_i|
        cos_alpha = -wavelength * q_sq / (2. * q_perp_norm)
    flag_reachable = numpy.abs(cos_alpha) <= 1.
    cos_alpha = numpy.clip(numpy.nan_to_num(cos_alpha, nan=-1.), -1., 1.)
    sin_alpha = numpy.sqrt(1. - numpy.square(cos_alpha))
    # left: (e_in x Q) . axis > 0, i.e. sin_alpha < 0 (see the module docs)
    if diffracted_beam.strip().lower() == "left":
        sin_alpha = -sin_alpha
    e_1 = numpy.nan_to_num(e_1)
    e_2 = numpy.cross(axis[:, na], e_1, axis=0)
    e_in = cos_alpha[na, :] * e_1 + sin_alpha[na, :] * e_2
    e_out = e_in + wavelength * q_ccs
    e_out = e_out / numpy.sqrt(numpy.square(e_out).sum(axis=0))[na, :]
    return e_in, e_out, flag_reachable


def calc_rotation_matrix_busing_levy(omega, chi, phi):
    """R = Omega Chi Phi of the four-circle diffractometer (Busing & Levy,
    Acta Cryst. 22 (1967) 457), angles in degrees, arrays (n,); returns
    (3, 3, n). Omega and phi turn about z, chi about y."""
    omega, chi, phi = (numpy.radians(numpy.asarray(hh, dtype=float)) for hh in (omega, chi, phi))
    zeros, ones = numpy.zeros_like(omega), numpy.ones_like(omega)

    def rot_z(angle):
        return numpy.array([[numpy.cos(angle), numpy.sin(angle), zeros],
                            [-numpy.sin(angle), numpy.cos(angle), zeros],
                            [zeros, zeros, ones]])

    m_chi = numpy.array([[numpy.cos(chi), zeros, numpy.sin(chi)],
                         [zeros, ones, zeros],
                         [-numpy.sin(chi), zeros, numpy.cos(chi)]])
    return numpy.einsum("ijn,jkn,kln->iln", rot_z(omega), m_chi, rot_z(phi))


def calc_beam_directions_by_angles(q_ccs, matrix_u, wavelength, omega, chi, phi):
    """Directions of the beams of a four-circle measurement.

    The rotation matrices of the circles are those of Busing & Levy (1967),
    R = Omega Chi Phi, but omega is the angle of the sample itself (as the
    motors of the ILL four-circle diffractometers and PyColl5 give it), not
    omega - theta of the bisecting convention: Q = R^-1 (k_f - k_i) in the phi
    frame. The incident beam is along +y of the laboratory system of Busing &
    Levy, which is x of the CCSL system of U (x_ccsl = y_BL, y_ccsl = -x_BL,
    z = z, as the type 'BusingLevy' of DiffrnOrientMatrix).

    Parameters
    ----------
    q_ccs : array (3, n)
        Scattering vectors in the Cartesian system of the crystal, 1/Ang.
    matrix_u : array (9,)
        Orientation matrix U (CCSL), row by row.
    wavelength : float or array (n,)
        Wavelength, Ang.
    omega, chi, phi : arrays (n,)
        Angles of the goniometer at the peak, degrees.

    Returns
    -------
    e_in, e_out : arrays (3, n)
        Unit vectors along the incident and diffracted beams in the
        Cartesian system of the crystal.
    misfit : array (n,)
        | |e_in + lambda Q| - 1 |: zero when the angles bring the reflection
        exactly into the diffraction condition.
    """
    q_ccs = numpy.asarray(q_ccs, dtype=float)
    m_u = numpy.asarray(matrix_u, dtype=float).reshape(3, 3)
    m_r = calc_rotation_matrix_busing_levy(omega, chi, phi)
    # incident beam in the phi frame of Busing & Levy: R^-1 (0, 1, 0) = R^T (0, 1, 0)
    e_in_bl = m_r[1, :, :]
    e_in_ccsl = numpy.stack([e_in_bl[1], -e_in_bl[0], e_in_bl[2]], axis=0)
    e_in = m_u.T @ e_in_ccsl
    k_f = e_in + numpy.asarray(wavelength, dtype=float) * q_ccs
    norm = numpy.sqrt(numpy.square(k_f).sum(axis=0))
    return e_in, k_f / norm[na, :], numpy.abs(norm - 1.)


def form_unit_ball_rule(n_r: int = 8, n_theta: int = 10, n_phi: int = 16):
    """Points (3, p) and weights (p,) of a product Gauss rule in the unit
    ball; the weights sum to 1 (a mean over the volume)."""
    x_r, w_r = numpy.polynomial.legendre.leggauss(n_r)
    radius = 0.5 * (x_r + 1.)
    w_r = 0.5 * w_r * numpy.square(radius) * 3.  # int_0^1 3 r^2 dr = 1
    cos_t, w_t = numpy.polynomial.legendre.leggauss(n_theta)
    w_t = 0.5 * w_t
    angle_phi = (numpy.arange(n_phi) + 0.5) * 2. * numpy.pi / n_phi
    w_p = numpy.full(n_phi, 1. / n_phi)
    rr, ct, pp = numpy.meshgrid(radius, cos_t, angle_phi, indexing="ij")
    st = numpy.sqrt(1. - numpy.square(ct))
    points = numpy.stack([(rr * st * numpy.cos(pp)).ravel(), (rr * st * numpy.sin(pp)).ravel(),
                          (rr * ct).ravel()], axis=0)
    weights = (w_r[:, na, na] * w_t[na, :, na] * w_p[na, na, :]).ravel()
    return points, weights


# 16 x 16 x 24 points: relative error of the transmission below 1e-3 for
# mu R <= 2 and about 1e-2 for mu R = 4 (sphere at 2theta = 0, the worst case)
UNIT_BALL_RULE = form_unit_ball_rule(16, 16, 24)
# used for mu R > 2: relative error about 3e-3 for mu R = 4
UNIT_BALL_RULE_FINE = form_unit_ball_rule(32, 20, 32)
MU_R_ACCURATE = 4.  # mu R above which a warning is given
MISFIT_TOLERANCE = 0.05  # | |k_f| lambda - 1 | of the four-circle angles


def calc_chord_in_unit_ball(points, vectors):
    """Parameter t >= 0 (n, p) at which points + t vectors leave the unit
    ball; points (3, p) inside it, vectors (3, n)."""
    se = vectors.T @ points
    ee = numpy.maximum(numpy.square(vectors).sum(axis=0)[:, na], 1e-30)
    ss = numpy.square(points).sum(axis=0)[na, :]
    return (-se + numpy.sqrt(numpy.maximum(se * se + ee * (1. - ss), 0.))) / ee


def calc_rotation_matrix_by_euler_zyz(euler_angles):
    """Rotation matrix (3, 3) of the Euler angles alpha, beta, gamma
    (degrees, z-y-z convention, R = Rz(alpha) Ry(beta) Rz(gamma)); its
    columns are the axes of the rotated frame."""
    alpha, beta, gamma = numpy.radians(numpy.asarray(euler_angles, dtype=float))

    def rot_z(angle):
        return numpy.array([[numpy.cos(angle), -numpy.sin(angle), 0.],
                            [numpy.sin(angle), numpy.cos(angle), 0.], [0., 0., 1.]])

    rot_y = numpy.array([[numpy.cos(beta), 0., numpy.sin(beta)], [0., 1., 0.],
                         [-numpy.sin(beta), 0., numpy.cos(beta)]])
    return rot_z(alpha) @ rot_y @ rot_z(gamma)


def calc_euler_zyz_by_rotation_matrix(matrix):
    """Euler angles alpha, beta, gamma (degrees, z-y-z) of a proper rotation
    matrix (3, 3); inverse of calc_rotation_matrix_by_euler_zyz."""
    m = numpy.asarray(matrix, dtype=float)
    beta = numpy.arccos(numpy.clip(m[2, 2], -1., 1.))
    if numpy.sin(beta) > 1e-8:
        alpha = numpy.arctan2(m[1, 2], m[0, 2])
        gamma = numpy.arctan2(m[2, 1], -m[2, 0])
    else:  # gimbal lock: only alpha + gamma (or alpha - gamma) is defined
        alpha = numpy.arctan2(m[1, 0], m[0, 0]) if m[2, 2] > 0 else numpy.arctan2(-m[1, 0], -m[0, 0])
        gamma = 0.
    return numpy.degrees(numpy.array([alpha, beta, gamma]))


def calc_transmission_ellipsoid(mu, semi_axes, matrix_axes, e_in, e_out, rule=None):
    """Transmission factor A of each reflection for an ellipsoidal crystal.

    Parameters
    ----------
    mu : float or array (n,)
        Linear absorption coefficient, 1/mm.
    semi_axes : array (3,)
        Semi-axes of the ellipsoid, mm (equal for a sphere).
    matrix_axes : array (3, 3)
        Directions of the semi-axes as columns, Cartesian system of the
        crystal.
    e_in, e_out : arrays (3, n)
        Unit vectors along the incident and diffracted beams.
    rule : tuple, optional
        Points and weights in the unit ball (form_unit_ball_rule).

    Returns
    -------
    array (n,)
        A = <exp(-mu (t_in + t_out))> over the volume, 0 < A <= 1.
    """
    points, weights = UNIT_BALL_RULE if rule is None else rule
    semi_axes = numpy.asarray(semi_axes, dtype=float)
    # x = M D p maps the unit ball onto the ellipsoid; a path t along e is
    # the path t' = t along D^-1 M^T e in the unit ball
    inv_l = (numpy.asarray(matrix_axes, dtype=float) / semi_axes[na, :]).T
    v_in = -(inv_l @ numpy.atleast_2d(e_in))  # back to where the beam came from
    v_out = inv_l @ numpy.atleast_2d(e_out)
    path = calc_chord_in_unit_ball(points, v_in) + calc_chord_in_unit_ball(points, v_out)
    mu = numpy.asarray(mu, dtype=float)
    if mu.ndim:
        mu = mu[:, na]
    return numpy.exp(-mu * path) @ weights


def calc_transmission_sphere_zero_angle(mu_radius):
    """Transmission factor of a sphere at 2theta = 0 for mu R:
    A = 3 / x^3 [2 - exp(-x) (x^2 + 2x + 2)] with x = mu D = 2 mu R (path in
    and out averaged over the volume; International Tables C, 6.3.3)."""
    x = 2. * numpy.asarray(mu_radius, dtype=float)  # mu D
    xs = numpy.maximum(x, 1e-2)
    exact = 3. / xs**3 * (2. - numpy.exp(-xs) * (xs**2 + 2. * xs + 2.))
    return numpy.where(x < 1e-2, 1. - 3. * x / 4. + 3. * x**2 / 10. - x**3 / 12., exact)


def calc_transmission_by_dictionary(
        dict_diffrn: dict, index_hkl, unit_cell_parameters, wavelength,
        flag_half_wavelength: bool = False, d_parameter: dict = None,
        flag_warn: bool = True, dict_crystal: dict = None):
    """Transmission factor of the measured reflections of a diffrn
    dictionary, or None when absorption is not corrected.

    mu is calculated from the composition of the measured crystal
    (dict_crystal, the phase of the experiment); 'absorption_mu' of the
    experiment, when given, replaces it (resonant absorbers). mu at lambda/2
    follows by the ratio of the composition.

    Parameters
    ----------
    dict_diffrn : dict
        Dictionary of the experiment with the keys of ExptlAbsorpt and
        ExptlCrystal ('absorption_correction_type', 'absorption_mu', ...),
        'matrix_u' and, optionally, 'angle_omega_chi_phi' (3, n).
    index_hkl : array (3, n)
        Measured reflections.
    unit_cell_parameters : array (6,)
        a, b, c (Ang), alpha, beta, gamma (radians).
    wavelength : float or array (n,)
        Wavelength, Ang.
    flag_half_wavelength : bool
        Also give the transmission of the lambda/2 contamination (the
        reflections 2h at lambda/2 follow the same paths).
    d_parameter : dict, optional
        Values replacing those of dict_diffrn (for numerical derivatives).
    flag_warn : bool
        Warn about doubtful input (False for repeated calls).
    dict_crystal : dict, optional
        Dictionary of the measured crystal, for mu from its composition.

    Returns
    -------
    None or (transmission, transmission_half)
        Arrays (n,); transmission_half is None without flag_half_wavelength.
    """
    correction_type = dict_diffrn.get("absorption_correction_type", "none")
    if correction_type not in ("sphere", "gaussian"):
        return None
    d_value = dict(dict_diffrn)
    if d_parameter is not None:
        d_value.update(d_parameter)
    wavelength_mean = float(numpy.mean(wavelength))
    try:
        mu_composition = calc_absorption_coefficients_by_dictionary(dict_crystal, wavelength_mean)
    except ValueError as error:
        if "absorption_mu" not in d_value:
            raise ValueError(f"Absorption: mu is not calculated from the crystal ({error:}).") from error
        mu_composition = None
    if "absorption_mu" in d_value:
        # the minimizer of the refinement may try negative values: the
        # transmission depends on |mu| and on the lengths of the axes only
        mu = abs(float(numpy.atleast_1d(d_value["absorption_mu"])[0]))
    elif mu_composition is not None:
        mu = mu_composition[0]
        if flag_warn:
            l_resonant = take_resonant_elements(dict_crystal["atom_type_symbol"])
            if l_resonant:
                warnings.warn(
                    "Absorption: mu from the composition follows the 1/v law, which is wrong for " +
                    ", ".join(l_resonant) + "; give _exptl_absorpt_coefficient_mu at this wavelength.",
                    UserWarning, stacklevel=2)
    else:
        if flag_warn:
            warnings.warn(
                f"Absorption correction '{correction_type:}' is not applied: no crystal structure "
                "of the phase to calculate mu (or give _exptl_absorpt_coefficient_mu).",
                UserWarning, stacklevel=2)
        return None
    ratio_half = 1.
    if mu_composition is not None and mu_composition[0] > 0.:
        ratio_half = mu_composition[1] / mu_composition[0]
    if correction_type == "sphere":
        if "absorption_radius" not in d_value:
            raise ValueError("Absorption correction 'sphere': give the radius of the crystal "
                             "_exptl_crystal_size_rad (mm).")
        radius = abs(float(numpy.atleast_1d(d_value["absorption_radius"])[0]))
        semi_axes = numpy.array([radius, radius, radius], dtype=float)
        matrix_axes = numpy.eye(3)
        l_unused = ["absorption_size", "absorption_euler"]
    else:
        if "absorption_size" not in d_value:
            raise ValueError("Absorption correction 'gaussian': give the sizes of the crystal "
                             "_exptl_crystal_size_max, _size_mid, _size_min (mm).")
        semi_axes = 0.5 * numpy.abs(numpy.asarray(d_value["absorption_size"], dtype=float))
        matrix_axes = calc_rotation_matrix_by_euler_zyz(d_value["absorption_euler"])
        l_unused = ["absorption_radius"]
    if numpy.any(semi_axes == 0.):
        raise ValueError("The sizes of the crystal for the absorption correction should be positive.")
    if flag_warn:
        if any(numpy.any(dict_diffrn.get(f"flags_{name:}", False)) for name in l_unused):
            warnings.warn(
                f"Absorption correction '{correction_type:}': the refined sizes or angles of "
                "exptl_crystal that it does not use (size_rad for 'gaussian', size_max/mid/min "
                "and the Euler angles for 'sphere') do not change the model.",
                UserWarning, stacklevel=2)
        if mu * semi_axes.max() > MU_R_ACCURATE:
            warnings.warn(
                f"Absorption: mu R = {mu * semi_axes.max():.3g} is large; the numerical "
                "transmission is accurate to about 1 % up to mu R = 4 only.",
                UserWarning, stacklevel=2)
        if flag_half_wavelength and mu_composition is None:
            warnings.warn(
                "Absorption with lambda/2 contamination: no crystal structure, mu at lambda/2 is "
                "taken equal to mu.", UserWarning, stacklevel=2)
    from cryspy.A_functions_base.unit_cell import (
        calc_eq_ccs_by_unit_cell_parameters, calc_sthovl_by_unit_cell_parameters)
    sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    q_ccs = 2. * sthovl[na, :] * eq_ccs
    if "angle_omega_chi_phi" in dict_diffrn:
        omega, chi, phi = dict_diffrn["angle_omega_chi_phi"]
        e_in, e_out, misfit = calc_beam_directions_by_angles(
            q_ccs, dict_diffrn["matrix_u"], wavelength, omega, chi, phi)
        if flag_warn and numpy.any(~(misfit < MISFIT_TOLERANCE)):
            warnings.warn(
                f"Absorption: {int((~(misfit < MISFIT_TOLERANCE)).sum()):} reflections are not in the "
                "diffraction condition at their angles omega, chi, phi (check U and the "
                "convention of the angles: omega is the angle of the sample, not omega - theta).",
                UserWarning, stacklevel=2)
    else:
        e_in, e_out = calc_beam_directions_normal_beam(
            q_ccs, dict_diffrn["matrix_u"][6:9], wavelength,
            str(dict_diffrn.get("absorption_diffracted_beam", "left")))[:2]
    rule = UNIT_BALL_RULE if mu * semi_axes.max() <= 2. else UNIT_BALL_RULE_FINE
    transmission = calc_transmission_ellipsoid(mu, semi_axes, matrix_axes, e_in, e_out, rule)
    transmission_half = None
    if flag_half_wavelength:
        # mu_abs is proportional to the wavelength, the scattering part is not
        mu_half = mu * ratio_half
        transmission_half = calc_transmission_ellipsoid(
            mu_half, semi_axes, matrix_axes, e_in, e_out, rule)
    return transmission, transmission_half


L_ABSORPTION_PARAMETERS = ("absorption_mu", "absorption_radius", "absorption_size",
                           "absorption_euler")


def calc_transmission_derivatives_by_dictionary(
        dict_diffrn: dict, index_hkl, unit_cell_parameters, wavelength,
        flag_half_wavelength: bool = False, step_relative: float = 1e-4,
        dict_crystal: dict = None) -> dict:
    """Numerical derivatives of the transmission over the refined
    parameters of the absorption correction: {name: (der, der_half)},
    arrays (n, n_component), for the names of L_ABSORPTION_PARAMETERS with
    flags; der_half (of the lambda/2 transmission) is None without
    flag_half_wavelength. Central differences, the step is relative to the
    value (0.01 deg for the Euler angles)."""
    d_der = {}
    for name in L_ABSORPTION_PARAMETERS:
        flags = dict_diffrn.get(f"flags_{name:}")
        if flags is None or not numpy.any(flags):
            continue
        value = numpy.array(numpy.atleast_1d(dict_diffrn[name]), dtype=float)
        der = numpy.zeros((index_hkl.shape[1], value.size), dtype=float)
        der_half = numpy.zeros_like(der) if flag_half_wavelength else None
        for i_comp in range(value.size):
            if not numpy.atleast_1d(flags)[i_comp]:
                continue
            step = 1e-2 if name == "absorption_euler" else max(abs(value[i_comp]) * step_relative, 1e-8)
            l_res = []
            for sign in (1., -1.):
                value_shift = value.copy()
                value_shift[i_comp] += sign * step
                l_res.append(calc_transmission_by_dictionary(
                    dict_diffrn, index_hkl, unit_cell_parameters, wavelength,
                    flag_half_wavelength=flag_half_wavelength,
                    d_parameter={name: value_shift}, flag_warn=False,
                    dict_crystal=dict_crystal))
            der[:, i_comp] = (l_res[0][0] - l_res[1][0]) / (2. * step)
            if flag_half_wavelength:
                der_half[:, i_comp] = (l_res[0][1] - l_res[1][1]) / (2. * step)
        d_der[name] = (der, der_half)
    return d_der
