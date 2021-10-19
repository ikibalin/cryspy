import numpy


def calc_extinction_sphere_primary(
        f_sq, radius, volume_unit_cell, cos_2theta, wavelength,
        flag_f_sq: bool = False, flag_radius: bool = False,
        flag_volume_unit_cell: bool = False, flag_cos_2theta: bool = False,
        flag_wavelength: bool = False):
    """Primary extinction for the sphere (Acta.Cryst.(1974), A30, 129)

    """
    x = 1.5 * f_sq * numpy.square(wavelength*radius/volume_unit_cell)
    A = 0.20 + 0.45 * cos_2theta
    B = 0.22 - 0.12 * (0.5-cos_2theta)**2 
    x_sq = numpy.square(x)
    y_p = 1./numpy.sqrt(1 + 2 * x + A * x_sq / (1.+ B*x))
    dder = {}
    der_yp_x = -0.5*numpy.power(y_p, 3)*(
        2. + 2.*A*x / (1.+B*x) - A*B*x_sq / numpy.square(1.+ B*x))
    if flag_f_sq:
        ones_f_sq = numpy.ones_like(f_sq)
        dder["f_sq"] = der_yp_x * 1.5 * numpy.square(
            wavelength*radius/volume_unit_cell)*ones_f_sq
    if flag_radius:
        ones_radius = numpy.ones_like(radius)
        dder["radius"] = der_yp_x * 3.0 * f_sq * radius*numpy.square(
            wavelength/volume_unit_cell) 

    if flag_volume_unit_cell:
        ones_volume_unit_cell = numpy.ones_like(volume_unit_cell)
        dder["volume_unit_cell"] = -3.0 * der_yp_x *  f_sq * numpy.square(
            wavelength*radius/volume_unit_cell) * \
            ones_volume_unit_cell / volume_unit_cell
        pass
    if flag_cos_2theta:
        ones_cos_2theta = numpy.ones_like(cos_2theta)
        der_A_cos_2theta = 0.45 
        der_B_cos_2theta = 0.24 * (0.5-cos_2theta) 
        der_yp_A = -0.5*numpy.power(y_p, 3)*x_sq / (1.+ B*x)
        der_yp_B = 0.5*numpy.power(y_p, 3)*A * x_sq*x / numpy.square(1.+ B*x)
        dder["cos_2theta"] = (der_yp_A*der_A_cos_2theta +
                              der_yp_B*der_B_cos_2theta) * ones_cos_2theta
    if flag_wavelength:
        ones_wavelength = numpy.ones_like(wavelength)
        dder["wavelength"] = der_yp_x * 3.0 * f_sq * wavelength*numpy.square(
            radius/volume_unit_cell) * ones_wavelength
    return y_p, dder


def calc_extinction_sphere_secondary_gauss(
        f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
        flag_f_sq: bool = False, flag_radius: bool = False,
        flag_mosaicity: bool = False,
        flag_volume_unit_cell: bool = False, flag_cos_2theta: bool = False,
        flag_wavelength: bool = False):
    """Secondary Gauss extinction for the sphere (Acta.Cryst.(1974), A30, 129)
    mosaicity are given in minutes
    """
    mosaicity_rad = mosaicity*numpy.pi/(180*60) #transfer minutes to radians
    if not(radius*mosaicity_rad>0):
        term_1 = numpy.zeros_like(cos_2theta)
        x = numpy.zeros_like(term_1)
    else:
        term_1 = mosaicity_rad / numpy.sqrt(numpy.square(mosaicity_rad * wavelength) +
                                            9./8. * numpy.square(radius) *
                                            (1.- numpy.square(cos_2theta)))
        x = 1.5 * f_sq * wavelength * \
            numpy.square(wavelength*radius/volume_unit_cell) * term_1
    A = 0.58 + 0.48 * cos_2theta + 0.24 * numpy.square(cos_2theta)
    B = 0.02 - 0.025 * cos_2theta
    x_sq = numpy.square(x)
    y_s = 1./numpy.sqrt(1 + 2.12 * x + A * x_sq / (1.+ B*x))
    y_s[numpy.isnan(y_s)]=1. 
    dder = {}
    der_ys_x = -0.5*numpy.power(y_s, 3)*(
        2.12 + 2.*A*x / (1.+B*x) - A*B*x_sq / numpy.square(1.+ B*x))
    if flag_f_sq:
        ones_f_sq = numpy.ones_like(f_sq)
        dder["f_sq"] = der_ys_x * 1.5 * wavelength*numpy.square(
            wavelength*radius/volume_unit_cell)*term_1*ones_f_sq
    if flag_radius:
        ones_radius = numpy.ones_like(radius)
        dder["radius"] = der_ys_x*(2*radius*1.5 * f_sq * wavelength * \
            numpy.square(wavelength/volume_unit_cell)*term_1+
            1.5 * f_sq * wavelength * \
            numpy.square(wavelength*radius/volume_unit_cell) * mosaicity_rad / \
            numpy.power(numpy.sqrt(numpy.square(mosaicity_rad * wavelength) + \
            9./8. * numpy.square(radius) * (1.- numpy.square(cos_2theta))),3) *
            9./8. * 2*radius * (1.- numpy.square(cos_2theta))
            )
    if flag_mosaicity:
        dder["mosaicity"] = numpy.pi/(180*60)*der_ys_x*1.5 * f_sq * wavelength * \
            numpy.square(wavelength*radius/volume_unit_cell)*(
                 1./ numpy.sqrt(numpy.square(mosaicity_rad * wavelength) +
                    9./8. * numpy.square(radius) *(1.- numpy.square(cos_2theta))) + 
                                            mosaicity_rad / \
            numpy.power(numpy.sqrt(numpy.square(mosaicity_rad * wavelength) + \
            9./8. * numpy.square(radius) * (1.- numpy.square(cos_2theta))),3)*2.*numpy.square(wavelength)*mosaicity_rad)
    if flag_volume_unit_cell:
        ones_volume_unit_cell = numpy.ones_like(volume_unit_cell)
        dder["volume_unit_cell"] = 0*ones_volume_unit_cell
    if flag_cos_2theta:
        ones_cos_2theta = numpy.ones_like(cos_2theta)
        dder["cos_2theta"] = 0*ones_cos_2theta
    if flag_wavelength:
        ones_wavelength = numpy.ones_like(wavelength)
        dder["wavelength"] = 0*ones_wavelength
    return y_s, dder


def calc_extinction_sphere_secondary_lorentz(
        f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
        flag_f_sq: bool = False, flag_radius: bool = False,
        flag_mosaicity: bool = False,
        flag_volume_unit_cell: bool = False, flag_cos_2theta: bool = False,
        flag_wavelength: bool = False):
    """Secondary Gauss extinction for the sphere (Acta.Cryst.(1974), A30, 129)
    mosaicity are given in minutes
    """
    mosaicity_rad = mosaicity*numpy.pi/(180*60) #transfer minutes to radians
    if not(radius*mosaicity_rad>0):
        term_1 = numpy.zeros_like(cos_2theta)
        x = numpy.zeros_like(term_1)
    else:
        term_1 = mosaicity_rad / (mosaicity_rad * wavelength + radius * numpy.sqrt(
            1.- numpy.square(cos_2theta)))
        x = 1.5 * f_sq * wavelength * \
            numpy.square(wavelength*radius/volume_unit_cell) * term_1
    A = 0.025 + 0.285 * cos_2theta
    B_1 = -0.45 * cos_2theta # cos_2theta <= 0.
    B_2 = 0.15 - 0.2 * numpy.square(0.75 - cos_2theta) # cos_2theta > 0.
    B = numpy.where(cos_2theta> 0., B_2, B_1)
    x_sq = numpy.square(x)
    y_s = 1./numpy.sqrt(1 + 2 * x + A * x_sq / (1.+ B*x))
    dder = {}
    der_ys_x = -0.5*numpy.power(y_s, 3)*(
        2. + 2.*A*x / (1.+B*x) - A*B*x_sq / numpy.square(1.+ B*x))
    if flag_f_sq:
        ones_f_sq = numpy.ones_like(f_sq)
        dder["f_sq"] = der_ys_x * 1.5 * wavelength*numpy.square(
            wavelength*radius/volume_unit_cell)*term_1*ones_f_sq
    if flag_radius:
        ones_radius = numpy.ones_like(radius)
        dder["radius"] = 0*ones_radius
    if flag_mosaicity:
        ones_mosaicity = numpy.ones_like(mosaicity)
        dder["mosaicity"] = 0*ones_mosaicity
    if flag_volume_unit_cell:
        ones_volume_unit_cell = numpy.ones_like(volume_unit_cell)
        dder["volume_unit_cell"] = 0*ones_volume_unit_cell
    if flag_cos_2theta:
        ones_cos_2theta = numpy.ones_like(cos_2theta)
        dder["cos_2theta"] = 0*ones_cos_2theta
    if flag_wavelength:
        ones_wavelength = numpy.ones_like(wavelength)
        dder["wavelength"] = 0*ones_wavelength
    return y_s, dder


def calc_extinction_sphere(
        f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
        model: str, flag_f_sq: bool = False, flag_radius: bool = False,
        flag_mosaicity: bool = False,
        flag_volume_unit_cell: bool = False, flag_cos_2theta: bool = False,
        flag_wavelength: bool = False):

    """Extinction correction for the sphere (Acta.Cryst.(1974), A30, 129)
    
    model: "gauss" or "lorentz"
    mosaicity are given in minutes
    """
    y_p, dder_y_p = calc_extinction_sphere_primary(
        f_sq, radius, volume_unit_cell, cos_2theta, wavelength,
        flag_f_sq=flag_f_sq, flag_radius=flag_radius,
        flag_volume_unit_cell=flag_volume_unit_cell,
        flag_cos_2theta=flag_cos_2theta, flag_wavelength=flag_wavelength)
    if model.lower().startswith("gauss"):
        y_s, dder_y_s = calc_extinction_sphere_secondary_gauss(
            f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
            flag_f_sq=flag_f_sq, flag_radius=flag_radius,
            flag_mosaicity=flag_mosaicity,
            flag_volume_unit_cell=flag_volume_unit_cell,
            flag_cos_2theta=flag_cos_2theta, flag_wavelength=flag_wavelength)
    elif model.lower().startswith("lorentz"):
        y_s, dder_y_s = calc_extinction_sphere_secondary_lorentz(
            f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
            flag_f_sq=flag_f_sq, flag_radius=flag_radius,
            flag_mosaicity=flag_mosaicity,
            flag_volume_unit_cell=flag_volume_unit_cell,
            flag_cos_2theta=flag_cos_2theta, flag_wavelength=flag_wavelength)
    else:
        y_s = numpy.ones(y_p.shape, dtype=float)
        dder_y_s = {
            "f_sq": numpy.zeros_like(f_sq),
            "radius": numpy.zeros_like(radius),
            "mosaicity": numpy.zeros_like(mosaicity),
            "volume_unit_cell": numpy.zeros_like(volume_unit_cell),
            "cos_2theta": numpy.zeros_like(cos_2theta),
            "wavelength": numpy.zeros_like(wavelength)}
    y = y_p * y_s
    dder = {}
    if flag_f_sq:
        dder["f_sq"] = y_p * dder_y_s["f_sq"] + dder_y_p["f_sq"] * y_s
    if flag_radius:
        dder["radius"] = y_p * dder_y_s["radius"] + dder_y_p["radius"] * y_s
    if flag_mosaicity:
        dder["mosaicity"] = y_p * dder_y_s["mosaicity"] 
    if flag_volume_unit_cell:
        dder["volume_unit_cell"] = y_p * dder_y_s["volume_unit_cell"] + \
            dder_y_p["volume_unit_cell"] * y_s
    if flag_cos_2theta:
        dder["cos_2theta"] = y_p * dder_y_s["cos_2theta"] + \
            dder_y_p["cos_2theta"] * y_s
    if flag_wavelength:
        dder["wavelength"] = y_p * dder_y_s["wavelength"] + \
            dder_y_p["wavelength"] * y_s
    return y, dder
