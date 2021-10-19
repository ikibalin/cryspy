import os
import numpy

from cryspy.A_functions_base.extinction import \
    calc_extinction_sphere, \
    calc_extinction_sphere_primary, \
    calc_extinction_sphere_secondary_gauss, \
    calc_extinction_sphere_secondary_lorentz



f_sq = numpy.array([4, 5], dtype=float)
radius = numpy.array([100], dtype=float)
mosaicity = numpy.array([200], dtype=float)*(180*60)/numpy.pi
volume_unit_cell = numpy.array([1000], dtype=float)
cos_2theta = numpy.array([0.4, 0.3], dtype=float)
wavelength = numpy.array([0.8], dtype=float)
model_extinction_g = "gauss"
model_extinction_l = "lorentz"

flag_f_sq = True
flag_radius = True
flag_mosaicity = True
flag_volume_unit_cell = True
flag_cos_2theta = True
flag_wavelength = True

y_g_s   = numpy.array([0.93118955, 0.91576117], dtype=float)
y_l_s   = numpy.array([0.94069951, 0.92734411], dtype=float)
y_s_p   = numpy.array([0.96343063, 0.95486817 ], dtype=float)
y_s_s_g = numpy.array([0.96653513, 0.9590446], dtype=float)
y_s_s_l = numpy.array([0.97640606, 0.97117502], dtype=float)


def test_calc_extinction_sphere():
    res, dder = calc_extinction_sphere(f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
        model_extinction_g, flag_f_sq=flag_f_sq, flag_radius=flag_radius,
        flag_mosaicity=flag_mosaicity,
        flag_volume_unit_cell=flag_volume_unit_cell, flag_cos_2theta=flag_cos_2theta,
        flag_wavelength=flag_wavelength)
    print("res: ", res)
    print(y_g_s)

    assert numpy.all(numpy.isclose(res, y_g_s))

    res, dder = calc_extinction_sphere(f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
        model_extinction_l, flag_f_sq=flag_f_sq, flag_radius=flag_radius,
        flag_mosaicity=flag_mosaicity,
        flag_volume_unit_cell=flag_volume_unit_cell, flag_cos_2theta=flag_cos_2theta,
        flag_wavelength=flag_wavelength)
    print("res: ", res)
    print(y_l_s)

    assert numpy.all(numpy.isclose(res, y_l_s))


def test_calc_extinction_sphere_primary():
    res, dder = calc_extinction_sphere_primary(f_sq, radius, volume_unit_cell, cos_2theta, wavelength,
    flag_f_sq=flag_f_sq, flag_radius=flag_radius, flag_volume_unit_cell=flag_volume_unit_cell,
    flag_cos_2theta=flag_cos_2theta, flag_wavelength=flag_wavelength)
    print("res: ", res)
    print(y_s_p)

    assert numpy.all(numpy.isclose(res, y_s_p))


def test_calc_extinction_sphere_secondary_gauss():
    res, dder = calc_extinction_sphere_secondary_gauss(f_sq, radius, mosaicity, volume_unit_cell, 
        cos_2theta, wavelength, flag_f_sq=flag_f_sq, flag_radius=flag_radius, flag_mosaicity=flag_mosaicity,
        flag_volume_unit_cell=flag_volume_unit_cell, flag_cos_2theta=flag_cos_2theta, flag_wavelength=flag_wavelength)
    print("res: ", res)
    print(y_s_s_g)

    assert numpy.all(numpy.isclose(res, y_s_s_g))


def test_calc_extinction_sphere_secondary_lorentz():
    res, dder = calc_extinction_sphere_secondary_lorentz(f_sq, radius, mosaicity, volume_unit_cell,
        cos_2theta, wavelength, flag_f_sq=flag_f_sq, flag_radius=flag_radius,
        flag_mosaicity=flag_mosaicity, flag_volume_unit_cell=flag_volume_unit_cell, 
        flag_cos_2theta=flag_cos_2theta, flag_wavelength=flag_wavelength)
    print("res: ", res)
    print(y_s_s_l)

    assert numpy.all(numpy.isclose(res, y_s_s_l))


