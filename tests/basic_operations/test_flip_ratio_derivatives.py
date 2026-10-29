"""Analytical derivatives of function_4_flip_ratio.calc_flip_ratio against
central finite differences."""
import numpy
import pytest

from cryspy.A_functions_base.function_4_flip_ratio import calc_flip_ratio

F_NUCL = numpy.array([1.0+0.5j, 2.0-0.3j, 0.4+0.0j])
F_MAG_PERP = (numpy.array([0.3+0.1j, 0.2j, -0.5+0.1j]),
              numpy.array([0.1+0.0j, 0.4+0.2j, 0.3-0.2j]),
              numpy.array([0.2-0.1j, 0.1+0.0j, 0.6+0.3j]))
E_UP = (numpy.array([0.1, 0.3, -0.2]), numpy.array([0.2, 0.1, 0.4]),
        numpy.array([0.97, 0.95, 0.89]))
STEP = 1e-6


def take_flip_ratio(f_nucl=F_NUCL, f_mag_perp=F_MAG_PERP, e_up=E_UP,
                    beam_polarization=0.9, flipper_efficiency=0.95,
                    model_extinction="gauss", **kwargs):
    return calc_flip_ratio(
        f_nucl, f_mag_perp, e_up, beam_polarization, flipper_efficiency,
        10., 20., model_extinction, 100., numpy.array([0.2, 0.3, 0.4]), 0.8,
        **kwargs)


def shift_component(vector, ind, delta):
    return tuple(v + delta if i == ind else v for i, v in enumerate(vector))


@pytest.mark.parametrize("model_extinction", ["gauss", "lorentz"])
def test_derivatives_flip_ratio(model_extinction):
    flags = dict(flag_f_nucl=True, flag_f_mag_perp=True, flag_e_up=True,
                 flag_beam_polarization=True, flag_flipper_efficiency=True)
    dder = take_flip_ratio(model_extinction=model_extinction, **flags)[1]

    def check(analytical, name, plus, minus):
        value_plus = take_flip_ratio(model_extinction=model_extinction, **{name: plus})[0]
        value_minus = take_flip_ratio(model_extinction=model_extinction, **{name: minus})[0]
        assert numpy.allclose(analytical, (value_plus-value_minus)/(2.*STEP), rtol=1e-5, atol=1e-7)

    check(dder["der_flip_ratio__beam_polarization"], "beam_polarization", 0.9+STEP, 0.9-STEP)
    check(dder["der_flip_ratio__flipper_efficiency"], "flipper_efficiency", 0.95+STEP, 0.95-STEP)
    check(dder["der_flip_ratio__f_nucl_re"], "f_nucl", F_NUCL+STEP, F_NUCL-STEP)
    check(dder["der_flip_ratio__f_nucl_im"], "f_nucl", F_NUCL+1j*STEP, F_NUCL-1j*STEP)
    for ind in range(3):
        check(dder["der_flip_ratio__f_mag_perp_re"][ind], "f_mag_perp",
              shift_component(F_MAG_PERP, ind, STEP), shift_component(F_MAG_PERP, ind, -STEP))
        check(dder["der_flip_ratio__f_mag_perp_im"][ind], "f_mag_perp",
              shift_component(F_MAG_PERP, ind, 1j*STEP), shift_component(F_MAG_PERP, ind, -1j*STEP))
        check(dder["der_flip_ratio__e_up"][ind], "e_up",
              shift_component(E_UP, ind, STEP), shift_component(E_UP, ind, -STEP))


def test_no_derivatives_by_default():
    flip_ratio, dder = take_flip_ratio()
    assert flip_ratio.shape == F_NUCL.shape
    assert dder == {}
