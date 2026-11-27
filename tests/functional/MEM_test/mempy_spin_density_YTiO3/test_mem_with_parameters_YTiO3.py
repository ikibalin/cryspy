import os

import numpy
import pytest

import cryspy

F_NAME = os.path.join(os.path.dirname(__file__), "mempy_spin_density_YTiO3.rcif")


def take_channel_col(rcif_obj):
    return [item for item in rcif_obj.items if isinstance(item, cryspy.ChannelColL)]


def test_mem_with_parameters_adds_channel_col(tmp_path, monkeypatch):
    # the .den files are written to the working directory
    monkeypatch.chdir(tmp_path)
    rcif_obj = cryspy.load_file(F_NAME)
    assert take_channel_col(rcif_obj) == []

    res_short = cryspy.mempy_reconstruction_with_parameters(rcif_obj, iteration_max=2)
    res_long = cryspy.mempy_reconstruction_with_parameters(rcif_obj, iteration_max=40)
    assert res_long["chi_sq"] < res_short["chi_sq"]

    # a second run replaces the loop instead of adding one more
    l_channel_col = take_channel_col(rcif_obj)
    assert len(l_channel_col) == 1
    channel_col = l_channel_col[0]
    density_plus = numpy.array(channel_col.density_plus)
    density_minus = numpy.array(channel_col.density_minus)
    point_multiplicity = numpy.array(channel_col.point_multiplicity)
    assert len(channel_col.items) == point_multiplicity.size
    assert numpy.all(density_plus >= 0.) and numpy.all(density_minus >= 0.)
    # both channels keep the same normalization over the unit cell
    # (relative tolerance for the rounding of the written densities)
    assert (density_plus * point_multiplicity).sum() == pytest.approx(
        (density_minus * point_multiplicity).sum(), rel=1e-4)
    assert not numpy.allclose(density_plus, density_minus)

    for diffrn in [item for item in rcif_obj.items if isinstance(item, cryspy.Diffrn)]:
        assert numpy.all(numpy.isfinite(diffrn.diffrn_refln.numpy_fr_calc))
        assert diffrn.refine_ls.number_reflns > 0
