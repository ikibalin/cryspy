"""Small parts of the FullProf conversion: types of atoms, labels, lines of
the .pcr file, background polynomial and data files."""
import numpy
import pytest

from cryspy.H_functions_global.fullprof import (
    FullProfFormatError, PcrLines, calc_fullprof_polynomial_background, form_fullprof_type,
    form_value, read_fullprof_data, take_label, take_scale_factor, take_type_symbol)


@pytest.mark.parametrize("fullprof_type, type_symbol", [
    ("PB", "Pb"), ("PB+2", "Pb2+"), ("O-2", "O2-"), ("s", "S"), ("NA+", "Na1+")])
def test_take_type_symbol(fullprof_type, type_symbol):
    assert take_type_symbol(fullprof_type) == type_symbol


def test_take_type_symbol_magnetic_form_factor():
    with pytest.raises(FullProfFormatError):
        take_type_symbol("MMN3")


def test_form_fullprof_type():
    assert form_fullprof_type("Pb2+", "X-rays") == "PB+2"
    assert form_fullprof_type("Pb2+", "neutrons") == "PB"
    assert form_fullprof_type("O", "X-rays") == "O"


def test_take_label():
    l_used = ["pbso4"]
    assert take_label("PbSO4 neutrons", l_used) == "pbso4_2"
    # 'data_' would start a new data block of CIF
    assert take_label("data_pd", []) == "pd"
    assert take_label("2phase", []) == "phase_2phase"


def test_form_value():
    assert form_value(1.25, 0.) == "1.25"
    assert form_value(1.25, 21.) == "1.25()"
    assert form_value(3., 0.) == "3.0"


def test_pcr_lines():
    lines = PcrLines("! comment\n\n  1.5  2  0.3  Number of refined\n  4 5 # text\n")
    assert lines.next_numbers(3) == [1.5, 2., 0.3]
    assert lines.next_integers() == [4, 5]
    with pytest.raises(FullProfFormatError):
        lines.next_text()


def test_polynomial_background():
    """B_0 + B_1 x + B_2 x^2 with x = 2theta/Bkpos - 1."""
    ttheta = numpy.array([40., 80., 120.])
    background = calc_fullprof_polynomial_background(ttheta, [10., 2., 1., 0., 0., 0.], 80.)
    x = ttheta / 80. - 1.
    assert numpy.allclose(background, 10. + 2. * x + x * x)


def test_scale_factor():
    assert take_scale_factor("X-rays") == 2.
    assert take_scale_factor("neutrons") == 1.


def test_read_xydata(tmp_path):
    file_data = tmp_path / "x.dat"
    file_data.write_text("XYDATA\nINTER 1 1 0\n10.0 100.0 10.0\n10.1 121.0 11.0\n")
    ttheta, intensity, sigma = read_fullprof_data(str(file_data), 10)
    assert numpy.allclose(ttheta, [10., 10.1])
    assert numpy.allclose(intensity, [100., 121.])
    assert numpy.allclose(sigma, [10., 11.])


def test_read_ins0(tmp_path):
    file_data = tmp_path / "x.dat"
    file_data.write_text("10.0 0.5 11.0\n 100 400 900\n")
    ttheta, intensity, sigma = read_fullprof_data(str(file_data), 0)
    assert numpy.allclose(ttheta, [10., 10.5, 11.])
    assert numpy.allclose(sigma, [10., 20., 30.])


def test_read_ins6_d1a(tmp_path):
    """D1A/D2B: pairs (I2, I6) of the number of averaged cells and the counts,
    the counts are averages over the cells, sigma = sqrt(counts / cells)."""
    file_data = tmp_path / "x.dat"
    file_data.write_text(
        "D1A test\n  1  1  0.05\n  10.00  0  0\n header\n"
        " 4   400 1   100\n"
        "-1000\n")
    ttheta, intensity, sigma = read_fullprof_data(str(file_data), 6)
    assert numpy.allclose(ttheta, [10., 10.05])
    assert numpy.allclose(intensity, [400., 100.])
    assert numpy.allclose(sigma, [10., 10.])


def test_unknown_data_format(tmp_path):
    file_data = tmp_path / "x.dat"
    file_data.write_text("1 2 3\n")
    with pytest.raises(FullProfFormatError):
        read_fullprof_data(str(file_data), 3)


def test_old_cubic_symbols():
    """FullProf accepts 'F m 3 m', CrysPy needs 'F m -3 m'; 'P 2 3' has no bar."""
    from cryspy.H_functions_global.fullprof import take_space_group_symbol
    assert take_space_group_symbol("F m 3 m") == "F m -3 m"
    assert take_space_group_symbol("I a 3") == "I a -3"
    assert take_space_group_symbol("P 2 3") == "P 2 3"
    assert take_space_group_symbol("P n m a") == "P n m a"


def test_chebyshev_background():
    """Nba = -5: sum B_m T_m(x), x from -1 at Thmin to 1 at Thmax."""
    from cryspy.H_functions_global.fullprof import calc_fullprof_chebyshev_background
    ttheta = numpy.array([10., 60., 110.])
    background = calc_fullprof_chebyshev_background(ttheta, [5., 2., 1.], 10., 110.)
    x = numpy.array([-1., 0., 1.])
    assert numpy.allclose(background, 5. + 2. * x + (2. * x * x - 1.))


def test_zero_counts_have_unit_sigma(tmp_path):
    """Data without sigma: sqrt(I), 1 for zero counts (no infinite weight)."""
    file_data = tmp_path / "x.dat"
    file_data.write_text("XYDATA\n10.0 0.0\n10.1 4.0\n")
    _, _, sigma = read_fullprof_data(str(file_data), 10)
    assert numpy.allclose(sigma, [1., 2.])


def test_read_ins0_with_title(tmp_path):
    """Ins = 0 files may start with a title and comment lines."""
    file_data = tmp_path / "x.dat"
    file_data.write_text("\n! Diffraction pattern\n#\n 10.0 0.5 11.0\n 100 400 900\n")
    ttheta, intensity, _ = read_fullprof_data(str(file_data), 0)
    assert numpy.allclose(ttheta, [10., 10.5, 11.])
    assert numpy.allclose(intensity, [100., 400., 900.])


def test_pcr_lines_skip_vary():
    """VARY and FIX commands of FullProf are not lines of parameters."""
    lines = PcrLines("VARY ratio\nFIX all\n 1.0 2.0\n")
    assert lines.next_numbers(2) == [1., 2.]


def test_more_lines_without_distances():
    """Jdi = 1 with no line of distances before the space group (dy.pcr of
    the examples of FullProf)."""
    from cryspy.H_functions_global.fullprof import read_more_lines
    lines = PcrLines("   3   1   0   0   0   0  1.0000  0.0000  0.0000  0.0000    1\nP 6/m m m\n")
    read_more_lines(lines, "phase")
    assert lines.next_text() == "P 6/m m m"
