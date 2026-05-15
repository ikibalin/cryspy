import os
import numpy
import cryspy

from cryspy.A_functions_base.structure_factor import (
    calc_f_m_perp_ordered_by_dictionary,
)

from cryspy.procedure_simulation.simulation import simulation_single_crystal


def test_tbte_mcif():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "TbTe_63_467m.mcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_rhochi = rhochi.get_dictionary()
    d_crystal = d_rhochi["crystal_tbte"]

    dict_in_out = {}
    dict_in_out["index_hkl"] = numpy.array(
        [
            [
                1,
            ],
            [
                1,
            ],
            [
                0,
            ],
        ],
        dtype=float,
    )
    f_m_perp_o, dder = calc_f_m_perp_ordered_by_dictionary(
        d_crystal, dict_in_out, flag_use_precalculated_data=False
    )
    flag = numpy.any(f_m_perp_o != 0.0)
    assert flag


def test_refine_mcif_2k():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "2k.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    cell = rhochi.crystal_s2k.cell
    assert cell.type_cell.startswith("o")


def test_refine_mcif_1k():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "1k.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    crystal = rhochi.crystal_s1k
    assert numpy.isclose(crystal.atom_site_moment["Gd_2"].crystalaxis_y, 0.0)
    d_out = cryspy.rhochi_rietveld_refinement(rhochi)
    assert numpy.isclose(
        crystal.atom_site_moment["Gd_2"].crystalaxis_x, -5.867102366251087
    )
    assert d_out["chi_sq"] <= 57741.0
