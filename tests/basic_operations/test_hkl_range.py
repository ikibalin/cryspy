import os
import math
import numpy

import cryspy
from cryspy.A_functions_base.structure_factor import (
    calc_index_hkl_multiplicity_in_range,
)


def test_calc_index_hkl_multiplicity_in_range():
    sthovl_min = 0.001
    sthovl_max = 0.5
    unit_cell_parameters = numpy.array(
        [5.0, 5.0, 5.0, numpy.pi / 2, numpy.pi / 2, numpy.pi / 2], dtype=float
    )
    flag_only_nuclear = False
    spgr = cryspy.SpaceGroup(name_hm_alt="F d -3 m")
    rsgs = spgr.reduced_space_group_symop
    np_rsgs = rsgs.get_symm_elems()
    shift = spgr.shift

    l_vals = []
    for ind in range(len(shift)):
        lcm = numpy.lcm.reduce([fr.denominator for fr in shift[ind]])
        vals = [int(fr.numerator * lcm / fr.denominator) for fr in shift[ind]]
        vals.append(lcm)
        l_vals.append(vals)
    translation_elems = numpy.array(l_vals, dtype=int)
    centrosymmetry = spgr.centrosymmetry
    hkl_1, _ = calc_index_hkl_multiplicity_in_range(
        sthovl_min,
        sthovl_max,
        unit_cell_parameters,
        np_rsgs,
        translation_elems,
        centrosymmetry,
        flag_only_nuclear=flag_only_nuclear,
    )
    flag_only_nuclear = True
    hkl_2, _ = calc_index_hkl_multiplicity_in_range(
        sthovl_min,
        sthovl_max,
        unit_cell_parameters,
        np_rsgs,
        translation_elems,
        centrosymmetry,
        flag_only_nuclear=flag_only_nuclear,
    )

    assert hkl_1.shape[1] != hkl_2.shape[1]
