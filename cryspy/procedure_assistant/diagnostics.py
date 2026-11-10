"""Diagnosis of the misfit of the experiments of an object (tool diagnose_fit).

The numbers come from ``A_functions_base.fit_diagnostics``; here the
measured and calculated data are taken from the data blocks after a
calculation: 1D powder (pd_proc), time of flight (tof_proc), single-crystal
intensities and flip ratios (diffrn_refln).
"""
from typing import List

import numpy

from cryspy.A_functions_base.fit_diagnostics import (
    diagnose_flip_ratios, diagnose_powder_profile, diagnose_single_crystal)
from cryspy.procedure_assistant import experiments


def take_columns(loop, names) -> dict:
    """Columns of a loop as float arrays (missing ones are left out)."""
    out = {}
    for name in names:
        if loop.is_attribute(name):
            out[name] = numpy.array(
                [numpy.nan if value is None else value for value in getattr(loop, name)], dtype=float)
    return out


def take_included(loop) -> numpy.ndarray:
    if loop.is_attribute("excluded"):
        return ~numpy.array([bool(value) for value in loop.excluded], dtype=bool)
    return numpy.ones(len(loop.items), dtype=bool)


L_PROFILE_COLUMNS = ("ttheta", "ttheta_corrected", "time", "intensity", "intensity_sigma",
                     "intensity_plus", "intensity_plus_sigma", "intensity_minus",
                     "intensity_minus_sigma", "intensity_plus_net", "intensity_minus_net",
                     "intensity_bkg_calc")


def diagnose_profile(block) -> List[tuple]:
    """Signatures of a 1D powder or time-of-flight experiment (the sum of
    up and down intensities for polarized neutrons)."""
    flag_time = block.is_attribute("tof_proc")
    loop = block.tof_proc if flag_time else block.pd_proc
    columns = take_columns(loop, L_PROFILE_COLUMNS)
    if "intensity" in columns:
        y_obs, sigma = columns["intensity"], columns["intensity_sigma"]
    else:
        y_obs = columns["intensity_plus"] + columns["intensity_minus"]
        sigma = numpy.sqrt(numpy.square(columns["intensity_plus_sigma"])
                           + numpy.square(columns["intensity_minus_sigma"]))
    net = columns["intensity_plus_net"] + columns["intensity_minus_net"]
    y_calc = net + columns["intensity_bkg_calc"]
    flag = take_included(loop) & numpy.isfinite(y_obs) & numpy.isfinite(y_calc)
    if flag_time:
        x, sthovl = columns["time"], None
    else:
        x = columns["ttheta"]
        ttheta = columns.get("ttheta_corrected", x)
        sthovl = numpy.sin(numpy.radians(ttheta) / 2.) / float(block.setup.wavelength)
        sthovl = sthovl[flag]
    return diagnose_powder_profile(x[flag], y_obs[flag], sigma[flag], y_calc[flag], net[flag],
                                   sthovl=sthovl, flag_time=flag_time)


def take_sthovl(cryspy_object, block, index_hkl) -> numpy.ndarray:
    l_crystal = experiments.take_experiment_crystals(cryspy_object, block)
    if len(l_crystal) == 0:
        return numpy.zeros(index_hkl.shape[1])
    return numpy.asarray(l_crystal[0].cell.calc_sthovl(*index_hkl), dtype=float)


def diagnose_diffrn(cryspy_object, block, n_parameter: int) -> List[tuple]:
    """Signatures of single-crystal intensities or flip ratios."""
    loop = block.diffrn_refln
    flag = take_included(loop)
    index_hkl = numpy.array([loop.index_h, loop.index_k, loop.index_l], dtype=float)[:, flag]
    sthovl = take_sthovl(cryspy_object, block, index_hkl)
    if loop.is_attribute("fr"):
        columns = take_columns(loop, ("fr", "fr_sigma", "fr_calc"))
        return diagnose_flip_ratios(index_hkl, columns["fr"][flag], columns["fr_sigma"][flag],
                                    columns["fr_calc"][flag], sthovl, n_parameter)
    columns = take_columns(loop, ("intensity", "intensity_sigma", "intensity_calc"))
    return diagnose_single_crystal(index_hkl, columns["intensity"][flag],
                                   columns["intensity_sigma"][flag],
                                   columns["intensity_calc"][flag], sthovl, n_parameter)


def diagnose_object(cryspy_object, experiment: str = "") -> str:
    """Text of the signatures of every experiment (or of the one named)."""
    n_parameter = len(cryspy_object.get_variable_names())
    ls_out = []
    for block in experiments.take_experiments(cryspy_object):
        name = block.get_name()
        if experiment and experiment not in (name, block.data_name):
            continue
        try:
            if block.is_attribute("pd_proc") or block.is_attribute("tof_proc"):
                l_sign = diagnose_profile(block)
            elif block.is_attribute("diffrn_refln"):
                l_sign = diagnose_diffrn(cryspy_object, block, n_parameter)
            else:
                ls_out.append(f"{name:}: no diagnosis for this kind of experiment.")
                continue
        except (AttributeError, KeyError, TypeError, ValueError) as error:
            ls_out.append(f"{name:}: no diagnosis ({error:}).")
            continue
        ls_out.append(f"{name:}:")
        ls_out.extend(f"- {key:}: {text:}" for key, text in l_sign)
    if not ls_out:
        return f"Error: no experiment '{experiment:}' in the object."
    return "\n".join(ls_out)
