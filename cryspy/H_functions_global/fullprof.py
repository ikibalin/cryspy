"""Conversion between the files of FullProf (.pcr and the data file) and rcif.

The conversion covers what both programs calculate in the same way: a
powder diffraction pattern measured at a constant wavelength (neutrons or
X-rays), one pattern, any number of nuclear phases with the
Thompson-Cox-Hastings pseudo-Voigt profile (FullProf Npr = 7, or Npr = 0
with X = Y = 0), isotropic or anisotropic atomic displacements, the
background by points or by the polynomial of FullProf (Nba = 0), excluded
regions, zero shift and SySin, the empirical asymmetry (Asy1-Asy4) and the
axial divergence (S_L, D_L). The zero shift of FullProf is offset_2theta
of CrysPy with the same sign; for X-rays the scale of CrysPy is twice the
one of FullProf (see take_scale_factor).

    >>> import cryspy
    >>> obj = cryspy.load_fullprof("pbso4.pcr")       # pbso4.dat is read too
    >>> cryspy.save_fullprof_files(obj, "copy/pbso4.pcr")  # and copy/pbso4.dat

What one of the programs cannot calculate (magnetic phases of FullProf,
polarized neutrons of CrysPy, TOF, single crystals, ...) is refused with a
clear message instead of a silent approximation. A value that the
conversion cannot keep exactly (a second wavelength, the limit of the
empirical asymmetry, constraints between parameters with the same code)
gives a warning.

References: J. Rodriguez-Carvajal, FullProf manual (fullprof.pcr
description, the input control file); Physica B 192 (1993) 55.
"""
import logging
import math
import os
import re
import warnings
from typing import List, Optional, Tuple

import numpy

from cryspy.exceptions import CrysPyWarning
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.H_functions_global.function_1_cryspy_objects import str_to_globaln

logger = logging.getLogger(__name__)

N_BACKGROUND_POINTS = 25  # points to describe a polynomial background of FullProf

# FullProf Job: 0, 2 X-rays; 1, 3 neutrons; 2, 3 only calculate the pattern
D_JOB_RADIATION = {0: "X-rays", 1: "neutrons", 2: "X-rays", 3: "neutrons"}
JOBS_SIMULATION = (2, 3)


def take_scale_factor(radiation: str) -> float:
    """Scale of CrysPy over the scale of FullProf.

    For X-rays the Lorentz-polarization factor of CrysPy has
    (1 - K + K cthm cos^2 2theta) with K = 0.5, that of FullProf
    (1 + Cthm cos^2 2theta): the scales differ by 2 (the pattern of
    examples/rhochi_unpolarized_powder_1d_xrays_PbSO4 calculated by
    FullProf is reproduced within 2 % with this factor).
    """
    return 2. if radiation == "X-rays" else 1.


class FullProfFormatError(ValueError):
    """The FullProf file has a form or an option the converter does not read."""


def warn_conversion(message: str):
    warnings.warn(message, CrysPyWarning, stacklevel=3)


class PcrLines:
    """The lines of a .pcr file without the comments (lines starting with '!')."""

    def __init__(self, text: str):
        self.lines = [
            line.rstrip() for line in text.splitlines()
            if line.strip() != "" and not line.lstrip().startswith("!")]
        self.position = 0

    def next_text(self) -> str:
        if self.position >= len(self.lines):
            raise FullProfFormatError("The .pcr file ends too early.")
        line = self.lines[self.position]
        self.position += 1
        return line

    def next_numbers(self, n_min: int = 1) -> List[float]:
        """The numbers of the next line (a text after '#' or after the
        numbers, like 'Number of refined parameters', is dropped)."""
        line = self.next_text().split("#")[0]
        l_value = []
        for word in line.split():
            try:
                l_value.append(float(word))
            except ValueError:
                break
        if len(l_value) < n_min:
            raise FullProfFormatError(
                f"A line with at least {n_min:} numbers is expected in the .pcr "
                f"file, found: '{line.strip():}'.")
        return l_value

    def next_integers(self, n_min: int = 1) -> List[int]:
        return [int(round(value)) for value in self.next_numbers(n_min)]


def is_refined(code: float) -> bool:
    return abs(code) > 1e-6


def form_value(value: float, code: float = 0., digits: int = 6) -> str:
    """A value of rcif, refined ones are marked by '()'."""
    s_value = f"{value:.{digits:}f}".rstrip("0")
    if s_value.endswith("."):
        s_value += "0"
    return s_value + ("()" if is_refined(code) else "")


def take_type_symbol(fullprof_type: str) -> str:
    """'PB' -> 'Pb', 'PB+2' -> 'Pb2+', 'O-2' -> 'O2-' (the type of the
    scattering factor of FullProf to type_symbol of CIF)."""
    match = re.fullmatch(r"([A-Za-z]{1,2})([+-])?(\d)?", fullprof_type.strip())
    if match is None:
        raise FullProfFormatError(
            f"The type of atom '{fullprof_type:}' is not an element or an ion "
            "(magnetic form factors and user scattering factors are not converted).")
    element = match.group(1).capitalize()
    sign, charge = match.group(2), match.group(3)
    if sign is None:
        return element
    return f"{element:}{charge or '1':}{sign:}"


def form_fullprof_type(type_symbol: str, radiation: str) -> str:
    """'Pb2+' -> 'PB+2' for X-rays and 'PB' for neutrons."""
    match = re.fullmatch(r"([A-Za-z]{1,2})(\d)?([+-])?", type_symbol.strip())
    if match is None:
        return type_symbol.upper()
    element = match.group(1).upper()
    if radiation == "neutrons" or match.group(3) is None:
        return element
    return f"{element:}{match.group(3):}{match.group(2) or '1':}"


def calc_fullprof_polynomial_background(
        ttheta: numpy.ndarray, coefficients: List[float], bkpos: float) -> numpy.ndarray:
    """Background of FullProf for Nba = 0, sum_m B_m (2theta/Bkpos - 1)^m,
    m = 0..5 (FullProf manual, background polynomial); 2theta and Bkpos in
    degrees."""
    x = numpy.asarray(ttheta, dtype=float) / bkpos - 1.
    return sum(coefficient * numpy.power(x, power)
               for power, coefficient in enumerate(coefficients))


def read_prf_columns(file_prf: str) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """2theta, Ycal and the background of a .prf file of FullProf written
    with Prf = 3 or -3 (columns 2Theta, Yobs, Ycal, Yobs-Ycal, Backg)."""
    l_row = []
    with open(file_prf, "r") as fid:
        flag_data = False
        for line in fid:
            if line.strip().startswith("2Theta"):
                flag_data = True
                continue
            if not flag_data:
                continue
            l_word = line.split()
            if len(l_word) < 5:
                continue
            try:
                l_row.append([float(word) for word in l_word[:5]])
            except ValueError:
                continue
    if len(l_row) == 0:
        raise FullProfFormatError(
            f"No columns '2Theta Yobs Ycal Yobs-Ycal Backg' in {file_prf:} "
            "(write it with Prf = 3).")
    array = numpy.array(l_row, dtype=float)
    return array[:, 0], array[:, 2], array[:, 4]


def read_xydata(file_data: str, l_line: List[str]) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Ins = 10: X Y [sigma] after the first line (XYDATA)."""
    l_row = []
    for line in l_line[1:]:
        l_word = line.split("!")[0].split()
        try:
            l_value = [float(word) for word in l_word]
        except ValueError:
            continue  # keywords of XYDATA (INTER, TEMP, ...)
        if len(l_value) >= 2:
            l_row.append(l_value[:3])
    if len(l_row) == 0:
        raise FullProfFormatError(f"No data in {file_data:}.")
    ttheta = numpy.array([row[0] for row in l_row], dtype=float)
    intensity = numpy.array([row[1] for row in l_row], dtype=float)
    if all(len(row) >= 3 for row in l_row):
        sigma = numpy.array([row[2] for row in l_row], dtype=float)
    else:
        sigma = numpy.sqrt(numpy.abs(intensity))
    return ttheta, intensity, sigma


def read_start_step_end(l_line: List[str]) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Ins = 0: '2theta_start step 2theta_end' then the intensities."""
    l_number = " ".join(l_line).split()
    start, step, end = [float(word) for word in l_number[:3]]
    intensity = numpy.array([float(word) for word in l_number[3:]], dtype=float)
    ttheta = start + step * numpy.arange(intensity.size)
    flag_range = ttheta <= end + 0.5 * step
    intensity = intensity[flag_range]
    return ttheta[flag_range], intensity, numpy.sqrt(numpy.abs(intensity))


def read_d1a(l_line: List[str]) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Ins = 6: D1A/D2B of ILL, fixed format: pairs (I2, I6) of the number of
    averaged cells and the counts, the data end with a negative number."""
    step = float(l_line[1].split()[2])
    start = float(l_line[2].split()[0])
    l_cell, l_count = [], []
    for line in l_line[4:]:
        if line.strip().startswith("-"):
            break
        for position in range(0, len(line.rstrip()), 8):
            pair = line[position:position + 8]
            if pair.strip() == "":
                continue
            l_cell.append(float(pair[:2] or 0.))
            l_count.append(float(pair[2:8]))
    n_cell = numpy.array(l_cell, dtype=float)
    intensity = numpy.array(l_count, dtype=float)
    ttheta = start + step * numpy.arange(intensity.size)
    # the counts are averaged over n_cell detector cells
    sigma = numpy.sqrt(numpy.abs(intensity) / numpy.where(n_cell > 0, n_cell, 1.))
    return ttheta, intensity, sigma


def read_fullprof_data(file_data: str, ins: int) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """2theta (degrees), intensity and its sigma of a data file of FullProf.

    Ins = 10: X Y [sigma] (XYDATA), Ins = 0: '2theta_start step
    2theta_end' then the intensities, Ins = 6: D1A/D2B of ILL (pairs of the
    number of averaged cells and the counts). Without a sigma it is sqrt(Y).
    """
    with open(file_data, "r") as fid:
        l_line = fid.read().splitlines()
    if ins == 10:
        return read_xydata(file_data, l_line)
    if ins == 0:
        return read_start_step_end(l_line)
    if ins == 6:
        return read_d1a(l_line)
    raise FullProfFormatError(
        f"The data format Ins = {ins:} of FullProf is not read; convert the data "
        "to Ins = 10 (2theta, intensity, sigma) or Ins = 0.")


def take_label(text: str, l_used: List[str]) -> str:
    """A data name of rcif from the name of a FullProf phase."""
    words = text.split()
    label = re.sub(r"\W", "_", words[0] if len(words) > 0 else "phase").lower()
    # reserved words of CIF ("data_pd" would start a data block)
    label = re.sub(r"^(data|loop|save|global|stop)_", "", label).strip("_")
    if label == "" or label[0].isdigit():
        label = f"phase_{label:}".strip("_")
    name, index = label, 1
    while name in l_used:
        index += 1
        name = f"{label:}_{index:}"
    return name


class FullProfPhase:
    """One phase of a .pcr file (only what the conversion uses)."""

    def __init__(self):
        self.label = ""
        self.space_group = ""
        self.atoms = []  # dicts
        self.scale, self.scale_code = 1., 0.
        self.resolution, self.resolution_codes = [0.] * 5, [0.] * 5
        self.igsize, self.igsize_code = 0., 0.
        self.cell, self.cell_codes = [], []
        self.pref, self.pref_codes = [0., 0.], [0., 0.]
        self.pref_axis = [0., 0., 1.]
        self.asymmetry, self.asymmetry_codes = [0.] * 4, [0.] * 4
        self.axial, self.axial_codes = [0., 0.], [0., 0.]


def read_fullprof_phase(lines: PcrLines, nor: int, l_used: List[str]) -> FullProfPhase:  # noqa: C901, PLR0912, PLR0915
    phase = FullProfPhase()
    phase.label = take_label(lines.next_text(), l_used)
    l_value = lines.next_numbers(11)
    nat, dis, ang = int(l_value[0]), int(l_value[1]), int(l_value[2])
    phase.pref_axis = l_value[3:6]
    jbt, irf, isy, i_str, furth = [int(value) for value in l_value[6:11]]
    nvk = int(l_value[12]) if len(l_value) > 12 else 0
    npr = int(l_value[13]) if len(l_value) > 13 else 7
    more = int(l_value[14]) if len(l_value) > 14 else 0
    if jbt != 0:
        raise FullProfFormatError(
            f"Phase '{phase.label:}' has Jbt = {jbt:}: only nuclear phases "
            "with atoms (Jbt = 0) are converted.")
    for name, value in (("Irf", irf), ("Isy", isy), ("Str", i_str),
                        ("Furth", furth), ("Nvk", nvk), ("Dis", dis), ("Ang", ang)):
        if value != 0:
            raise FullProfFormatError(
                f"Phase '{phase.label:}' has {name:} = {value:}: this option of "
                "FullProf is not converted.")
    if npr not in (0, 7):
        raise FullProfFormatError(
            f"Phase '{phase.label:}' has the profile Npr = {npr:}: only the "
            "Thompson-Cox-Hastings pseudo-Voigt (Npr = 7) is converted.")
    if more == 1:
        l_more = lines.next_integers(6)
        jdi, hel, sol, mom, ter = l_more[1:6]
        for name, value in (("Hel", hel), ("Sol", sol), ("Mom", mom), ("Ter", ter)):
            if value != 0:
                raise FullProfFormatError(
                    f"Phase '{phase.label:}' has {name:} = {value:}: this option "
                    "of FullProf is not converted.")
        if jdi in (1, 2, 3, 4):
            line = lines.next_text()  # distances, angles, BVS
            if "BVS" in line.upper():
                n_cation, n_anion = lines.next_integers(2)[:2]
                lines.next_text()  # names of the cations
                if n_anion > 0:
                    lines.next_text()  # names of the anions
    phase.space_group = lines.next_text().split("<--")[0].strip()

    for _ in range(nat):
        l_word = lines.next_text().split("#")[0].split()
        if len(l_word) < 10:
            raise FullProfFormatError(f"Wrong line of an atom: '{' '.join(l_word):}'.")
        atom = {
            "label": l_word[0], "type": take_type_symbol(l_word[1]),
            "xyz": [float(word) for word in l_word[2:5]],
            "b_iso": float(l_word[5]), "occ": float(l_word[6]),
            "n_t": int(l_word[9])}
        atom["codes"] = (lines.next_numbers(5)[:5])
        if atom["n_t"] == 2:
            atom["beta"] = lines.next_numbers(6)[:6]
            atom["beta_codes"] = lines.next_numbers(6)[:6]
        elif atom["n_t"] != 0:
            raise FullProfFormatError(
                f"Atom '{atom['label']:}' has N_t = {atom['n_t']:}: only isotropic "
                "(0) and anisotropic (2) atoms are converted.")
        phase.atoms.append(atom)

    l_value = lines.next_numbers(1)
    l_code = lines.next_numbers(1)
    phase.scale, phase.scale_code = l_value[0], l_code[0]
    if any(abs(value) > 1e-8 for value in l_value[1:6]):
        warn_conversion(
            f"Phase '{phase.label:}': Shape1, Bov and the strain parameters of "
            "FullProf are not converted.")
    l_value = lines.next_numbers(5)
    l_code = lines.next_numbers(5)
    phase.resolution, phase.resolution_codes = l_value[:5], l_code[:5]
    if len(l_value) > 5:
        phase.igsize = l_value[5]
        phase.igsize_code = l_code[5] if len(l_code) > 5 else 0.
    if len(l_value) > 6 and abs(l_value[6]) > 1e-8:
        warn_conversion(f"Phase '{phase.label:}': LorSiz of FullProf is not converted.")
    phase.cell = lines.next_numbers(6)[:6]
    phase.cell_codes = lines.next_numbers(6)[:6]
    l_value = lines.next_numbers(8)
    l_code = lines.next_numbers(8)
    phase.pref, phase.pref_codes = l_value[0:2], l_code[0:2]
    phase.asymmetry, phase.asymmetry_codes = l_value[2:6], l_code[2:6]
    phase.axial, phase.axial_codes = l_value[6:8], l_code[6:8]
    flag_pref = (nor == 1 and abs(phase.pref[0] - 1.) > 1e-8) or (
        nor == 0 and abs(phase.pref[0]) > 1e-8)
    if flag_pref:
        # FullProf (Nor = 1) uses the angle between the scattering vector and
        # the axis, CrysPy (texture) its complement: the functions differ
        warn_conversion(
            f"Phase '{phase.label:}': the preferred orientation of FullProf is not "
            "converted (the March-Dollase function of CrysPy takes sin^2 of the angle "
            "to the texture axis where FullProf takes cos^2).")
    return phase


class FullProfPattern:
    """What a .pcr file says about one powder pattern."""

    def __init__(self):
        self.title = ""
        self.job = 1
        self.nba = 0
        self.ins = 10
        self.nor = 0
        self.wavelength = 1.
        self.cthm = 0.
        self.bkpos = 0.
        self.asylim = 180.
        self.ttheta_min, self.ttheta_step, self.ttheta_max = 0., 0., 0.
        self.background_points = []  # (2theta, intensity, code)
        self.background_polynomial, self.background_codes = [], []
        self.excluded = []
        self.zero, self.zero_code = 0., 0.
        self.sysin, self.sysin_code = 0., 0.
        self.file_data = ""
        self.phases = []


def read_fullprof_pcr(file_pcr: str) -> FullProfPattern:  # noqa: C901, PLR0912, PLR0915
    """Read a .pcr file of FullProf (one pattern, constant wavelength powder).

    Raise FullProfFormatError for the options the conversion does not
    take, warn about the values it cannot keep.
    """
    with open(file_pcr, "r", errors="replace") as fid:
        text = fid.read()
    pattern = FullProfPattern()
    match = re.search(r"DAT-file:\s*([^,\s]+)", text)
    if match is not None:
        pattern.file_data = match.group(1)
    lines = PcrLines(text)

    line = lines.next_text()
    if line.upper().startswith("NPATT"):
        raise FullProfFormatError("Multi-pattern .pcr files of FullProf are not converted.")
    pattern.title = line[4:].strip() if line.upper().startswith("COMM") else line.strip()
    if not line.upper().startswith("COMM"):
        raise FullProfFormatError(
            "The .pcr file should start with the line 'COMM title' "
            "(the old fixed format of FullProf is not converted).")

    l_job = lines.next_integers(15)
    (job, npr, nph, nba, nex, nsc, nor, _, _, _, _, _, _, nre, cry) = l_job[:15]
    pattern.job, pattern.nba, pattern.nor = job, nba, nor
    if job not in D_JOB_RADIATION:
        raise FullProfFormatError(
            f"Job = {job:} of FullProf is not converted: only the powder at a constant "
            "wavelength (Job = 0, 1, 2, 3) is (TOF is not yet).")
    if cry != 0:
        raise FullProfFormatError("Single crystal data (Cry != 0) of FullProf are not converted.")
    if nsc != 0:
        raise FullProfFormatError(
            "User scattering factors (Nsc != 0) of FullProf are not converted.")
    if nba not in (0, -1) and nba < 2:
        raise FullProfFormatError(
            f"The background Nba = {nba:} of FullProf is not converted "
            "(points Nba >= 2 or the polynomial Nba = 0 are).")
    l_ipr = lines.next_integers(11)
    pattern.ins = l_ipr[10]

    l_value = lines.next_numbers(6)
    lambda_1, lambda_2, ratio = l_value[0:3]
    pattern.wavelength = lambda_1
    pattern.bkpos, pattern.cthm = l_value[3], l_value[5]
    if len(l_value) > 6 and abs(l_value[6]) > 1e-8:
        warn_conversion("The absorption muR of FullProf is not converted.")
    if len(l_value) > 7:
        pattern.asylim = l_value[7]
    if abs(ratio) > 1e-8 and abs(lambda_2 - lambda_1) > 1e-6:
        warn_conversion(
            f"The second wavelength {lambda_2:} of FullProf (Ratio = {ratio:}) is "
            "not converted: CrysPy calculates the pattern for lambda1 only.")
    if len(l_value) > 8 and abs(l_value[8]) > 1e-8:
        warn_conversion("Rpolarz of FullProf is not converted.")

    l_value = lines.next_numbers(9)
    pattern.ttheta_min, pattern.ttheta_step, pattern.ttheta_max = l_value[6:9]

    if nba >= 2:
        for _ in range(nba):
            l_value = lines.next_numbers(2)
            pattern.background_points.append(
                (l_value[0], l_value[1], l_value[2] if len(l_value) > 2 else 0.))
    for _ in range(nex):
        pattern.excluded.append(lines.next_numbers(2)[:2])
    lines.next_numbers(1)  # number of refined parameters

    l_value = lines.next_numbers(8)
    pattern.zero, pattern.zero_code = l_value[0], l_value[1]
    if abs(l_value[2]) > 1e-8:
        warn_conversion(
            "SyCos of FullProf (shift SyCos*cos(theta)) is not converted: "
            "offset_SyCos of CrysPy multiplies cos(2theta).")
    pattern.sysin, pattern.sysin_code = l_value[4], l_value[5]
    if len(l_value) > 8 and int(l_value[8]) == 1:
        lines.next_text()  # microabsorption
        warn_conversion("The microabsorption (MORE = 1) of FullProf is not converted.")
    if nba == 0:
        pattern.background_polynomial = lines.next_numbers(6)[:6]
        pattern.background_codes = lines.next_numbers(6)[:6]
    elif nba == -1:
        for _ in range(6):
            lines.next_text()  # Debye-like and polynomial terms, taken from the .prf file

    l_used = []
    for _ in range(nph):
        phase = read_fullprof_phase(lines, nor, l_used)
        l_used.append(phase.label)
        pattern.phases.append(phase)
    if nre > 0:
        logger.info("The limits of %d parameters in the .pcr file are not converted.", nre)
    return pattern


def form_crystal_rcif(phase: FullProfPhase, l_occupancy: Optional[List[float]] = None) -> str:
    """rcif text of the crystal of a phase (occupancies of CIF, 1 by default)."""
    if l_occupancy is None:
        l_occupancy = [1.] * len(phase.atoms)
    ls_out = [f"data_{phase.label:}", ""]
    for name, value, code in zip(
            ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"),
            phase.cell, phase.cell_codes):
        ls_out.append(f"_cell_{name:} {form_value(value, code):}")
    ls_out.extend(["", f'_space_group_name_H-M_alt "{phase.space_group:}"', ""])
    ls_out.extend(["loop_", "_atom_site_label", "_atom_site_type_symbol",
                   "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
                   "_atom_site_occupancy", "_atom_site_adp_type",
                   "_atom_site_B_iso_or_equiv"])
    for atom, occupancy in zip(phase.atoms, l_occupancy):
        codes = atom["codes"]
        adp_type = "Uani" if atom["n_t"] == 2 else "Biso"
        ls_out.append(" ".join([
            atom["label"], atom["type"],
            *[form_value(value, code) for value, code in zip(atom["xyz"], codes[:3])],
            form_value(occupancy, codes[4]), adp_type, form_value(atom["b_iso"], codes[3])]))
    l_aniso = [atom for atom in phase.atoms if atom["n_t"] == 2]
    if len(l_aniso) > 0:
        ls_out.extend(["", "loop_", "_atom_site_aniso_label"] + [
            f"_atom_site_aniso_U_{ij:}" for ij in ("11", "22", "33", "12", "13", "23")])
        a_r = calc_reciprocal_lengths(phase.cell)
        for atom in l_aniso:
            l_b = []
            for (i, j), beta, code in zip(
                    ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)),
                    atom["beta"], atom["beta_codes"]):
                # beta_ij = 2 pi^2 U_ij a*_i a*_j
                l_b.append(form_value(beta / (2. * numpy.pi**2 * a_r[i] * a_r[j]), code))
            ls_out.append(" ".join([atom["label"]] + l_b))
            if abs(atom["b_iso"]) > 1e-8:
                warn_conversion(
                    f"Atom '{atom['label']:}' has betas and Biso = {atom['b_iso']:}; "
                    "only the betas are converted.")
    return "\n".join(ls_out)


def calc_reciprocal_lengths(cell: List[float]) -> numpy.ndarray:
    """a*, b*, c* (1/Å) of a cell (a, b, c in Å, angles in degrees)."""
    length = numpy.array(cell[:3], dtype=float)
    cos_angle = numpy.cos(numpy.radians(numpy.array(cell[3:6], dtype=float)))
    sin_angle = numpy.sqrt(1. - numpy.square(cos_angle))
    volume = numpy.prod(length) * math.sqrt(
        1. - numpy.square(cos_angle).sum() + 2. * numpy.prod(cos_angle))
    return numpy.array([
        length[1] * length[2] * sin_angle[0], length[0] * length[2] * sin_angle[1],
        length[0] * length[1] * sin_angle[2]]) / volume


def calc_occupancies(phase: FullProfPhase) -> List[float]:
    """Occupancies of CIF: FullProf gives Occ = occupancy * m_site /
    m_general (FullProf manual)."""
    crystal = str_to_globaln(form_crystal_rcif(phase)).items[0]
    dict_crystal = crystal.get_dictionary()
    n_general = dict_crystal["full_symm_elems"].shape[1]
    return [atom["occ"] * n_general / multiplicity
            for atom, multiplicity in zip(phase.atoms, dict_crystal["atom_multiplicity"])]


def form_pd_rcif(pattern: FullProfPattern, name: str,  # noqa: C901, PLR0912, PLR0915
                 ttheta: numpy.ndarray, intensity: numpy.ndarray, sigma: numpy.ndarray,
                 background: Optional[Tuple[numpy.ndarray, numpy.ndarray]]) -> str:
    radiation = D_JOB_RADIATION[pattern.job]
    phase_1 = pattern.phases[0]
    ls_out = [f"data_{name:}", ""]
    ls_out.append(f"_setup_wavelength {form_value(pattern.wavelength, 0., 5):}")
    ls_out.append("_setup_field 0.0")
    ls_out.append(f"_setup_offset_2theta {form_value(pattern.zero, pattern.zero_code, 5):}")
    if abs(pattern.sysin) > 1e-8 or is_refined(pattern.sysin_code):
        ls_out.append(f"_setup_offset_SySin {form_value(pattern.sysin, pattern.sysin_code, 5):}")
    ls_out.append(f"_setup_radiation {radiation:}")
    ls_out.append(f"_setup_K {0.5 if radiation == 'X-rays' else 0.:.1f}")
    ls_out.append(f"_setup_cthm {pattern.cthm:.5f}")
    ls_out.append("")
    for name_uvw, value, code in zip("UVWXY", phase_1.resolution, phase_1.resolution_codes):
        ls_out.append(f"_pd_instr_resolution_{name_uvw:} {form_value(value, code):}")
    ls_out.append("")

    flag_asymmetry = any(abs(value) > 1e-8 for value in phase_1.asymmetry)
    if flag_asymmetry:
        for index, (value, code) in enumerate(zip(phase_1.asymmetry, phase_1.asymmetry_codes)):
            ls_out.append(f"_pd_instr_reflex_asymmetry_p{index + 1:} {form_value(value, code):}")
        if pattern.asylim < pattern.ttheta_max:
            warn_conversion(
                f"AsyLim = {pattern.asylim:} of FullProf is not converted: CrysPy applies "
                "the empirical asymmetry at all angles.")
        ls_out.append("")
    if any(abs(value) > 1e-8 for value in phase_1.axial):
        for name_sd, value, code in zip(("s_l", "d_l"), phase_1.axial, phase_1.axial_codes):
            ls_out.append(f"_pd_instr_axial_divergence_{name_sd:} {form_value(value, code):}")
        if pattern.asylim < 180.:
            ls_out.append(f"_pd_instr_axial_divergence_ttheta_max {pattern.asylim:.2f}")
        ls_out.append("")
    for phase in pattern.phases[1:]:
        if phase.asymmetry != phase_1.asymmetry or phase.axial != phase_1.axial:
            warn_conversion(
                f"Phase '{phase.label:}' has another asymmetry than phase "
                f"'{phase_1.label:}': CrysPy keeps one for the pattern, the first one.")

    ls_out.append(f"_range_2theta_min {pattern.ttheta_min:.4f}")
    ls_out.append(f"_range_2theta_max {pattern.ttheta_max:.4f}")
    ls_out.append("")

    ls_out.extend(["loop_", "_phase_label", "_phase_scale", "_phase_igsize",
                   "_phase_U", "_phase_V", "_phase_W", "_phase_X", "_phase_Y"])
    for phase in pattern.phases:
        # the resolution of CrysPy is the one of the instrument plus the one of the phase
        l_res = [
            form_value(value - value_1, code if phase is not phase_1 else 0.)
            for value, value_1, code in zip(
                phase.resolution, phase_1.resolution, phase.resolution_codes)]
        ls_out.append(" ".join([
            phase.label, form_value(phase.scale * take_scale_factor(radiation),
                                    phase.scale_code, 8),
            form_value(phase.igsize, phase.igsize_code)] + l_res))
    ls_out.append("")

    if len(pattern.excluded) > 0:
        ls_out.extend(["loop_", "_exclude_id", "_exclude_2theta_min", "_exclude_2theta_max"])
        for index, (low, high) in enumerate(pattern.excluded):
            ls_out.append(f"{index + 1:} {low:.4f} {high:.4f}")
        ls_out.append("")

    ls_out.extend(["loop_", "_pd_background_2theta", "_pd_background_intensity"])
    background_ttheta, background_intensity, background_refined = background
    for value_ttheta, value, flag in zip(background_ttheta, background_intensity, background_refined):
        ls_out.append(f"{value_ttheta:.4f} {form_value(value, 1. if flag else 0., 4):}")
    ls_out.append("")

    ls_out.extend(["loop_", "_pd_meas_2theta", "_pd_meas_intensity", "_pd_meas_intensity_sigma"])
    for value_ttheta, value, value_sigma in zip(ttheta, intensity, sigma):
        ls_out.append(f"{value_ttheta:.4f} {value:.4f} {value_sigma:.4f}")
    return "\n".join(ls_out)


def take_background(pattern: FullProfPattern, file_pcr: str
                    ) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Points of the background: the points of FullProf, the polynomial
    sampled at N_BACKGROUND_POINTS points, or (Nba = -1) the background of
    the .prf file next to the .pcr file."""
    if pattern.nba >= 2:
        array = numpy.array(pattern.background_points, dtype=float)
        return array[:, 0], array[:, 1], numpy.abs(array[:, 2]) > 1e-6
    ttheta = numpy.linspace(pattern.ttheta_min, pattern.ttheta_max, N_BACKGROUND_POINTS)
    if pattern.nba == 0:
        if pattern.bkpos <= 0.:
            raise FullProfFormatError("Bkpos of the polynomial background should be positive.")
        intensity = calc_fullprof_polynomial_background(
            ttheta, pattern.background_polynomial, pattern.bkpos)
        flag_refined = any(is_refined(code) for code in pattern.background_codes)
        return ttheta, intensity, numpy.full(ttheta.shape, flag_refined)
    file_prf = os.path.splitext(file_pcr)[0] + ".prf"
    if not os.path.isfile(file_prf):
        raise FullProfFormatError(
            f"The background Nba = {pattern.nba:} of FullProf is taken from the .prf "
            f"file, {file_prf:} is not found (write it with Prf = 3).")
    prf_ttheta, _, prf_background = read_prf_columns(file_prf)
    warn_conversion(
        f"The background Nba = {pattern.nba:} is described by {N_BACKGROUND_POINTS:} "
        "points of the background of the .prf file.")
    return ttheta, numpy.interp(ttheta, prf_ttheta, prf_background), numpy.zeros(ttheta.shape, bool)


def find_data_file(pattern: FullProfPattern, file_pcr: str) -> str:
    directory = os.path.dirname(os.path.abspath(file_pcr))
    l_name = []
    if pattern.file_data != "":
        l_name.append(pattern.file_data)
    l_name.append(os.path.splitext(os.path.basename(file_pcr))[0])
    for name in l_name:
        for extension in ("", ".dat", ".xys", ".xye", ".DAT"):
            file_data = os.path.join(directory, name + extension)
            if os.path.isfile(file_data):
                return file_data
    raise FullProfFormatError(
        f"The data file of {file_pcr:} is not found; give it with file_data.")


def warn_shared_codes(pattern: FullProfPattern):
    """FullProf constrains the parameters with the same code number; CrysPy
    constrains special positions itself, other links are lost."""
    d_owner = {}
    for phase in pattern.phases:
        for atom in phase.atoms:
            for code in atom["codes"] + atom.get("beta_codes", []):
                if is_refined(code):
                    d_owner.setdefault(int(abs(code)) // 10, set()).add(
                        f"{phase.label:}.{atom['label']:}")
        for code in [phase.scale_code] + phase.cell_codes + phase.resolution_codes:
            if is_refined(code):
                d_owner.setdefault(int(abs(code)) // 10, set()).add(phase.label)
    l_shared = sorted(
        f"{number:}: {', '.join(sorted(owner)):}" for number, owner in d_owner.items()
        if len(owner) > 1)
    if len(l_shared) > 0:
        warn_conversion(
            "These FullProf parameters share a code (a constraint), CrysPy refines them "
            "independently; add a constraint if they should stay equal: "
            + "; ".join(l_shared) + ".")


def load_fullprof(file_pcr: str, file_data: str = "") -> GlobalN:
    """Read a FullProf .pcr file and its data file into a CrysPy object.

    The crystals of the nuclear phases and the powder experiment at a
    constant wavelength (neutrons or X-rays) are converted: cell, space
    group, atoms with occupancies and isotropic or anisotropic
    displacements, scales, U V W X Y, zero shift, background, excluded
    regions, asymmetry and axial divergence. Refined parameters of FullProf (code != 0) are
    refined in CrysPy too. For a simulation (Job = 2, 3) the pattern has
    zero intensities.

    Menu: -

    Parameters
    ----------
    file_pcr : str
        The .pcr file of FullProf.
    file_data : str
        The data file; by default the DAT-file named in the .pcr file or
        the file with the same name and the extension .dat.

    Returns
    -------
    GlobalN
        The object with the crystals and the powder experiment.
    """
    pattern = read_fullprof_pcr(file_pcr)
    if pattern.job in JOBS_SIMULATION:
        ttheta = numpy.arange(pattern.ttheta_min, pattern.ttheta_max + 0.5 * pattern.ttheta_step,
                              pattern.ttheta_step)
        intensity, sigma = numpy.zeros_like(ttheta), numpy.ones_like(ttheta)
    else:
        if file_data == "":
            file_data = find_data_file(pattern, file_pcr)
        ttheta, intensity, sigma = read_fullprof_data(file_data, pattern.ins)
    background = take_background(pattern, file_pcr)
    warn_shared_codes(pattern)

    name_pd = take_label(os.path.splitext(os.path.basename(file_pcr))[0],
                         [phase.label for phase in pattern.phases])
    ls_text = ["global_", ""]
    for phase in pattern.phases:
        ls_text.append(form_crystal_rcif(phase, calc_occupancies(phase)))
        ls_text.append("")
    ls_text.append(form_pd_rcif(pattern, name_pd, ttheta, intensity, sigma, background))
    obj = str_to_globaln("\n".join(ls_text))
    logger.info("FullProf file %s is converted: %d phases, %d points.",
                file_pcr, len(pattern.phases), ttheta.size)
    return obj


def add_fullprof_files(rcif_object: GlobalN, file_pcr: str, file_data: str = ""):
    """Add the crystals and the powder experiment of a FullProf .pcr file
    and its data file to the object, to check a FullProf refinement in
    CrysPy.

    The nuclear phases and the powder pattern at a constant wavelength
    (neutrons or X-rays) are converted; options of FullProf that CrysPy
    does not calculate (magnetic phases, TOF, single crystals, other
    profiles) are refused with a message, values that cannot be kept give
    a warning.

    Menu: FullProf

    Parameters
    ----------
    rcif_object : GlobalN
        The object to add to.
    file_pcr : str
        The .pcr file of FullProf.
    file_data : str
        The data file; by default the DAT-file named in the .pcr file.

    Returns
    -------
    None
    """
    obj = load_fullprof(file_pcr, file_data=file_data)
    l_used = [str(getattr(item, "data_name", "")).lower() for item in rcif_object.items]
    l_new = [item for item in obj.items if str(getattr(item, "data_name", "")).lower()
             not in l_used]
    if len(l_new) < len(obj.items):
        raise ValueError(
            "The object already has data blocks with the names of the FullProf "
            "phases or pattern; start from an empty object.")
    rcif_object.add_items(l_new)


class CodeCounter:
    """Codes of FullProf for the refined parameters: 10 n + 1, the same n
    for parameters that CrysPy keeps equal (the resolution of the
    instrument for every phase)."""

    def __init__(self):
        self.number = 0
        self.d_shared = {}

    def take(self, flag_refined: bool, shared_key: Optional[str] = None) -> float:
        if not flag_refined:
            return 0.
        if shared_key is not None and shared_key in self.d_shared:
            return self.d_shared[shared_key]
        self.number += 1
        code = 10. * self.number + 1.
        if shared_key is not None:
            self.d_shared[shared_key] = code
        return code


def take_refined(item, name: str) -> bool:
    return bool(getattr(item, f"{name:}_refinement", False)) and not bool(
        getattr(item, f"{name:}_constraint", False))


def take_value(item, name: str, default: float = 0.) -> float:
    value = getattr(item, name, None) if item is not None else None
    return default if value is None else float(value)


def take_optional(container, name: str):
    return getattr(container, name) if container.is_attribute(name) else None


def form_numbers(values, width: int = 10, digits: int = 5) -> str:
    return "".join(f"{value:{width:}.{digits:}f}" for value in values)


def form_atom_lines(crystal, radiation: str, counter: CodeCounter) -> Tuple[List[str], bool]:  # noqa: PLR0915
    """Lines of the atoms of a crystal in a .pcr file and whether some
    refined atom is on a special position."""
    dict_crystal = crystal.get_dictionary()
    n_general = dict_crystal["full_symm_elems"].shape[1]
    atom_multiplicity = dict_crystal["atom_multiplicity"]
    cell = crystal.cell
    a_r = calc_reciprocal_lengths([cell.length_a, cell.length_b, cell.length_c,
                                   cell.angle_alpha, cell.angle_beta, cell.angle_gamma])
    atom_site_aniso = take_optional(crystal, "atom_site_aniso")
    d_aniso = {} if atom_site_aniso is None else {
        item.label: item for item in atom_site_aniso.items}
    ls_out = []
    flag_special_refined = False
    for atom, multiplicity in zip(crystal.atom_site.items, atom_multiplicity):
        adp_type = str(getattr(atom, "adp_type", "Biso") or "Biso")
        if adp_type in ("Uiso", "Biso"):
            if adp_type == "Uiso":
                b_iso = 8. * numpy.pi**2 * take_value(atom, "u_iso_or_equiv")
                code_b = counter.take(take_refined(atom, "u_iso_or_equiv"))
            else:
                b_iso = take_value(atom, "b_iso_or_equiv")
                code_b = counter.take(take_refined(atom, "b_iso_or_equiv"))
            n_t = 0
        else:
            b_iso, code_b, n_t = 0., 0., 2
        occ = take_value(atom, "occupancy", 1.) * multiplicity / n_general
        l_code = [counter.take(take_refined(atom, name))
                  for name in ("fract_x", "fract_y", "fract_z")]
        code_occ = counter.take(take_refined(atom, "occupancy"))
        if multiplicity < n_general and any(is_refined(code) for code in l_code):
            flag_special_refined = True
        ls_out.append(
            f"{atom.label:<6} {form_fullprof_type(atom.type_symbol, radiation):<6}"
            + form_numbers([atom.fract_x, atom.fract_y, atom.fract_z, b_iso, occ], 9, 5)
            + f"   0   0 {n_t:3}    0")
        ls_out.append(" " * 14 + form_numbers(l_code + [code_b, code_occ], 9, 2))
        if n_t == 2:
            aniso = d_aniso.get(atom.label)
            if aniso is None:
                raise ValueError(f"Atom '{atom.label:}' is anisotropic without atom_site_aniso.")
            l_beta, l_code_beta = [], []
            for (i, j), ij in zip(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)),
                                  ("11", "22", "33", "12", "13", "23")):
                if adp_type == "Bani":
                    value = take_value(aniso, f"b_{ij:}") / 4.
                    flag = take_refined(aniso, f"b_{ij:}")
                else:
                    value = 2. * numpy.pi**2 * take_value(aniso, f"u_{ij:}")
                    flag = take_refined(aniso, f"u_{ij:}")
                l_beta.append(value * a_r[i] * a_r[j])
                code = counter.take(flag)
                l_code_beta.append(code)
                if multiplicity < n_general and is_refined(code):
                    flag_special_refined = True
            ls_out.append("    " + form_numbers(l_beta, 11, 7))
            ls_out.append("    " + form_numbers(l_code_beta, 11, 2))
    return ls_out, flag_special_refined


def form_fullprof_pcr(pd, l_crystal: list, file_data_name: str,  # noqa: C901, PLR0912, PLR0915
                      ttheta: numpy.ndarray) -> str:
    """Text of a .pcr file for a powder experiment and its crystals."""
    setup = pd.setup
    radiation = str(getattr(setup, "radiation", "neutrons") or "neutrons")
    job = 0 if radiation == "X-rays" else 1
    counter = CodeCounter()

    background = pd.pd_background
    l_bkg = []
    for item in background.items:
        l_bkg.append((item.ttheta, item.intensity,
                      counter.take(take_refined(item, "intensity"))))
    if len(l_bkg) < 2:
        raise ValueError("FullProf needs at least two background points.")
    exclude = take_optional(pd, "exclude")
    l_excluded = [] if exclude is None else [
        (item.ttheta_min, item.ttheta_max) for item in exclude.items
        if getattr(item, "ttheta_min", None) is not None]

    code_zero = counter.take(take_refined(setup, "offset_ttheta"))
    sysin = take_value(setup, "offset_sysin")
    code_sysin = counter.take(take_refined(setup, "offset_sysin"))
    if abs(take_value(setup, "offset_sycos")) > 1e-8:
        warn_conversion(
            "offset_SyCos of CrysPy (shift with cos(2theta)) is not written: "
            "SyCos of FullProf multiplies cos(theta).")

    resolution = pd.pd_instr_resolution
    asymmetry = take_optional(pd, "pd_instr_reflex_asymmetry")
    axial = take_optional(pd, "pd_instr_axial_divergence")
    asylim = 180.
    if axial is not None and getattr(axial, "ttheta_max", None) is not None:
        asylim = float(axial.ttheta_max)
        if asymmetry is not None and asylim < 180.:
            warn_conversion(
                "AsyLim of FullProf is the limit of the axial divergence "
                f"({asylim:}); it limits the empirical asymmetry too.")
    for name in ("texture", "texture_sh", "magnetic_broadening", "stacking_fault"):
        if pd.is_attribute(name):
            warn_conversion(f"The loop {name:} of CrysPy is not written to FullProf.")

    d_crystal = {str(crystal.data_name).lower(): crystal for crystal in l_crystal}
    l_phase = []
    flag_special = False
    for phase in pd.phase.items:
        crystal = d_crystal.get(str(phase.label).lower())
        if crystal is None:
            raise ValueError(f"The crystal of the phase '{phase.label:}' is not found.")
        for name in ("atom_site_susceptibility", "atom_site_moment"):
            if crystal.is_attribute(name):
                warn_conversion(
                    f"The magnetic atoms ({name:}) of the crystal '{crystal.data_name:}' are not "
                    "written to FullProf: only the nuclear structure is compared.")
        l_phase.append((phase, crystal))

    ls_phase = []
    for index, (phase, crystal) in enumerate(l_phase):
        l_atom_line, flag = form_atom_lines(crystal, radiation, counter)
        flag_special = flag_special or flag
        nat = len(crystal.atom_site.items)
        ls = ["!-------------------------------------------------------------------------------",
              f"!  Data for PHASE number: {index + 1:3}",
              "!-------------------------------------------------------------------------------",
              f"{crystal.data_name:}",
              "!",
              "!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More",
              f"{nat:4}   0   0 0.0 0.0 1.0   0   0   0   0   0         0.000   0   7   0",
              "!",
              f"{crystal.space_group.name_hm_alt:<24} <--Space group symbol",
              "!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes"]
        ls.extend(l_atom_line)
        scale = take_value(phase, "scale", 1.) / take_scale_factor(radiation)
        ls.extend([
            f"!-------> Profile Parameters for Pattern #  1  ----> Phase # {index + 1:3}",
            "!  Scale        Shape1      Bov      Str1      Str2      Str3   Strain-Model",
            f"  {scale:.7e}   0.00000   0.00000   0.00000   0.00000   0.00000       0",
            f"  {counter.take(take_refined(phase, 'scale')):10.2f}     0.000     0.000     0.000     0.000     0.000"])
        l_uvwxy, l_code_uvwxy = [], []
        for name in ("u", "v", "w", "x", "y"):
            value = take_value(resolution, name) + take_value(phase, name)
            l_uvwxy.append(value)
            # the instrument is one for all phases: the same code number
            flag_instr = take_refined(resolution, name)
            flag_phase = take_refined(phase, name)
            if flag_instr:
                l_code_uvwxy.append(counter.take(True, f"resolution_{name:}"))
            else:
                l_code_uvwxy.append(counter.take(flag_phase))
        igsize = take_value(phase, "igsize")
        ls.extend([
            "!       U         V          W           X          Y        GauSiz   LorSiz Size-Model",
            form_numbers(l_uvwxy + [igsize, 0.], 11, 6) + "    0",
            form_numbers(l_code_uvwxy + [counter.take(take_refined(phase, "igsize")), 0.], 11, 3)])
        cell = crystal.cell
        l_cell_name = ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")
        ls.extend([
            "!     a          b         c        alpha      beta       gamma      #Cell Info",
            form_numbers([getattr(cell, name) for name in l_cell_name], 11, 6),
            form_numbers([counter.take(take_refined(cell, name)) for name in l_cell_name], 11, 5)])
        l_asy, l_code_asy = [0.] * 4, [0.] * 4
        if asymmetry is not None:
            l_asy = [take_value(asymmetry, f"p{i:}") for i in range(1, 5)]
            l_code_asy = [counter.take(take_refined(asymmetry, f"p{i:}"), f"asymmetry_p{i:}")
                          for i in range(1, 5)]
        l_axial, l_code_axial = [0., 0.], [0., 0.]
        if axial is not None:
            l_axial = [take_value(axial, "s_l"), take_value(axial, "d_l")]
            l_code_axial = [counter.take(take_refined(axial, name), f"axial_{name:}")
                            for name in ("s_l", "d_l")]
        ls.extend([
            "!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4      S_L      D_L",
            form_numbers([0., 0.] + l_asy + l_axial, 9, 5),
            form_numbers([0., 0.] + l_code_asy + l_code_axial, 9, 2)])
        ls_phase.extend(ls)
    if flag_special:
        warn_conversion(
            "Refined coordinates or betas of atoms on special positions are written "
            "with independent codes: check the codes of FullProf for these atoms.")

    step = float(numpy.mean(numpy.diff(ttheta))) if ttheta.size > 1 else 0.
    if ttheta.size > 2 and numpy.ptp(numpy.diff(ttheta)) > 1e-3 * abs(step):
        warn_conversion("The 2theta step of the data is not constant; FullProf expects a constant step.")
    l_range = take_optional(pd, "range")
    ttheta_min = take_value(l_range, "ttheta_min", float(ttheta.min()))
    ttheta_max = take_value(l_range, "ttheta_max", float(ttheta.max()))

    nba = len(l_bkg)
    ls_out = [
        f"COMM {pd.data_name:} (written by CrysPy)",
        f"! Files => DAT-file: {file_data_name:},  PCR-file: {file_data_name:}",
        "!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut",
        f"{job:4}   7{len(l_phase):4}{nba:4}{len(l_excluded):4}   0   0   0   0   0   0   0   0   0   0   0   0   0   1",
        "!",
        "!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 NLI Prf Ins Rpa Sym Hkl Fou Sho Ana",
        "   0   0   1   0   1   0   4   0   0   3  10   0   0   0   0   0   0",
        "!",
        "! Lambda1  Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz  2nd-muR -> Patt# 1",
        f" {setup.wavelength:.6f} {setup.wavelength:.6f}  0.00000   90.000  8.0000"
        f" {take_value(setup, 'cthm'):7.4f}  0.0000  {asylim:7.2f}    0.0000  0.0000",
        "!",
        "!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0",
        f" 10  0.10  1.00  1.00  1.00  1.00  {ttheta_min:10.4f} {step:10.6f} {ttheta_max:10.4f}   0.000   0.000",
        "!",
        "!2Theta/TOF/E(Kev)   Background  for Pattern#  1"]
    for value_ttheta, value, code in l_bkg:
        ls_out.append(f"{value_ttheta:15.4f}{value:15.4f}{code:10.2f}")
    if len(l_excluded) > 0:
        ls_out.append("! Excluded regions (LowT  HighT) for Pattern#  1")
        for low, high in l_excluded:
            ls_out.append(f"{low:12.2f}{high:12.2f}")
    ls_out.extend([
        "!",
        "NREFINED",  # known when all codes are given
        "!",
        "!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1",
        f" {take_value(setup, 'offset_ttheta'):8.5f} {code_zero:6.2f}  0.00000   0.00"
        f" {sysin:8.5f} {code_sysin:6.2f} 0.000000    0.00   0"])
    ls_out.extend(ls_phase)
    ls_out.extend(["!  2Th1/TOF1    2Th2/TOF2  Pattern to plot",
                   f"{ttheta_min:12.3f}{ttheta_max:12.3f}       1"])
    text = "\n".join(ls_out) + "\n"
    return text.replace("NREFINED", f"{counter.number:8}    !Number of refined parameters")


def save_fullprof_files(rcif_object: GlobalN, file_pcr: str):
    """Write the crystals and the powder experiments of the object as
    FullProf files (.pcr and the data file .dat, Ins = 10), to check the
    calculation in FullProf.

    One .pcr file is written for each powder experiment at a constant
    wavelength; with several experiments the name of the experiment is
    added to the file name. Refined parameters get codes of FullProf, the
    resolution of the instrument has the same code in every phase.
    Categories that FullProf does not calculate in the same way (texture,
    polarized neutrons, magnetic structures) are not written, with a
    warning.

    Menu: FullProf

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    file_pcr : str
        The .pcr file to write.

    Returns
    -------
    List[str]
        The names of the written .pcr files.
    """
    from cryspy.E_data_classes.cl_1_crystal import Crystal
    from cryspy.E_data_classes.cl_2_pd import Pd

    l_crystal = [item for item in rcif_object.items if isinstance(item, Crystal)]
    l_pd = [item for item in rcif_object.items if isinstance(item, Pd)]
    if len(l_pd) == 0:
        raise ValueError(
            "The object has no powder experiment at a constant wavelength (pd) to write "
            "to FullProf.")
    l_other = [type(item).__name__ for item in rcif_object.items
               if not isinstance(item, (Crystal, Pd)) and hasattr(item, "data_name")]
    if len(l_other) > 0:
        warn_conversion(f"These data blocks are not written to FullProf: {', '.join(l_other):}.")
    base, _ = os.path.splitext(file_pcr)
    l_file = []
    for pd in l_pd:
        file_base = base if len(l_pd) == 1 else f"{base:}_{pd.data_name:}"
        ttheta = numpy.array(pd.pd_meas.ttheta, dtype=float)
        if pd.pd_meas.is_attribute("intensity") and getattr(pd.pd_meas, "intensity", None) is not None:
            intensity = numpy.array(pd.pd_meas.intensity, dtype=float)
            sigma = numpy.array(pd.pd_meas.intensity_sigma, dtype=float)
        else:
            intensity = numpy.array(pd.pd_meas.intensity_plus, dtype=float) + numpy.array(
                pd.pd_meas.intensity_minus, dtype=float)
            sigma = numpy.sqrt(numpy.square(numpy.array(pd.pd_meas.intensity_plus_sigma, dtype=float))
                               + numpy.square(numpy.array(pd.pd_meas.intensity_minus_sigma, dtype=float)))
            warn_conversion(
                f"Experiment '{pd.data_name:}' is polarized: FullProf gets the sum of the "
                "intensities up and down, the magnetic part is not converted.")
        name_data = os.path.basename(file_base)
        text_pcr = form_fullprof_pcr(pd, l_crystal, name_data, ttheta)
        directory = os.path.dirname(os.path.abspath(file_base))
        os.makedirs(directory, exist_ok=True)
        with open(file_base + ".pcr", "w") as fid:
            fid.write(text_pcr)
        with open(file_base + ".dat", "w") as fid:
            fid.write(f"XYDATA {pd.data_name:} (written by CrysPy)\n")
            for value_ttheta, value, value_sigma in zip(ttheta, intensity, sigma):
                fid.write(f"{value_ttheta:12.4f}{value:14.4f}{value_sigma:12.4f}\n")
        l_file.append(file_base + ".pcr")
        logger.info("FullProf files %s.pcr and %s.dat are written.", file_base, file_base)
    return l_file
